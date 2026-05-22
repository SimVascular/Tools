"""Combined VTP -> capped MDL, driven by an existing mesh-surfaces/ partition.

Requires Python with the `vtk` package installed.

The script does not depend on the SimVascular Python API (`sv`) and can be
run in any Python environment with VTK available. Using SimVascular's bundled
Python may still be convenient because it already includes VTK.

MacBook example:
  /Applications/SimVascular.app/Contents/Resources/simvascular --python -- \
      vtp_to_capped_mdl.py --mesh-exterior mesh.exterior.vtp --mesh-surfaces-dir mesh-surfaces \
      --model-name model

Linux example:
  /usr/local/sv/2025-12-21/simvascular --python -- \
      vtp_to_capped_mdl.py --mesh-exterior mesh.exterior.vtp --mesh-surfaces-dir mesh-surfaces \
      --model-name model

Alternatively, if you have a Python environment with VTK set up, you can run the script directly:
  python vtp_to_capped_mdl.py \
    --mesh-exterior mesh.exterior.vtp --mesh-surfaces-dir mesh-surfaces --model-name model
"""
import argparse
import os
import sys
import xml.etree.ElementTree as ET
import vtk

FACE_ID_ARRAY_NAME = "ModelFaceID"
NODE_ID_ARRAY_NAME = "GlobalNodeID"
WALL_PREFIXES = ("wall_", "wall_blend_")


def load_vtp(filepath):
    if not os.path.exists(filepath):
        print(f"Error: File not found {filepath}")
        sys.exit(1)
    reader = vtk.vtkXMLPolyDataReader()
    reader.SetFileName(filepath)
    reader.Update()
    return reader.GetOutput()


def is_wall_file(fname):
    return fname.startswith(WALL_PREFIXES)


def clean_cap_name(filename):
    base = os.path.splitext(filename)[0]
    if base == "r_pa_x_2" or base == "inflow":
        return "inflow"
    if base.endswith("_x"):
        base = base[:-2]
    if "_" in base:
        parts = base.split("_", 1)
        base = parts[0] + parts[1]
    return base


def face_name_and_type(fname):
    """Cap files get a cleaned anatomical name; wall files keep their basename."""
    if is_wall_file(fname):
        return os.path.splitext(fname)[0], "wall"
    return clean_cap_name(fname), "cap"


def triangle_keys(polydata):
    """Yield (cell_index, frozenset of the cell's GlobalNodeIDs)."""
    gnid = polydata.GetPointData().GetArray(NODE_ID_ARRAY_NAME)
    if gnid is None:
        print(f"CRITICAL ERROR: '{NODE_ID_ARRAY_NAME}' point array missing on "
              f"{polydata.GetNumberOfPoints()}-point polydata. Cannot match "
              "faces without it.")
        sys.exit(1)
    for c in range(polydata.GetNumberOfCells()):
        pid = polydata.GetCell(c).GetPointIds()
        key = frozenset(int(gnid.GetTuple1(pid.GetId(k)))
                        for k in range(pid.GetNumberOfIds()))
        yield c, key


def group_faces(surfaces_dir, merge_walls):
    """Group surface files into faces. Returns a list of dicts with keys
    'name', 'type', 'files' (source filenames), ordered so face id = index + 1.

    Without merge_walls each file is its own face. With merge_walls every
    wall_* / wall_blend_* file collapses into a single 'wall' face; caps stay
    one face each."""
    files = sorted(f for f in os.listdir(surfaces_dir) if f.endswith(".vtp"))
    if not files:
        print(f"Error: no .vtp files in {surfaces_dir}")
        sys.exit(1)

    faces = []
    if merge_walls:
        wall_files = []
        for fname in files:
            if is_wall_file(fname):
                wall_files.append(fname)
            else:
                name, ftype = face_name_and_type(fname)
                faces.append({"name": name, "type": ftype, "files": [fname]})
        if wall_files:
            faces.append({"name": "wall", "type": "wall", "files": wall_files})
    else:
        for fname in files:
            name, ftype = face_name_and_type(fname)
            faces.append({"name": name, "type": ftype, "files": [fname]})

    guard_unique_names(faces)
    return faces


def guard_unique_names(faces):
    """Abort if two distinct faces resolve to the same name (e.g. clean_cap_name
    maps two different cap files onto one anatomical name). A duplicate would
    silently merge faces in the GUI, so refuse rather than guess. A merged
    'wall' face is a single entry owning many files and never trips this."""
    counts = {}
    for face in faces:
        counts[face["name"]] = counts.get(face["name"], 0) + 1
    dupes = {}
    for face in faces:
        if counts[face["name"]] > 1:
            dupes.setdefault(face["name"], []).extend(face["files"])
    if dupes:
        print("CRITICAL ERROR: face name collision(s) detected. These distinct "
              "surface files resolve to the same face name:")
        for name, files in sorted(dupes.items()):
            print(f"  '{name}' <- {', '.join(sorted(files))}")
        print("Rename the source .vtp files or adjust clean_cap_name() so each "
              "face name is unique, then re-run.")
        sys.exit(1)


def build_triple_map(surfaces_dir, faces):
    """Map each triangle's GlobalNodeID triple to its 1-based face id."""
    triple_to_face = {}
    collisions = 0
    total = sum(len(f["files"]) for f in faces)
    print(f"--- Building face partition: {total} surface files -> "
          f"{len(faces)} faces ---")
    for fid, face in enumerate(faces, start=1):
        for fname in face["files"]:
            pd = load_vtp(os.path.join(surfaces_dir, fname))
            for _, key in triangle_keys(pd):
                if key in triple_to_face and triple_to_face[key] != fid:
                    collisions += 1
                triple_to_face[key] = fid
    if collisions:
        print(f"Warning: {collisions} triangle(s) claimed by more than one "
              "face; later files win. Partition is not clean.")
    return triple_to_face


def label_model(global_pd, triple_to_face):
    """Add a ModelFaceID cell array to global_pd from the triple map.
    Returns the number of unmatched cells."""
    n = global_pd.GetNumberOfCells()
    arr = vtk.vtkIntArray()
    arr.SetName(FACE_ID_ARRAY_NAME)
    arr.SetNumberOfComponents(1)
    arr.SetNumberOfTuples(n)
    unmatched = 0
    for c, key in triangle_keys(global_pd):
        fid = triple_to_face.get(key)
        if fid is None:
            fid = 0
            unmatched += 1
        arr.SetTuple1(c, fid)
    global_pd.GetCellData().AddArray(arr)
    global_pd.GetCellData().SetActiveScalars(FACE_ID_ARRAY_NAME)
    return unmatched


def write_vtp(polydata, filepath):
    w = vtk.vtkXMLPolyDataWriter()
    w.SetFileName(filepath)
    w.SetInputData(polydata)
    w.SetDataModeToAppended()
    w.SetCompressorTypeToZLib()
    w.Write()


def build_mdl_tree(faces):
    """Build the MDL XML tree: one <face> per grouped face, ids 1..N.
    Schema verified against existing .mdl files and sv4gui_ModelIO."""
    root = ET.Element("model", {"type": "PolyData", "version": "1.0"})
    timestep = ET.SubElement(root, "timestep", {"id": "0"})
    me = ET.SubElement(timestep, "model_element", {
        "type": "PolyData",
        "num_sampling": "0",
        "use_uniform": "0",
    })
    ET.SubElement(me, "segmentations")
    faces_el = ET.SubElement(me, "faces")
    caps = walls = 0
    for idx, face in enumerate(faces):
        if face["type"] == "cap":
            caps += 1
        else:
            walls += 1
        ET.SubElement(faces_el, "face", {
            "id": str(idx + 1),
            "name": face["name"],
            "type": face["type"],
            "visible": "true",
            "opacity": "1",
            "color1": "1",
            "color2": "1",
            "color3": "1",
        })
    ET.SubElement(me, "blend_radii")
    ET.SubElement(me, "blend_param", {
        "blend_iters": "2",
        "sub_blend_iters": "3",
        "cstr_smooth_iters": "2",
        "lap_smooth_iters": "50",
        "subdivision_iters": "1",
        "decimation": "0.01",
    })
    return ET.ElementTree(root), caps, walls


def main(mesh_exterior_path, mesh_surfaces_dir, mdl_out, vtp_out, merge_walls):
    print("--- Starting Processing ---")
    print(f"Mesh Exterior: {mesh_exterior_path}")
    print(f"Mesh Surfaces Dir: {mesh_surfaces_dir}")
    print(f"Merge walls:  {merge_walls}")

    faces = group_faces(mesh_surfaces_dir, merge_walls)
    triple_to_face = build_triple_map(mesh_surfaces_dir, faces)

    global_pd = load_vtp(mesh_exterior_path)
    unmatched = label_model(global_pd, triple_to_face)
    print(f"Model cells: {global_pd.GetNumberOfCells()}, unmatched: {unmatched}")
    if unmatched:
        print("Warning: unmatched cells were assigned ModelFaceID 0. The model "
              "surface and mesh-surfaces partition may be inconsistent.")

    write_vtp(global_pd, vtp_out)
    print(f"Wrote labeled model VTP -> {vtp_out}")

    tree, caps, walls = build_mdl_tree(faces)
    ET.indent(tree, space="    ")
    with open(mdl_out, "wb") as f:
        f.write(b'<?xml version="1.0" encoding="UTF-8"?>\n')
        tree.write(f, encoding="UTF-8", xml_declaration=False)
    print(f"Wrote MDL -> {mdl_out}")
    print(f"\nSuccess! {len(faces)} faces: {caps} caps, {walls} walls.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Generate a capped SimVascular .mdl + labeled .vtp from a "
                    "model VTP and a mesh-surfaces/ folder, assigning faces "
                    "exactly from the surface partition (no feature angle).")
    parser.add_argument("--mesh-exterior", required=True, help="Path to mesh exterior .vtp")
    parser.add_argument("--mesh-surfaces-dir", required=True,
                        help="Path to mesh-surfaces/ folder")
    parser.add_argument("--model-name", required=True,
                        help="Base output model name, without .mdl or .vtp")
    parser.add_argument("--merge-walls", action="store_true",
                        help="Collapse all wall_* / wall_blend_* files into a "
                             "single 'wall' face. Default: each surface file "
                             "is its own face.")
    args    = parser.parse_args()
    mdl_out = args.model_name + ".mdl"
    out_vtp = args.model_name + ".vtp"
    main(args.mesh_exterior, args.mesh_surfaces_dir, mdl_out, out_vtp, args.merge_walls)
