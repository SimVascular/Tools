# vtp_to_capped_mdl

Generates a SimVascular model — a labeled `.vtp` **and** its `.mdl` — from a model `.vtp` plus its `mesh-surfaces/` folder in one shot.

This script assigns a `ModelFaceID` to every model triangle directly from the `.vtp` files in the `mesh-surfaces/` folder, using the triangle's `GlobalNodeID` triple as an exact integer join key (no geometry, no tolerance).

### Two outputs, not one

The input model may not carry a `ModelFaceID` array (e.g., older `model.vtp` files in the Vascular Model Repository does not — it only has `GlobalElementID`, `GlobalNodeID`, `ModelRegionID`). The `.mdl` references faces by id, and those ids must exist as a `ModelFaceID` cell array on the model polydata. So the script writes **both**.

The script writes both outputs from `--model-name`:

- `<model-name>.vtp` — the model polydata with a freshly computed `ModelFaceID` array.
- `<model-name>.mdl` — the `.mdl` whose `<face>` ids match that array.

After the script is run, put both in your SimVascular project's `Models/` directory.

## Run

The script does not depend on the SimVascular Python API (`sv`) and can be
run in any Python environment with VTK available. Using SimVascular's bundled
Python may still be convenient because it already includes VTK.

On macOS:

```bash
/Applications/SimVascular.app/Contents/Resources/simvascular --python -- \
    vtp_to_capped_mdl.py \
    --mesh-exterior  mesh.exterior.vtp \
    --mesh-surfaces-dir mesh-surfaces \
    --model-name model
```

Or on Linux:
```
/usr/local/sv/2025-12-21/simvascular --python -- \
    vtp_to_capped_mdl.py \
    --mesh-exterior mesh.exterior.vtp \
    --mesh-surfaces-dir mesh-surfaces \
    --model-name model
```

Use the `Resources/simvascular` wrapper, not the `bin/simvascular` binary directly — the wrapper sets `DYLD_LIBRARY_PATH` so `lib_simvascular_post.dylib` resolves. On Linux/Windows: substitute the equivalent `simvascular` launcher from your install.

Alternatively, if you have a Python environment with VTK set up, you can run the script directly:
```
python vtp_to_capped_mdl.py \
    --mesh-exterior mesh.exterior.vtp --mesh-surfaces-dir mesh-surfaces --model-name model
```

## Arguments

- `--mesh-exterior` — model `.vtp`. Needs a `GlobalNodeID` point array (any SimVascular-meshed surface has one). A `ModelFaceID` array is **not** required; the script creates it.
- `--mesh-surfaces-dir` — folder of the `mesh-surfaces` folder, which includes the per-face `.vtp` files. Files starting with `wall_` or `wall_blend_` are walls; everything else is a cap.
- `--model-name` — base output path/name without extension. For example, `--model-name out/model` writes `out/model.mdl` and `out/model.vtp`.
- `--merge-walls` — collapse all wall files into a single `wall` face. Off by default (each surface file becomes its own face).

## Face naming

Currently, this follows the typical naming convention for pulmonary models.

- **Caps** get a cleaned anatomical name via `clean_cap_name` (`l_pa_4_1_x.vtp` → `lpa_4_1`); `inflow.vtp` and `r_pa_x_2.vtp` both map to `inflow`.
- **Walls** keep their filename stem (`wall_RPA_07_01.vtp` → `wall_RPA_07_01`), unless `--merge-walls`, in which case they all become one face named `wall`.

### Collision guard

If two distinct surface files resolve to the same face name (e.g. `clean_cap_name` maps two cap files onto one name), the script **aborts** with a `CRITICAL ERROR` listing the offending files rather than silently merging them. Rename the source `.vtp` files or adjust `clean_cap_name()` and re-run. (A `--merge-walls` `wall` face owns many files by design and does not trip this.)

## Example: model

**Inputs**

```
mesh.exterior.vtp   # model VTP (139,948 triangles), no ModelFaceID
mesh-surfaces/      # 272 .vtp files: 92 caps + 180 wall_* / wall_blend_*
```

**Command (default — each file is its own face)**

```bash
/Applications/SimVascular.app/Contents/Resources/simvascular --python -- \
    vtp_to_capped_mdl.py \
    --mesh-exterior mesh.exterior.vtp --mesh-surfaces-dir mesh-surfaces \
    --model-name out/model
```

Tail of stdout:

```
--- Building face partition: 272 surface files -> 272 faces ---
Model cells: 139948, unmatched: 0
Wrote labeled model VTP -> out/model.vtp
Wrote MDL -> out/model.mdl

Success! 272 faces: 92 caps, 180 walls.
```

**Command (merged walls)**

Add `--merge-walls` to collapse the 180 wall pieces into one face:

```
--- Building face partition: 272 surface files -> 93 faces ---
Success! 93 faces: 92 caps, 1 walls.
```

**Output `.mdl` (first lines)**

```xml
<?xml version="1.0" encoding="UTF-8"?>
<model type="PolyData" version="1.0">
    <timestep id="0">
        <model_element type="PolyData" num_sampling="0" use_uniform="0">
            <segmentations />
            <faces>
                <face id="1" name="LPA"     type="cap"  visible="true" opacity="1" color1="1" color2="1" color3="1" />
                <face id="2" name="LPA01"   type="cap"  visible="true" opacity="1" color1="1" color2="1" color3="1" />
                ...
                <face id="93" name="wall_LPA" type="wall" visible="true" opacity="1" color1="1" color2="1" color3="1" />
                ...
```

Caps show with `type=cap`, walls with `type=wall` — ready to assign BCs in the Model tab.
