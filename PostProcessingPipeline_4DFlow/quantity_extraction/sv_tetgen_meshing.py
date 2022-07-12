'''Test the TetGen class interface.
    Meshing a surface from a .vtk file with no ModelFaceID array.
    Usage:
        mesh-vtk-surf.py FILE_NAME.vtk
    Writes out FILE_NAME.vtu volume and FILE_NAME.vtk surface mesh file.
'''

import os
from pathlib import Path
import sv
import sys
import vtk
import time
import pdb

def main(edge_length, input_file, output_file):
    # Create a model and compute faces.
    modeler = sv.modeling.Modeler(sv.modeling.Kernel.POLYDATA)
    model = modeler.read(input_file)
    face_ids = model.compute_boundary_faces(angle=60.0)

    # Create a TetGen mesher.
    mesher = sv.meshing.TetGen()

    # Set the model for the mesher.
    mesher.set_model(model)

    # Set meshing options.
    options = sv.meshing.TetGenOptions(global_edge_size=edge_length, surface_mesh_flag=True, volume_mesh_flag=True)

    # Generate the mesh. 
    mesher.set_walls(face_ids)
    mesher.generate_mesh(options)

    # Get the mesh as a vtkUnstructuredGrid. 
    mesh = mesher.get_mesh()
    print("Mesh: %s" % output_file.split('/')[-1]);
    print("  Number of nodes: {0:d}".format(mesh.GetNumberOfPoints()))
    print("  Number of elements: {0:d}".format(mesh.GetNumberOfCells()))

    # Write volume mesh
    mesh_vol_file = output_file
    mesher.write_mesh(file_name=str(mesh_vol_file))

    # Write surface mesh
    mesh_surf_file = output_file[:-4] + '_surfacemesh.vtk'
    mesh_surface = mesher.get_surface()
    writer = vtk.vtkPolyDataWriter()
    writer.SetFileName(mesh_surf_file)
    writer.SetInputData(model)
    writer.Update()
    writer.Write()


"""
USER CALL
"""
# main()
# MAIN FUNCTION
if __name__ == "__main__":

    if len(sys.argv) != 4:
        print("Error: Arguments Missing")
        print("Usage:")
        print(" python {} [edge_length] [input_file] [output_file]".format(sys.argv[0]))
        sys.exit(0)

    edge_length = float(sys.argv[1])
    input_file = sys.argv[2]
    output_file = sys.argv[3]

    sys.exit(main(edge_length, input_file, output_file))