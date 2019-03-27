# Contributed by Aslak W. Bergersen (aslak.bergersen@gmail.com)

"""
This is an example script of how one could convert a mesh output from
SimVascular into a format which is redeble in FEniCS.
"""

from os import path
from argparse import ArgumentParser

try:
    import vtk
except:
    raise ImportError("Could not find vtk, please install using pip or conda")

try:
    import meshio
except:
    raise ImportError("Could not find meshio, please install by running" + \
                      " pip install meshio")


def read_data(filename):
    """
    Load the given file, and return a vtkPolyData object for it.

    Args:
        filename (str): Path to input file.
        datatype (str): Additional parameter for vtkIdList objects.

    Returns:
        Data (vtkSTL/vtkPolyData/vtkXMLStructured/
                    vtkXMLRectilinear/vtkXMLPolydata/vtkXMLUnstructured/
                    vtkXMLImage/Tecplot): Output data.
    """
    # Check if file exists
    if not path.exists(filename):
        raise NameError("Could not find file: %s" % filename)

    # Check filename format
    file_type = filename.split(".")[-1]
    if file_type == '':
        raise NameError('The file does not have an extension')

    # Get reader
    if file_type == 'vtp':
        reader = vtk.vtkXMLPolyDataReader()
    elif file_type == "stl":
        reader = vtk.vtkSTLReader()
    elif file_type == 'vtu':
        reader = vtk.vtkXMLUnstructuredGridReader()
    else:
        raise NameError('Unknown file type %s' % fileType)

    # Read
    reader.SetFileName(filename)
    reader.Update()
    data = reader.GetOutput()

    return data


def write_data(input_data, filename):
    """
    Write the given input data based on the file name extension.

    Args:
        input_data (vtkSTL/vtkPolyData/vtkXMLPolydata/vtkXMLUnstructured): Input data.
        filename (str): Save path location.
    """
    # Check filename format
    fileType = filename.split(".")[-1]
    if fileType == '':
        raise NameError('The file does not have an extension')

    # Get writer
    if fileType == 'stl':
        writer = vtk.vtkSTLWriter()
    elif fileType == 'vtp':
        writer = vtk.vtkXMLPolyDataWriter()
    elif fileType == 'vtu':
        writer = vtk.vtkXMLUnstructuredGridWriter()
    else:
        raise NameError('Unknown file type %s' % fileType)

    # Set filename and input
    writer.SetFileName(filename)
    writer.SetInputData(input_data)
    writer.Update()

    # Write
    writer.Write()


def convert_vtp_to_vtu(vtp):
    """Wrapper around the vtkAppendFilter for appending data sets into on vtkUnstructedGrid.
    Here it is used to convert only one vtkPolyData object into a vtkUnstructuedGrid.

    Args:
        vtp (vtkPolyData): Surface to convert.

    Returns:
        vtu (vtkUnstructuredGrid): Converted surface.
    """
    append_filter = vtk.vtkAppendFilter()
    append_filter.AddInputData(vtp)
    append_filter.Update()

    return append_filter.GetOutput()


def convert_to_dolfin_xdmf(vtp_path, vtu_path, output_path, global_node_ID_name, face_ID_name):
    """Convert the vtu file in output_path to two .xdmf files which can be read in from FEniCS.

    Args:
        vtp_path (str): Path to the vtp file outputed from SimVascular.
        vtu_path (str): Path to the vtu file outputed from SimVascular.
        output_path (str): Path to the merged file.
        global_node_ID_name (str): Name of the array holding the global cell ID node
        face_ID_name (str): Name of the array holding the boundary markers.
    """
    # Convert the vtkPolyData to an vtkUnstructuredGrid
    vtp = read_data(vtp_path)
    vtu_surface = convert_vtp_to_vtu(vtp)
    write_data(vtu_surface, vtp_path.replace(".vtp", "_converted_vtp.vtu"))

    msh_vtu = meshio.read(vtu_path)
    msh_vtp = meshio.read(vtp_path.replace(".vtp", "_converted_vtp.vtu"))

    # Convert volume mesh, equivalent of the commandline:
    # > meshio-convert inputfile.vtu outputfile.xdmf
    meshio.write(output_path, msh_vtu)

    # Node map from the vtp file to the vtu mesh
    vtp_global_node_ID = msh_vtp.point_data[global_node_ID_name] - 1

    # Map the values from the surface to the mesh
    N = msh_vtp.cells["triangle"].shape[0]
    cells = msh_vtp.cells["triangle"].flatten()
    mapped_cells = vtp_global_node_ID[cells].reshape((N, 3))
    cell_data = msh_vtp.cell_data["triangle"][face_ID_name]

    # Write the boundary markers
    facet_path = output_path.replace(".xdmf", "_facet_markers.xdmf")
    meshio.write(facet_path, meshio.Mesh(points=msh_vtu.points,
                                         cells={"triangle": mapped_cells},
                                         cell_data={"triangle": {face_ID_name: cell_data}}))


def read_command_line():
    description = """Convert the output mesh from SimVascular to a format
    compatible with FEniCS."""
    parser = ArgumentParser(description=description)

    # Required arguments
    parser.add_argument("-ivtp", "--vtp-path", type=str, required=True,
                        help="Path to the vtp file")
    parser.add_argument("-ivtu", "--vtu-path", type=str, required=True,
                        help="Path to the vtu file")
    parser.add_argument("-o", "--output-path", type=str, required=True,
                        help="Path to store the new xdmf file")

    # Optional argments
    parser.add_argument("-n", "--global-node-ID-name", type=str, default="GlobalNodeID",
                        help="Name of the global cell ID array")
    parser.add_argument("-f", "--face-ID-name", type=str, default="ModelFaceID",
                        help="Name of the face ID array")

    args = parser.parse_args()

    if not args.output_path.endswith(".xdmf"):
        raise NameError("The output path to have an .xdmf extension")
    if not args.vtp_path.endswith(".vtp"):
        raise NameError("The input path to the VTP file has to have an .vtp extension")
    if not args.vtu_path.endswith(".vtu"):
        raise NameError("The input path to the VTU file has to have an .vtu extension")

    return args.__dict__


if __name__ == "__main__":
    convert_to_dolfin_xdmf(**read_command_line())
