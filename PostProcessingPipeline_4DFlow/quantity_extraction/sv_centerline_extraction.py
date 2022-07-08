'''Centerline Extraction using SimVascular Algorithm
'''
import os
import sys

import pdb
import sv
import vtk


def main(mdir, outputdir, total_timepts):
    # Create a modeler.
    kernel = sv.modeling.Kernel.POLYDATA
    modeler = sv.modeling.Modeler(kernel)

    # Obtain Face IDS from mdl files
    MPA_face_ids = []
    LPA_face_ids = []
    RPA_face_ids = []
    for i in range(0, total_timepts):
        mdl_file = 'def_%02d.mdl' % i
        with open(os.path.join(mdir, mdl_file), 'r') as f:
            for line in f:
                if "<face id" in line:
                    x = line.split()
                    id_num = int(x[1][4:-1])
                    if x[2][6:-1] == 'MPA':
                        MPA_face_ids.append(id_num)
                    if x[2][6:-1] == 'LPA':
                        LPA_face_ids.append(id_num)
                    if x[2][6:-1] == 'RPA':
                        RPA_face_ids.append(id_num)
        print('%s: %d %d %d' % (mdl_file, MPA_face_ids[-1], LPA_face_ids[-1], RPA_face_ids[-1]))

    print("Read surface model file ...")
    for i in range(0, total_timepts):
        # Read model geometry.
        file_name = 'def_%02d.vtp' % i
        model = modeler.read(os.path.join(mdir, file_name))
        model_polydata = model.get_polydata()
        print(file_name)
        print("\tModel: num nodes: {0:d}".format(model_polydata.GetNumberOfPoints()))
        #
        # Use node IDs.
        inlet_ids = [MPA_face_ids[i]]
        outlet_ids = [LPA_face_ids[i], RPA_face_ids[i]]
        centerlines_polydata = sv.vmtk.centerlines(model_polydata, inlet_ids, outlet_ids, use_face_ids=True)

        ## Write the centerlines
        if not (os.path.exists(os.path.join(outputdir, 'cl_sv'))):
            os.mkdir(os.path.join(outputdir, 'cl_sv'))
        file_name = os.path.join(outputdir, 'cl_sv', "centerlines-result_%02d.vtp" % i)
        writer = vtk.vtkXMLPolyDataWriter()
        writer.SetFileName(file_name)
        writer.SetInputData(centerlines_polydata)
        writer.Update()
        writer.Write()

        print("Centerlines: num nodes: {0:d}".format(centerlines_polydata.GetNumberOfPoints()))


"""
USER CALL
"""
# main()
# MAIN FUNCTION
if __name__ == "__main__":

    if len(sys.argv) < 4:
        print("Error: Arguments Missing")
        print("Usage:")
        print(" python {} [Model Dir] [Output Dir] [Total Timepoints]".format(sys.argv[0]))
        sys.exit(0)

    mdir = sys.argv[1]
    outputdir = sys.argv[2]
    total_timepts = sys.argv[3]
    if os.path.isdir(mdir):
        print('Model directory does not exist')
        sys.exit(0)
    if os.path.isdir(outputdir):
        print('Output directory does not exist.')
        sys.exit(0)

    sys.exit(main(mdir, outputdir, total_timepts))