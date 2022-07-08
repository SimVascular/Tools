"""Helper Functions for 4DMRI analysis - Shape registration, remeshing, centerline extraction, and vtk
Requires the following additional python packages:
- deformetrica
=================================================================================================================================
---Created--- | ------Owner------- | Notes---------------------------------------------------------------------------------------
=================================================================================================================================
   11-23-2020   Melody Dong          Created
=================================================================================================================================
"""
import os
from os import listdir, path
from os.path import isfile, join
import sys

import deformetrica as dfca
import multiprocessing as mp
import numpy as np
import pandas as pd
import pdb
import subprocess
import time
import vtk
from vmtk import vtkvmtk, vmtkscripts
from vtk.util import numpy_support
import vtk.util.numpy_support as nps
import warnings
warnings.filterwarnings("ignore")

# Add SimVascular to path
sys.path.insert(0, '/usr/local/sv/simvascular/2021-02-06/simvascular')


"""PUBLIC FUNCTIONS"""
def centerline_extraction(model_name, cl_outputfile):
	"""Compute the centerlines of the model using VMTK's centerline extraction. If model given is a volume
	mesh, convert to a surface mesh for visualization in VMTK's gui. User must specify the source and target
	points for each model provided.
	Args:
		model_name: full path and file of the geometry to extract centerlines from
		cl_outputfile: full path and file for saving the centerline output to
	Return:
		centerlines_output: object containing vmtk centerlines
	"""
	model = read_polydata(model_name)
	# Convert vtu to vtk for vmtk centerline extraction
	if 'vtu' in model_name:
		vtp_filter = vtk.vtkGeometryFilter()
		vtp_filter.SetInputData(model)
		vtp_filter.Update()
		model_vtk = model_name[:-4] + '_surfacemesh.vtk'
		write_polydata(vtp_filter.GetOutput(), model_vtk)
		model = read_polydata(model_vtk)

	# Centerline Extraction
	centerlines = vmtkscripts.vmtkCenterlines()
	centerlines.Surface = model
	centerlines.SeedSelectorName = "pickpoint"
	centerlines.AppendEndPoints = 1
	centerlines.Resampling = 1
	centerlines.ResamplingStepLength = 0.01 # distance between points in resampled line. 10% of outlet diameter A=0.02 cm^2 D=0.159
	centerlines.Execute()
	centerlines_output = centerlines.Centerlines

	# Centerline Smoothing
	cl_smoothing = vmtkscripts.vmtkCenterlineSmoothing()
	cl_smoothing.Centerlines = centerlines_output
	cl_smoothing.NumberofSmoothingIterations = 500
	cl_smoothing.SmoothFactor = 0.2
	cl_smoothing.Execute()
	centerlines_output = cl_smoothing.Centerlines

	# Branch identification - bifurcations
	branchextractor = vmtkscripts.vmtkBranchExtractor()
	branchextractor.Centerlines = centerlines_output
	branchextractor.Execute()
	centerlines_output = branchextractor.Centerlines
	print("WRITING OUT CENTERLINE FILE")

	# Write out centerlines to file
	cl_geometry = vmtkscripts.vmtkCenterlineGeometry()
	cl_geometry.Centerlines = centerlines_output
	cl_geometry.Execute()
	centerlines_output = cl_geometry.Centerlines
	write_polydata(centerlines_output, cl_outputfile)

	return centerlines_output


def compute_vortical_quantities(model_name):
	"""Compute other quantities of interest such as vorticity, Q-criterion, Helicity density
	Use vtkGradientFilter() to calculate the gradients of velocity, vorticity, Q-criterion. Should
	default to computing on an Unstructured Grid DataSet if model is a vtu
	"""
	model = read_polydata(model_name)

	# Compute the gradients, vorticity, and Q-criterion
	grad = vtk.vtkGradientFilter()
	grad.SetInputData(model)
	grad.SetInputScalars(0, 'velocity') # 0=pointdata, 'velocity'=input to compute gradients
	grad.SetResultArrayName('gradient_u')
	grad.ComputeVorticityOn()
	grad.ComputeQCriterionOn()
	grad.SetVorticityArrayName('vorticity')
	grad.SetQCriterionArrayName('Q-criterion')
	grad.Update()

	# Convert results to vtkArray and update model
	gradu_array = vtk.vtkDoubleArray()
	vort_array = vtk.vtkDoubleArray()
	qcrit_array = vtk.vtkDoubleArray()
	gradu_array = grad.GetOutput().GetPointData().GetArray('gradient_u')
	vort_array = grad.GetOutput().GetPointData().GetArray('vorticity')
	qcrit_array = grad.GetOutput().GetPointData().GetArray('Q-criterion')
	model.GetPointData().AddArray(gradu_array)
	model.GetPointData().AddArray(vort_array)
	model.GetPointData().AddArray(qcrit_array)

	# Compute Helicity Density and Relative Helicity Density
	velocity_np = nps.vtk_to_numpy(model.GetPointData().GetArray('velocity'))
	vort_np = nps.vtk_to_numpy(vort_array)
	helicity_dens_array = vtk.vtkFloatArray()
	helicity_dens_array.SetNumberOfComponents(1)
	helicity_dens_array.SetName("helicity_density")
	rel_helicity_array = vtk.vtkFloatArray()
	rel_helicity_array.SetNumberOfComponents(1)
	rel_helicity_array.SetName("relative_helicity_density")
	for pointId in range(model.GetNumberOfPoints()):
		helicity_dens_array.InsertNextValue(np.dot(velocity_np[pointId], vort_np[pointId]))
		rel_helicity_array.InsertNextValue(np.dot(velocity_np[pointId], vort_np[pointId])/(np.linalg.norm(velocity_np[pointId])*np.linalg.norm(vort_np[pointId])))
	model.GetPointData().AddArray(helicity_dens_array)
	model.GetPointData().AddArray(rel_helicity_array)

	# Viscous Dissipation - for calculating energy loss
	gradu_np = nps.vtk_to_numpy(gradu_array)
	div_vel = gradu_np[:, 0] + gradu_np[:, 4] + gradu_np[:, 8]
	viscous_diss = 0.5*((gradu_np[:, 0] - (2/3)*div_vel)**2 + 2*(gradu_np[:,3] + gradu_np[:,1])**2 + 
						(gradu_np[:, 4] - (2/3)*div_vel)**2 + 2*(gradu_np[:,5] + gradu_np[:,7])**2 + 
						(gradu_np[:, 8] - (2/3)*div_vel)**2 + 2*(gradu_np[:,6] + gradu_np[:,2])**2)
	viscous_diss_array = vtk.vtkFloatArray()
	viscous_diss_array = nps.numpy_to_vtk(viscous_diss)
	viscous_diss_array.SetName("viscous_dissipation")
	model.GetPointData().AddArray(viscous_diss_array)

	write_polydata(model, model_name)


def extract_qoi(ModelName, centerlines_output_file, seg_location, cl_type):
	"""Extract quantities of interest from cross-sectional slices based on the centerlines for each group segment
	Quantities of interest include: Flow, Vorticity, Q-Criterion, and Area
	Args:
		ModelName: full path file to the 3D data set that includes all hemodynamic quantities
		centerlines_output: centerlines vtp file containing groupIDs categorizing vessel segments
		seg_location: location from 0 to 1 down length of vessel segment where slice is extracted
	"""
	centerlines_output = read_polydata(centerlines_output_file)

	# Resample SimVascular centerlines because not refined enough - uses splines
	if cl_type.lower() == 'sv':
		cl_resample = vmtkscripts.vmtkCenterlineResampling()
		cl_resample.Centerlines = centerlines_output
		cl_resample.Length = 0.01
		cl_resample.Execute()
		centerlines_output = cl_resample.Centerlines
		write_polydata(centerlines_output, ModelName[0:-4] + "_clresample.vtp")

	# Read Centerlines
	num_group, num_path, num_cells, group_elems = read_centerlines(centerlines_output, cl_type)
	group_roi_center, group_roi_tan, group_roi_maxR = calculate_group_roi(centerlines_output, num_group, group_elems, seg_location, cl_type)
	
	# Integrate and extract info
	results = cut_integrate_slice(ModelName, num_group, group_roi_center, group_roi_tan, group_roi_maxR, seg_location, cl_type)
	cl_results = cl_smallcut_integrate_slice(ModelName, num_group, group_roi_center, group_roi_tan, group_roi_maxR, seg_location, cl_type)
	cl_results = cl_results.drop(['group_Num'], axis=1)
	results = pd.concat([results, cl_results], axis=1)

	return results

def integrate_slice(vtu_slice_file):
	"""Integrate a given slice and return integrated outputs"""
	vtu_slice = read_polydata(vtu_slice_file)

	# Integrate quantities of interest over slice
	vn = [0.0,0.0,0.0]
	prev_vn = vn
	slice_area = 0.0
	slice_flow = 0.0
	slice_maxvel = 0.0
	slice_vort_avg = 0.0
	slice_q_avg = 0.0
	slice_hd_avg = 0.0
	slice_relhd_avg = 0.0
	slice_uorth_avg = 0.0
	sep_bubble_area = 0.0
	q_crit_vel_area = 0.0
	velocity_array = vtu_slice.GetPointData().GetArray("velocity")
	vorticity_array = vtu_slice.GetPointData().GetArray("vorticity")
	qcrit_array = vtu_slice.GetPointData().GetArray("Q-criterion")
	hel_dens_array = vtu_slice.GetPointData().GetArray("helicity_density")
	rel_hel_dens_array = vtu_slice.GetPointData().GetArray("relative_helicity_density")
	u_orth_array = vtu_slice.GetPointData().GetArray("u.n")

	missing_cell = []
	for cellId in range(vtu_slice.GetNumberOfCells()):
		cell = vtu_slice.GetCell(cellId)
		p0 = cell.GetPoints().GetPoint(0)
		p1 = cell.GetPoints().GetPoint(1)
		# For 4DMRI data, if slice goes through centers slices, it will only include full cells with 3 points
		try:
			p2 = cell.GetPoints().GetPoint(2) 
		except:
			# print('ERROR: Point 2 not found in cell %2d' % cellId)
			missing_cell.append(cellId)
			continue
		cell_area = vtk.vtkTriangle().TriangleArea(p0, p1, p2)
		slice_area = slice_area + cell_area
		nodeids = cell.GetPointIds()
		# Compute Flow
		v1 = np.array(velocity_array.GetTuple3(nodeids.GetId(0)))
		v2 = np.array(velocity_array.GetTuple3(nodeids.GetId(1)))
		v3 = np.array(velocity_array.GetTuple3(nodeids.GetId(2)))
		vtk.vtkTriangle().ComputeNormal(p0, p1, p2, vn)
		if prev_vn != vn:
			print("Previous cell does not have same normal as current cell: " + str(vn))
		prev_vn = vn
		slice_flow = slice_flow + sum((v1+v2+v3)/3.0*vn)*cell_area # 1 point Gauss quad rule
		# Maximum velocity
		v1_mag = np.linalg.norm(v1)
		v2_mag = np.linalg.norm(v2)
		v3_mag = np.linalg.norm(v3)
		if np.mean([v1_mag, v2_mag, v3_mag]) > slice_maxvel:
			slice_maxvel =  np.mean([v1_mag, v2_mag, v3_mag])
		# Compute Vorticity Magnitude
		w1 = np.array(vorticity_array.GetTuple3(nodeids.GetId(0)))
		w2 = np.array(vorticity_array.GetTuple3(nodeids.GetId(1)))
		w3 = np.array(vorticity_array.GetTuple3(nodeids.GetId(2)))
		w1_mag = np.linalg.norm(w1)
		w2_mag = np.linalg.norm(w2)
		w3_mag = np.linalg.norm(w3)
		slice_vort_avg = slice_vort_avg + np.mean([w1_mag, w2_mag, w3_mag])*cell_area
		# Compute Q-criterion
		q1 = np.array(qcrit_array.GetTuple(nodeids.GetId(0)))
		q2 = np.array(qcrit_array.GetTuple(nodeids.GetId(1)))
		q3 = np.array(qcrit_array.GetTuple(nodeids.GetId(2)))
		slice_q_avg = slice_q_avg + np.mean([q1, q2, q3])*cell_area
		# Compute Helicity Density
		hd1 = np.array(hel_dens_array.GetTuple(nodeids.GetId(0)))
		hd2 = np.array(hel_dens_array.GetTuple(nodeids.GetId(1)))
		hd3 = np.array(hel_dens_array.GetTuple(nodeids.GetId(2)))
		slice_hd_avg = slice_hd_avg + np.mean([hd1, hd2, hd3])*cell_area
		# Compute Relative Helicity density
		relhd1 = np.array(rel_hel_dens_array.GetTuple(nodeids.GetId(0)))
		relhd2 = np.array(rel_hel_dens_array.GetTuple(nodeids.GetId(1)))
		relhd3 = np.array(rel_hel_dens_array.GetTuple(nodeids.GetId(2)))
		slice_relhd_avg = slice_relhd_avg + np.mean([relhd1, relhd2, relhd3])*cell_area
		# Compute u orthogonal (u . n) - if velocity vector is orthogonal to forward flow direction
		uorth1 = np.array(u_orth_array.GetTuple(nodeids.GetId(0)))
		uorth2 = np.array(u_orth_array.GetTuple(nodeids.GetId(1)))
		uorth3 = np.array(u_orth_array.GetTuple(nodeids.GetId(2)))
		slice_uorth_avg = slice_uorth_avg + np.mean([uorth1, uorth2, uorth3])*cell_area
		# Threshold u orthogonal (u.n) - if velocity vector is orthogonal or negative to forward flow direction (<1 to account for error)
		if np.mean([uorth1, uorth2, uorth3]) <= 1:
			sep_bubble_area += cell_area
		# Threshold Q-criterion - if greater than 0, presence of a vortex
		if np.mean([q1, q2, q3]) > 0:
			q_crit_vel_area += cell_area


	# print("WARNING: The %d out of %d total cells were not included in the slice calculation due to a missing point" % (len(missing_cell), vtu_slice.GetNumberOfCells()))
	# print(str(missing_cell))
	integrated_results = [slice_area, slice_flow, slice_vort_avg, slice_q_avg, slice_hd_avg, 
				   		  slice_relhd_avg, slice_uorth_avg, sep_bubble_area, q_crit_vel_area, slice_maxvel]

	return integrated_results


def registration(name, subj_id, tmpl_name, data_name, max_iterations, timepoints, kernel_width):
	""" Registration - Creates geometry interpolate between template and target (data_name) using Deformetrica's shape analysis"""
	# Convert vtp to vtk for Deformetrica
	if tmpl_name[-3:] == 'vtp' or data_name[-3:] == 'vtp':
		tmpl_model = read_polydata(tmpl_name)
		data_model = read_polydata(data_name)
		vtp_filter = vtk.vtkGeometryFilter()
		vtp_filter.SetInputData(tmpl_model)
		vtp_filter.Update()
		tmpl_vtk = tmpl_name[:-4] + '_surfacemesh.vtk'
		write_polydata(vtp_filter.GetOutput(), tmpl_vtk)
		vtp_filter.SetInputData(data_model)
		vtp_filter.Update()
		data_vtk = data_name[:-4] + '_surfacemesh.vtk'
		write_polydata(vtp_filter.GetOutput(), data_vtk)
		data_name = data_vtk
		tmpl_name = tmpl_vtk

	dataset_specifications = {
		'dataset_filenames': [[{name: data_name}]],
		'subject_ids': [subj_id]
	}

	template_specifications = {
		name: {'deformable_object_type': 'SurfaceMesh',
			   'noise_std': 0.05,
			   'kernel_type': 'keops',
			   'kernel_width': kernel_width,
			   'filename': tmpl_name,
			   'attachment_type': 'current'}
	}
	estimator_options = {'optimization_method_type': 'ScipyLBFGS', 'max_iterations': max_iterations}

	data_path = '/'.join(data_name.split('/')[0:-1])
	deformetrica = dfca.Deformetrica(output_dir=os.path.join(data_path, 'registration_output'), verbosity='INFO')

	deformetrica.estimate_registration(
		template_specifications, dataset_specifications,
		estimator_options=estimator_options,
		model_options={'deformation_kernel_type': 'keops', 'deformation_kernel_width': kernel_width, 'number_of_time_points': timepoints})


def resample_4dmr(vtu_file, flowmri_file, data_path):
	"""Resample 4DMRI full scan with velocity to constrain only to MPA vtu mesh file"""
	source_vtk = read_polydata(flowmri_file, 'Structured')
	dest_mesh = read_polydata(vtu_file)
	resample = vtk.vtkResampleWithDataSet()
	resample.AddInputData(dest_mesh)
	resample.SetSourceData(source_vtk)
	resample.Update()
	output_file = flowmri_file.split('/')[-1]
	write_polydata(resample.GetOutput(), os.path.join(data_path, output_file[:-4] + '_resampled.vtu'))

	return(flowmri_file)


def volume_remesh(vtk_file, edge_length, output_file, data_path, mesh_type):
	"""Uses vmtk volume mesher to remesh the surface vtk outputted from Deformetrica
	Args:
		vtk_file (str): Path to the input vtk file
		edge_length (float): Target edge length for meshing
		output_file (str): Path to the output vtu file to write out remeshed geometry
		data_path (str): Path to where to write output vtu file to
		mesh_type (str): "vmtk" or "sv" for which mesher to call
	"""
	print(output_file)
	if mesh_type == 'vmtk':
		model = read_polydata(os.path.join(data_path, 'registration_output', vtk_file))
		vmtk_mesh = vmtkscripts.vmtkMeshGenerator()
		vmtk_mesh.Surface = model
		vmtk_mesh.TargetEdgeLength = edge_length
		vmtk_mesh.Execute()
		# If error occurs and only surface mesh will be generated
		if vmtk_mesh.Mesh.GetNumberOfCells() < 50000:
			print(f'ERROR: Surface mesh generated only - {vmtk_mesh.Mesh.GetNumberOfCells()}\n\n')
			return "ERROR"
		write_polydata(vmtk_mesh.Mesh, os.path.join(data_path, 'registration_output', output_file))

		vtp_filter = vtk.vtkGeometryFilter()
		vtp_filter.SetInputData(vmtk_mesh.Mesh)
		vtp_filter.Update()
		model_vtk = output_file[:-4] + '_surfacemesh.vtk'
		write_polydata(vtp_filter.GetOutput(), os.path.join(data_path, 'registration_output', model_vtk))

	elif mesh_type == 'sv':
		input_vtk_path = os.path.join(data_path, 'registration_output', vtk_file)
		output_vtk_path = os.path.join(data_path, 'registration_output', output_file)
		mesher_script = os.path.join(sys.path[1], 'sv_tetgen_meshing.py')
		sv_cmd = [sys.path[0], "--python", "--", mesher_script, str(edge_length), input_vtk_path, output_vtk_path]

		# Call SimVascular externally and write outputs to file
		FNULL = open(os.devnull, 'w')
		subprocess.call(sv_cmd, stdout=FNULL, stderr=subprocess.STDOUT)
		FNULL.close()

		# Check that mesh was outputted
		if path.exists(output_file):
			vtk_file = vtk_file
		else:
			vtk_file = "ERROR"

	return vtk_file


def read_polydata(filename, *datatype):
	"""
	Load the given file, and return a vtkPolyData object for it.
	Args:
		filename (str): Path to input file.
		datatype (str): Additional parameter for vtkIdList objects.
	Returns:
		polyData (vtkSTL/vtkPolyData/vtkXMLStructured/
				  vtkXMLRectilinear/vtkXMLPolydata/vtkXMLUnstructured/
				  vtkXMLImage/Tecplot): Output data.
	"""

	# Check if file exists
	if not path.exists(filename):
		raise RuntimeError("Could not find file: %s" % filename)

	# Check filename format
	fileType = filename.split(".")[-1]
	if fileType == '':
		raise RuntimeError('The file does not have an extension')

	if datatype == ():
		datatype = ('', '')

	# Get reader
	if fileType == 'stl':
		reader = vtk.vtkSTLReader()
		reader.MergingOn()
	elif fileType == 'vtk' and datatype[0].lower() == 'structured':
		reader = vtk.vtkStructuredGridReader()
	elif fileType == 'vtk':
		reader = vtk.vtkPolyDataReader()
	elif fileType == 'vtp':
		reader = vtk.vtkXMLPolyDataReader()
	elif fileType == 'vts':
		reader = vtk.vtkXMinkorporereLStructuredGridReader()
	elif fileType == 'vtr':
		reader = vtk.vtkXMLRectilinearGridReader()
	elif fileType == 'vtu':
		reader = vtk.vtkXMLUnstructuredGridReader()
	elif fileType == "vti":
		reader = vtk.vtkXMLImageDataReader()
	elif fileType == "np" and datatype[0] == "vtkIdList":
		result = np.load(filename).astype(np.int)
		id_list = vtk.vtkIdList()
		id_list.SetNumberOfIds(result.shape[0])
		for i in range(result.shape[0]):
			id_list.SetId(i, result[i])
		return id_list
	else:
		raise RuntimeError('Unknown file type %s' % fileType)

	# Read
	reader.SetFileName(filename)
	reader.Update()
	polydata = reader.GetOutput()

	return polydata


def write_polydata(input_data, filename, *datatype):
	"""
	Write the given input data based on the file name extension.
	Args:
		input_data (vtkSTL/vtkPolyData/vtkXMLStructured/
					vtkXMLRectilinear/vtkXMLPolydata/vtkXMLUnstructured/
					vtkXMLImage/Tecplot): Input data.
		filename (str): Save path location.
		datatype (str): Additional parameter for vtkIdList objects.
	"""
	# Check filename format
	fileType = filename.split(".")[-1]
	if fileType == '':
		raise RuntimeError('The file does not have an extension')

	if datatype == ():
		datatype = ('', '')

	# Get writer
	if fileType == 'stl':
		writer = vtk.vtkSTLWriter()
	elif fileType == 'vtk' and datatype[0].lower() == 'structured':
		writer = vtk.vtkStructuredGridWriter()
	elif fileType == 'vtk':
		writer = vtk.vtkPolyDataWriter()
	elif fileType == 'vts':
		writer = vtk.vtkXMLStructuredGridWriter()
	elif fileType == 'vtr':
		writer = vtk.vtkXMLRectilinearGridWriter()
	elif fileType == 'vtp':
		writer = vtk.vtkXMLPolyDataWriter()
	elif fileType == 'vtu':
		writer = vtk.vtkXMLUnstructuredGridWriter()
	elif fileType == "vti":
		writer = vtk.vtkXMLImageDataWriter()
	elif fileType == "np" and datatype == "vtkIdList":
		output_data = np.zeros(input_data.GetNumberOfIds())
		for i in range(input_data.GetNumberOfIds()):
			output_data[i] = input_data.GetId(i)
		output_data.dump(filename)
		return
	else:
		raise RuntimeError('Unknown file type %s' % fileType)

	# Set filename and input
	writer.SetFileName(filename)
	writer.SetInputData(input_data)
	writer.Update()

	# Write
	writer.Write()


def merged_centerline(model_name, cl_sv_file, data_path):
	"""Created merged centerlines"""
	model = read_polydata(model_name[:-4] + '_surfacemesh.vtk')
	cl_sv = read_polydata(cl_sv_file)

	# Create VMTK centerlines based on the inlet and outlet points from SimVascular centerlines
	# Find inlet and outlet coordinates for seeding
	num_group, num_path, num_cells, group_elems = read_centerlines(cl_sv, 'sv')
	num_pts = cl_sv.GetNumberOfPoints()
	group_list = nps.vtk_to_numpy(cl_sv.GetPointData().GetArray('BranchId'))
	path_group = nps.vtk_to_numpy(cl_sv.GetPointData().GetArray('Path'))
	center_pts = []
	for i in range(0, num_group):
		ids = [j for j in range(0, num_pts) if group_list[j] == group_elems[i]]
		path_sort = path_group[ids]
		# Sort ids based on a previously calculated path since resampled ids out of order
		ids = [x for _,x in sorted(zip(path_sort, ids))]
		num_ids = len(ids)
		if i == 0:
			center_pts.append(np.array(cl_sv.GetPoints().GetPoint(ids[0])))
		else:
			center_pts.append(np.array(cl_sv.GetPoints().GetPoint(ids[-1])))
	# print(center_pts)

	# Generate VMTK centerlines
	cl_vmtk = vmtkscripts.vmtkCenterlines()
	cl_vmtk.Surface = model
	cl_vmtk.SeedSelectorName = "pointlist"
	cl_vmtk.AppendEndPoints = 1
	cl_vmtk.SourcePoints = center_pts[0].tolist()
	cl_vmtk.TargetPoints = center_pts[1].tolist() + center_pts[2].tolist()
	cl_vmtk.Execute()
	cl_vmtk_output = cl_vmtk.Centerlines

	# Smooth Centerlines
	cl_smoothing = vmtkscripts.vmtkCenterlineSmoothing()
	cl_smoothing.Centerlines = cl_vmtk_output
	cl_smoothing.NumberofSmoothingIterations = 500
	cl_smoothing.SmoothFactor = 0.5
	cl_smoothing.Execute()
	cl_vmtk_output = cl_smoothing.Centerlines

	# Branch Extractor
	branchextractor = vmtkscripts.vmtkBranchExtractor()
	branchextractor.Centerlines = cl_vmtk_output
	branchextractor.Execute()
	cl_vmtk_output = branchextractor.Centerlines

	# Create Merged Centerlines
	cl_vmtk_merge = vmtkscripts.vmtkCenterlineMerge()
	cl_vmtk_merge.Centerlines = cl_vmtk_output
	cl_vmtk_merge.Execute()
	cl_vmtk_output = cl_vmtk_merge.Centerlines

	# Resampling
	cl_resample = vmtkscripts.vmtkCenterlineResampling()
	cl_resample.Centerlines = cl_vmtk_output
	cl_resample.Length = 0.01
	cl_resample.Execute()
	cl_vmtk_output = cl_resample.Centerlines

	# Compute  Centerline Geometry for Tangents
	cl_vmtk_geom = vmtkscripts.vmtkCenterlineGeometry()
	cl_vmtk_geom.Centerlines = cl_vmtk_output
	cl_vmtk_geom.Execute()
	cl_vmtk_output = cl_vmtk_geom.Centerlines
	cl_vmtk_output_filename = os.path.join(data_path, model_name.split('/')[-1][:-4] + '_clmerge.vtp')
	write_polydata(cl_vmtk_output, cl_vmtk_output_filename)

	return(cl_vmtk_output_filename)


def volume_qoi(model_name, cl_vmtk_file, cl_sv_file):
	"""Create volume clips of the regions of interest using the merged centerlines from VMTK and extract qoi"""
	model = read_polydata(model_name)
	cl_vmtk = read_polydata(cl_vmtk_file)

	# Resample SimVascular centerlines because not refined enough - uses splines
	cl_sv = read_polydata(cl_sv_file)
	cl_resample = vmtkscripts.vmtkCenterlineResampling()
	cl_resample.Centerlines = cl_sv
	cl_resample.Length = 0.01
	cl_resample.Execute()
	cl_sv = cl_resample.Centerlines

	save_path = os.path.join('/'.join(model_name.split('/')[:-1]), model_name.split('/')[-1][:-4])
	
	# Clip Volumes for MPA, LPA, RPA and integrate
	results = volume_clipping(model, cl_vmtk, cl_sv, save_path)
	return results








"""PRIVATE FUNCTIONS"""
def calculate_group_roi(centerlines_output, num_group, group_elems, seg_location, cl_type):
	"""Find the center and sphere for each group's region of interest from the produced centerlines. 
	This will serve as the center for slices taken at the specification seg_location
	Args:
		centerlines_output: centerline vtk object
		num_group (int): number of groups (number of segments) along centerline
		group_elems (list): the cell (for vmtk cl) or the points (for sv cl) belonging to each group
		seg_location (float): location from 0 to 1 at which to take the RoI for each group
		cl_type (str): 'sv' or 'vmtk' to specify what type of centerline
	Returns:
		group_roi_center (list): coordinate points of the center of the slice
		group_roi_tan (list): tangent vector along centerline at location of slice
		group_roi_maxR (list): maximum inscribed radius at that slice
	"""
	pointdata = centerlines_output.GetPointData().GetArray("MaximumInscribedSphereRadius")
	points_maxR = nps.vtk_to_numpy(pointdata)

	#Calculate Frenet tangential vector
	if cl_type.lower() == 'vmtk':
		pointdata = centerlines_output.GetPointData().GetArray("FrenetTangent")
	elif cl_type.lower() == 'sv':
		pointdata = centerlines_output.GetPointData().GetArray("CenterlineSectionNormal")
		num_pts = centerlines_output.GetNumberOfPoints()
		group_list = nps.vtk_to_numpy(centerlines_output.GetPointData().GetArray('BranchId'))
		path_group = nps.vtk_to_numpy(centerlines_output.GetPointData().GetArray('Path'))
	points_tangent = nps.vtk_to_numpy(pointdata)

	#For each group, the center of roi is recorded
	group_roi_center = [0]*num_group
	group_roi_maxR = [0]*num_group
	group_roi_tan = [0]*num_group
	for i in range(0, num_group):
		if cl_type.lower() == 'vmtk':
			ids = vtk.vtkIdList()
			if not group_elems[i]:
				continue
			centerlines_output.GetCellPoints(group_elems[i][0],ids)
			num_ids = ids.GetNumberOfIds()
		elif cl_type.lower() == 'sv':
			ids = [j for j in range(0, num_pts) if group_list[j] == group_elems[i]]
			path_sort = path_group[ids]
			# Sort ids based on a previoulys calculated path since resampled ids out of order
			ids = [x for _,x in sorted(zip(path_sort, ids))]
			num_ids = len(ids)
		print("group = %d, num of points = %d" % (i, num_ids))

		# Calculate distance along path
		path_dist = np.zeros(num_ids)
		for pathptid in range(1, num_ids):
			if cl_type.lower() == 'vmtk':
				id1 = ids.GetId(pathptid)
				id2 = ids.GetId(pathptid-1)
			elif cl_type.lower() == 'sv':
				id1 = ids[pathptid - 1]
				id2 = ids[pathptid]
			pt1 = np.array(centerlines_output.GetPoints().GetPoint(id1))
			pt2 = np.array(centerlines_output.GetPoints().GetPoint(id2))
			pt_dist = np.linalg.norm(pt1-pt2)
			if pt_dist > 0.01*10: #based on cl resampling length
				path_dist[pathptid] = path_dist[pathptid-1]
			else:
				path_dist[pathptid] = path_dist[pathptid-1] + np.linalg.norm(pt1-pt2)

		# Calculate RoI based on point closest to seg_location
		dist2target = abs(path_dist-path_dist[num_ids - 1]*seg_location)
		index = np.argmin(dist2target)
		if cl_type.lower() == 'vmtk':
			tmpid = ids.GetId(index)
		elif cl_type.lower() == 'sv':
			tmpid = ids[np.argmin(dist2target)]
		tmpR = points_maxR[tmpid]
		pt1 = np.array(centerlines_output.GetPoints().GetPoint(tmpid))

		# Append the center and the max inscribed radius of that slice for each group
		group_roi_center[i] = pt1
		group_roi_maxR[i] = tmpR
		# Calculates a tangent normal if the tangent stored in the centerlines is an artifact
		if np.linalg.norm(points_tangent[tmpid]) < 1e-6:
			if index < num_ids-1:
				pt2 = np.array(centerlines_output.GetPoints().GetPoint(ids.GetId(index+1)))
			else:
				pt2 = np.array(centerlines_output.GetPoints().GetPoint(ids.GetId(index-1)))
			dx = np.linalg.norm(pt2-pt1)
			tmptan = (pt2-pt1)/dx
			print("tangent finite diff",tmptan)
			group_roi_tan[i] = tmptan
		else:
			group_roi_tan[i] = points_tangent[tmpid]

	return group_roi_center, group_roi_tan, group_roi_maxR


def centroid(infile):
	"""Calculate the centroid of a vtp file"""
	poly_data = read_polydata(infile)
	x_list = []
	y_list = []
	z_list = []
	for i in range(poly_data.GetNumberOfPoints()):
		x_list.append(float(poly_data.GetPoints().GetPoint(i)[0]))
		y_list.append(float(poly_data.GetPoints().GetPoint(i)[1]))
		z_list.append(float(poly_data.GetPoints().GetPoint(i)[2]))

	return [np.mean(x_list), np.mean(y_list), np.mean(z_list)]


def cut_integrate_slice(ModelName, num_group, group_roi_center, group_roi_tan, group_roi_maxR, seg_location, cl_type):
	"""Create a slice through volume data and integrate quantities of interest
	Args:
		ModelName (str): Full path to the volume data that contains the velocity and is confined anatomically
		num_group (int): number of segments in the volume data (based on centerlines)
		group_roi_center (list): center coordinates of the slice for each group
		group_roi_tan (list): tangent vector of the centerline at the slice for each group
		group_roi_maxR (list): maximum inscribed radius from the centerline to the volume surface for each slice
		seg_location (float): Float form 0 to 1 to create slice down length of vessel segments
		cl_type (str): 'sv' or 'vmtk'
	Returns: Integrated and spatially averaged values of the following for each slice in a group
		group_flow
		group_area
		group_vort
		group_qcrit
	"""
	# Initialize results arrays
	group_qthresh = [[]]*num_group
	group_sepbubble = [[]]*num_group
	group_uorth = [[]]*num_group
	group_relhd = [[]]*num_group
	group_hd = [[]]*num_group
	group_flow = [[]]*num_group
	group_area = [[]]*num_group
	group_vort = [[]]*num_group
	group_qcrit = [[]]*num_group
	group_maxvel = [[]]*num_group

	# Calling for each output
	filename = ModelName 
	vtu_data = read_polydata(filename)

	# Cut a slice and integrate for each group
	for j in range(0, num_group):
		# Determine coordinates and angle of the plane at which to take a slice
		plane = vtk.vtkPlane()
		plane.SetOrigin(group_roi_center[j][0], group_roi_center[j][1], group_roi_center[j][2])
		plane.SetNormal(group_roi_tan[j][0], group_roi_tan[j][1], group_roi_tan[j][2])

		# Perform slice through entire volume data
		cutter = vtk.vtkCutter()
		cutter.SetCutFunction(plane)
		cutter.SetInputData(vtu_data)
		cutter.Update()
		vtu_slice = cutter.GetOutput()
		
		# Confine slice to solely the maximum inscribed radius plus a little extra
		sphere = vtk.vtkSphereSource()
		sphere.SetCenter(group_roi_center[j][0], group_roi_center[j][1], group_roi_center[j][2])
		sphere.SetRadius(group_roi_maxR[j]*3.0)
		sphere.Update()

		implicitPolyDataDistance = vtk.vtkImplicitPolyDataDistance()
		implicitPolyDataDistance.SetInput(sphere.GetOutput())

		# Create additional arrays of information
		signedDistances = vtk.vtkFloatArray()
		signedDistances.SetNumberOfComponents(1)
		signedDistances.SetName("SignedDistances")
		# print"cutter num points",vtu_slice.GetNumberOfPoints()
		if vtu_slice.GetNumberOfPoints() == 0:
			print("0 sliced point for group = ",j)
			exit()
		for pointId in range(vtu_slice.GetNumberOfPoints()):
			tmp_point = vtu_slice.GetPoint(pointId)
			signedDistance = implicitPolyDataDistance.EvaluateFunction(tmp_point)
			signedDistances.InsertNextValue(signedDistance)

		# Compute cross product of velocity and slice normal (to obtain if velocity vector is perpendicular to forward direction)
		velocity_array = vtu_slice.GetPointData().GetArray("velocity")
		uxn_array = nps.numpy_to_vtk(np.dot(velocity_array, group_roi_tan[0]))
		uxn_array.SetName("u.n")
		vtu_slice.GetPointData().AddArray(uxn_array)

		# Clip the slice
		image_array = vtu_slice.GetPointData().GetArray("image")
		vtu_slice.GetPointData().SetScalars(signedDistances) 
		vtu_slice.GetPointData().AddArray(image_array)
		clipper = vtk.vtkClipDataSet()
		clipper.SetInputData(vtu_slice)
		clipper.InsideOutOn()
		clipper.SetValue(0.0)
		clipper.Update()

		# Filter to only slice region in proximity of centerline (sometimes will get slice of the RPA)
		cnnct_filter = vtk.vtkConnectivityFilter()
		cnnct_filter.SetInputData(clipper.GetOutput())
		cnnct_filter.SetExtractionModeToClosestPointRegion()
		cnnct_filter.SetClosestPoint(group_roi_center[j][0], group_roi_center[j][1], group_roi_center[j][2])
		cnnct_filter.Update()
		
		vtu_slice = cnnct_filter.GetOutput()
		print("number of points clip",vtu_slice.GetNumberOfPoints())
		if vtu_slice.GetNumberOfPoints() == 0:
			print("0 clipped point for group=", j)
			exit()

		# Save Slice
		if cl_type == 'sv':
			vtu_slice_file = ModelName[0:-4] + "_group" + str(j) + "_slice" + str(seg_location*10) + "sv.vtu"
		elif cl_type == 'vmtk':
			vtu_slice_file = ModelName[0:-4] + "_group" + str(j) + "_slice" + str(seg_location*10) + "vmtk.vtu"
		write_polydata(vtu_slice, vtu_slice_file)

		# Integrate Slice
		integrated_results = integrate_slice(vtu_slice_file)

		group_qthresh[j] = integrated_results[8]
		group_sepbubble[j] = integrated_results[7]
		group_uorth[j] = integrated_results[6]/integrated_results[0]
		group_relhd[j] = integrated_results[5]
		group_hd[j] = integrated_results[4]
		group_qcrit[j] = integrated_results[3]/integrated_results[0]
		group_vort[j] = integrated_results[2]
		group_flow[j] = integrated_results[1]
		group_area[j] = integrated_results[0]
		group_maxvel[j] = integrated_results[9]

	results = pd.DataFrame()
	results['group_Num'] = range(0, len(group_flow))
	results['Flow'] = group_flow
	results['Area'] = group_area
	results['Vorticity'] = group_vort
	results['Q-criterion'] = group_qcrit
	results['Helicity_Density'] = group_hd
	results['Relative_Helicity_Density'] = group_relhd
	results['u.n'] = group_uorth
	results['Separation_Bubble_Area'] = group_sepbubble
	results['Q-criterion_Threshold_Area'] = group_qthresh
	results['max_Vel'] = group_maxvel

	return results


def estimator_callback(status_dict):
	iteration_status_dictionaries.append(status_dict)
	return True


def read_centerlines(centerlines_output, cl_type):
	"""Read the cell data from the centerline output
	Args:
		centerlines_output: vtk polydata object of the centerlines (read in previously)
		cl_type: "sv" or "vmtk" to specify the data arrays
	"""
	num_pts = centerlines_output.GetNumberOfPoints()
	print("Number of Points:", centerlines_output.GetNumberOfPoints())
	print("Number of Cell Arrays:", centerlines_output.GetCellData().GetNumberOfArrays())
	print("Number of Point Arrays:", centerlines_output.GetPointData().GetNumberOfArrays())

	# A cell/element (named as LINE) is a segment/line that consists of n points.
	# A centerline consists of m cells, m=number of tract ids, the length of a cell/line is an approximation of a group. 
	# In Cell Data, lines are listed from 0 to m. For each line, the first number is the number of points for this line followed by 1st point to the last point.
	num_cells = centerlines_output.GetNumberOfCells()
	print("Number of Cells:", centerlines_output.GetNumberOfCells())

	# Initiate outputs
	num_group = 0
	num_path = 0
	group_elems = []

	if cl_type.lower() == 'sv':
		centerline_list = nps.vtk_to_numpy(centerlines_output.GetPointData().GetArray('CenterlineId'))
		blank_list = nps.vtk_to_numpy(centerlines_output.GetPointData().GetArray('BifurcationId'))
		group_list = nps.vtk_to_numpy(centerlines_output.GetPointData().GetArray('BranchId'))

		num_path = centerline_list[-1] + 1
		print('number of paths = ', num_path)
		num_group = len(set(group_list))
		print('number of groups = ', num_group)

		#group_elems[i] records the unique group numbers
		group_elems = list(set(group_list))

	if cl_type.lower() == 'vmtk':
		# Read cell data, for each cell (line), the lists record its centerline id (n lines starting from the inlet to outlets), blanking (0 non bifurcation, 1 bifurcation), 
		# group ids (vessels are splitted into single segments/branches and a bifurcation region, 
		# if a vessel shared by multiple centerlines, there multiple elements sharing the same group id
		# if a segment is the terminal segment without any child elements, its group id is unique.
		# for a bifurcation region, the elemens' blanking is 1.
		centerline_list = nps.vtk_to_numpy( centerlines_output.GetCellData().GetArray('CenterlineIds'))
		blank_list = nps.vtk_to_numpy(centerlines_output.GetCellData().GetArray('Blanking'))
		group_list = nps.vtk_to_numpy(centerlines_output.GetCellData().GetArray('GroupIds'))
		tract_list = nps.vtk_to_numpy(centerlines_output.GetCellData().GetArray('TractIds'))

		num_path = centerline_list[-1] + 1 
		print("number of paths=", num_path)
		# num_group = max(group_list) + 1
		num_group = len(set(group_list))
		print("number of groups=",num_group)

		#group_elems[i] records the element(line) indices for group id=i
		group_elems = []
		for i in range(0,num_group):
			group_elems.append([])

		for i in range(0,num_cells):
		#  print "cellid=",i
		#  print "groupid=",group_list[i]
			group_elems[i].append(i)
			# group_elems[group_list[i]].append(i)

	return num_group, num_path, num_cells, group_elems


def volume_clipping(model, cl_vmtk, cl_sv, save_path):
	"""Create a volume clip for the MPA, LPA, and RPA based on vmtk and sv centerlines
	Args:
		model: already read in as a vtk object
		cl_vmtk: vmtk merged centerline (vtk object)
		cl_sv: simvascular centerline (vtk object)
	"""
	# Process VMTK centerlines
	vmtk_num_group, num_path, num_cells, group_elems = read_centerlines(cl_vmtk, 'vmtk')
	vmtk_group_roi_center, vmtk_group_roi_tan, vmtk_group_roi_maxR = calculate_group_roi(cl_vmtk, vmtk_num_group, group_elems, 1.0, 'vmtk')

	# Process SV Centerlines
	sv_num_group, num_path, num_cells, group_elems = read_centerlines(cl_sv, 'sv')
	sv_group_roi_center, sv_group_roi_tan, sv_group_roi_maxR = calculate_group_roi(cl_sv, sv_num_group, group_elems, 0.1, 'sv')

	# Initialize results arrays
	group_qthresh = [[]]*vmtk_num_group
	group_poshd = [[]]*vmtk_num_group
	group_neghd = [[]]*vmtk_num_group
	group_relhd = [[]]*vmtk_num_group
	group_hd = [[]]*vmtk_num_group
	group_vol = [[]]*vmtk_num_group
	group_vort = [[]]*vmtk_num_group
	group_qcrit = [[]]*vmtk_num_group
	group_ke = [[]]*vmtk_num_group
	group_el = [[]]*vmtk_num_group
	group_qvthresh = [[]]*vmtk_num_group
	group_qvthresh_11 = [[]]*vmtk_num_group
	group_qvthresh_12 = [[]]*vmtk_num_group
	group_qvthresh_13 = [[]]*vmtk_num_group
	group_qvthresh_14 = [[]]*vmtk_num_group
	group_qvthresh_21 = [[]]*vmtk_num_group
	group_qvthresh_22 = [[]]*vmtk_num_group
	group_qvthresh_23 = [[]]*vmtk_num_group
	group_qvthresh_24 = [[]]*vmtk_num_group
	group_absrelhd = [[]]*vmtk_num_group

	# Create Volume Clipping
	# Use vtkconnectivity filter to ensure that volume clips do not overlap between MPA/LPA/RPA
	plane = vtk.vtkPlane()
	extractor = vtk.vtkExtractGeometry()
	extractor.SetInputData(model)
	for i in range(0, vmtk_num_group):
		# (OLD CLIPPING) For MPA (equiv to group 0), take a volume clip from the last slice of group0 from the VMTK centerlines towards anterior
		# NEW CLIPPING: For MPA (equiv to gorup 0), take all volumes of PA that is not included in group1 and group2 clipping (assumes only MPA, lobar LPA/RPA included in model)
		if i == 0:
			# Clip everything not in LPA
			plane.SetOrigin(sv_group_roi_center[1][0], sv_group_roi_center[1][1], sv_group_roi_center[1][2])
			plane.SetNormal(-sv_group_roi_tan[1][0], -sv_group_roi_tan[1][1], -sv_group_roi_tan[1][2])
			extractor.SetExtractInside(False)
			extractor.SetImplicitFunction(plane)
			extractor.SetExtractBoundaryCells(True)
			extractor.Update()

			# Clip everything not in RPA AND not in LPA
			mpa_extractor = vtk.vtkExtractGeometry()
			mpa_extractor.SetInputData(extractor.GetOutput())
			plane.SetOrigin(sv_group_roi_center[2][0], sv_group_roi_center[2][1], sv_group_roi_center[2][2])
			plane.SetNormal(-sv_group_roi_tan[2][0], -sv_group_roi_tan[2][1], -sv_group_roi_tan[2][2])
			mpa_extractor.SetExtractInside(False)
			mpa_extractor.SetImplicitFunction(plane)
			mpa_extractor.SetExtractBoundaryCells(True)
			mpa_extractor.Update()

			cnnct_filter = vtk.vtkConnectivityFilter()
			cnnct_filter.SetInputData(mpa_extractor.GetOutput())
			cnnct_filter.SetExtractionModeToSpecifiedRegions()
			cnnct_filter.AddSpecifiedRegion(0)
			cnnct_filter.Update()

			clipped_vol = cnnct_filter.GetOutput()
		
		# For LPA/RPA (equiv to group 1/2), take a volume clip from the first slice of group1 and 2 of the SV centerlines towards the end points centerlines
		else:
			plane.SetOrigin(sv_group_roi_center[i][0], sv_group_roi_center[i][1], sv_group_roi_center[i][2])
			plane.SetNormal(sv_group_roi_tan[i][0], sv_group_roi_tan[i][1], sv_group_roi_tan[i][2])

			extractor.SetExtractInside(False)
			extractor.SetImplicitFunction(plane)
			extractor.SetExtractBoundaryCells(True)
			extractor.Update()

			cnnct_filter = vtk.vtkConnectivityFilter()
			cnnct_filter.SetInputData(extractor.GetOutput())
			cnnct_filter.SetExtractionModeToClosestPointRegion()
			cnnct_filter.SetClosestPoint(sv_group_roi_center[i][0], sv_group_roi_center[i][1], sv_group_roi_center[i][2])
			cnnct_filter.Update()

			clipped_vol = cnnct_filter.GetOutput()
		
		clip_file = save_path + '_' + str(i) + '_clip.vtu'
		write_polydata(clipped_vol, clip_file)

		# Integrate volume
		start_intvol = time.time()
		integrated_results = integrate_volume(clipped_vol)
		end_intvol = time.time()
		print("Time Elapsed to integrate volume: %d" % (end_intvol - start_intvol))
		group_qthresh[i] = integrated_results[5]
		group_poshd[i] = integrated_results[7]
		group_neghd[i] = integrated_results[6]
		group_relhd[i] = integrated_results[4]
		group_hd[i] = integrated_results[3]
		group_qcrit[i] = integrated_results[2]/integrated_results[0]
		group_vort[i] = integrated_results[1]
		group_vol[i] = integrated_results[0]
		group_ke[i] = integrated_results[8]
		group_el[i] = integrated_results[9]
		group_qvthresh[i] = integrated_results[10]
		group_qvthresh_11[i] = integrated_results[11]
		group_qvthresh_12[i] = integrated_results[12]
		group_qvthresh_13[i] = integrated_results[13]
		group_qvthresh_14[i] = integrated_results[14]
		group_qvthresh_21[i] = integrated_results[15]
		group_qvthresh_22[i] = integrated_results[16]
		group_qvthresh_23[i] = integrated_results[17]
		group_qvthresh_24[i] = integrated_results[18]
		group_absrelhd[i] = integrated_results[19]

	results = pd.DataFrame()
	results['group_Num'] = range(0, len(group_vol))
	results['Volume'] = group_vol
	results['Vorticity'] = group_vort
	results['Q-criterion'] = group_qcrit
	results['Helicity_Density'] = group_hd
	results['Relative_Helicity_Density'] = group_relhd
	results['Q-criterion_Threshold_Vol'] = group_qthresh
	results['Positive_Helicity_Vol'] = group_poshd
	results['Negative_Helicity_Vol'] = group_neghd
	results['Kinetic_Energy'] = group_ke
	results['Rate_Viscous_Energy_Loss'] = group_el
	results['Q-criterion_Velocity_Threshold_Vol'] = group_qvthresh
	results['Q-criterion_Velocity_Threshold_11_Vol'] = group_qvthresh_11
	results['Q-criterion_Velocity_Threshold_12_Vol'] = group_qvthresh_12
	results['Q-criterion_Velocity_Threshold_13_Vol'] = group_qvthresh_13
	results['Q-criterion_Velocity_Threshold_14_Vol'] = group_qvthresh_14
	results['Q-criterion_Velocity_Threshold_21_Vol'] = group_qvthresh_21
	results['Q-criterion_Velocity_Threshold_22_Vol'] = group_qvthresh_22
	results['Q-criterion_Velocity_Threshold_23_Vol'] = group_qvthresh_23
	results['Q-criterion_Velocity_Threshold_24_Vol'] = group_qvthresh_24
	results['Absolute_Relative_Helicity_Density'] = group_absrelhd

	return results


def integrate_volume(clip):
	"""Integrate variables in volume clip (one per time and branch) and extract quantities of interest with thresholding
		Integrates by averaging point data for each cell
	Args:
		vol_clip: vtu volume clip of MPA/LPA/RPA
	Return: Dataframe that includes the following quantities
		Volume
		Vorticity (spatially averaged from volume)
		Q-criterion
		Helicity
		Relative Helicity Density
		Q-criterion thresholded volume - volume with Q-criterion > 0
		Positive Helicity thresholded volume - volume with Helicity > 0
		Negative Helicity Thresholded volume - volume with Helicity < 0
	"""
	velocity_array = clip.GetPointData().GetArray("velocity")
	vorticity_array = clip.GetPointData().GetArray("vorticity")
	qcrit_array = clip.GetPointData().GetArray("Q-criterion")
	hel_dens_array = clip.GetPointData().GetArray("helicity_density")
	rel_hel_dens_array = clip.GetPointData().GetArray("relative_helicity_density")
	visc_diss_array = clip.GetPointData().GetArray("viscous_dissipation")
	q_crit_thresh_vol = 0.0
	q_crit_vel_thresh_vol = 0.0
	q_crit_vel_thresh_11_vol = 0.0
	q_crit_vel_thresh_12_vol = 0.0
	q_crit_vel_thresh_13_vol = 0.0
	q_crit_vel_thresh_14_vol = 0.0
	q_crit_vel_thresh_21_vol = 0.0
	q_crit_vel_thresh_22_vol = 0.0
	q_crit_vel_thresh_23_vol = 0.0
	q_crit_vel_thresh_24_vol = 0.0
	pos_hd_vol = 0.0
	neg_hd_vol = 0.0
	clip_vort_avg = 0.0
	clip_q_avg = 0.0
	clip_hd_avg = 0.0
	clip_relhd_avg = 0.0
	clip_absrelhd_avg = 0.0
	clip_ke = 0.0 #kinetic energy over a volumetric domain of interest
	clip_el = 0.0 #rate of viscous energy loss (power loss) over a volumetric domain of interest

	num_cells = clip.GetNumberOfCells()
	missing_cell = []
	clip_vol = 0.0
	for i in range(num_cells):
		cell = clip.GetCell(i)
		p0 = cell.GetPoints().GetPoint(0)
		p1 = cell.GetPoints().GetPoint(1)
		p2 = cell.GetPoints().GetPoint(2)
		# Clip may not go incorporate full cells. If only partial tetrahedron, ignore
		try:
			p3 = cell.GetPoints().GetPoint(3)
		except:
			# print('ERROR: Point 3 not found in cell %2d' % i)
			missing_cell.append(i)
			continue
		cell_vol = vtk.vtkMeshQuality.TetVolume(cell)
		clip_vol += cell_vol
		nodeids = cell.GetPointIds()
		# Compute Vorticity Magnitude
		w1 = np.array(vorticity_array.GetTuple3(nodeids.GetId(0)))
		w2 = np.array(vorticity_array.GetTuple3(nodeids.GetId(1)))
		w3 = np.array(vorticity_array.GetTuple3(nodeids.GetId(2)))
		w4 = np.array(vorticity_array.GetTuple3(nodeids.GetId(3)))
		wx_mean = np.mean([w1[0], w2[0], w3[0], w4[0]])
		wy_mean = np.mean([w1[1], w2[1], w3[1], w4[1]])
		wz_mean = np.mean([w1[2], w2[2], w3[2], w4[2]])
		clip_vort_avg = clip_vort_avg + np.linalg.norm([wx_mean, wy_mean, wz_mean])*cell_vol
		# Compute Q-criterion
		q1 = np.array(qcrit_array.GetTuple(nodeids.GetId(0)))
		q2 = np.array(qcrit_array.GetTuple(nodeids.GetId(1)))
		q3 = np.array(qcrit_array.GetTuple(nodeids.GetId(2)))
		q4 = np.array(qcrit_array.GetTuple(nodeids.GetId(3)))
		clip_q_avg = clip_q_avg + np.mean([q1, q2, q3, q4])*cell_vol
		# Compute Helicity Density
		hd1 = np.array(hel_dens_array.GetTuple(nodeids.GetId(0)))
		hd2 = np.array(hel_dens_array.GetTuple(nodeids.GetId(1)))
		hd3 = np.array(hel_dens_array.GetTuple(nodeids.GetId(2)))
		hd4 = np.array(hel_dens_array.GetTuple(nodeids.GetId(3)))
		clip_hd_avg = clip_hd_avg + np.mean([hd1, hd2, hd3, hd4])*cell_vol
		# Compute Relative Helicity density
		relhd1 = np.array(rel_hel_dens_array.GetTuple(nodeids.GetId(0)))
		relhd2 = np.array(rel_hel_dens_array.GetTuple(nodeids.GetId(1)))
		relhd3 = np.array(rel_hel_dens_array.GetTuple(nodeids.GetId(2)))
		relhd4 = np.array(rel_hel_dens_array.GetTuple(nodeids.GetId(3)))
		clip_relhd_avg = clip_relhd_avg + np.mean([relhd1, relhd2, relhd3, relhd4])*cell_vol
		# Absolute Relative Helicity Density
		clip_absrelhd_avg = clip_absrelhd_avg + abs(np.mean([relhd1, relhd2, relhd3, relhd4]))*cell_vol
		# Kinetic Energy
		density = 1.06 # g/cm3
		v1 = np.array(velocity_array.GetTuple3(nodeids.GetId(0)))
		v2 = np.array(velocity_array.GetTuple3(nodeids.GetId(1)))
		v3 = np.array(velocity_array.GetTuple3(nodeids.GetId(2)))
		v4 = np.array(velocity_array.GetTuple3(nodeids.GetId(3)))
		ke1 = np.linalg.norm(v1)**2
		ke2 = np.linalg.norm(v2)**2
		ke3 = np.linalg.norm(v3)**2
		ke4 = np.linalg.norm(v4)**2
		clip_ke = clip_ke + (density/2)*np.mean([ke1, ke2, ke3, ke4])*cell_vol
		# Rate of Viscous Energy Loss
		viscosity = 0.04 # Poise
		visc1 = np.array(visc_diss_array.GetTuple(nodeids.GetId(0)))
		visc2 = np.array(visc_diss_array.GetTuple(nodeids.GetId(1)))
		visc3 = np.array(visc_diss_array.GetTuple(nodeids.GetId(2)))
		visc4 = np.array(visc_diss_array.GetTuple(nodeids.GetId(3)))
		clip_el = clip_el + viscosity*(np.mean([visc1, visc2, visc3, visc4])*cell_vol)
		# Threshold Q-criterion - if greater than 0, presence of a vortex
		if np.mean([q1, q2, q3, q4]) > 0:
			q_crit_thresh_vol += cell_vol
			# Threshold volume by low velocity to get vorticity that visually looks like a vortex
			v1_mag = np.linalg.norm(v1)
			v2_mag = np.linalg.norm(v2)
			v3_mag = np.linalg.norm(v3)
			v4_mag = np.linalg.norm(v4)
			if np.mean([v1_mag, v2_mag, v3_mag, v4_mag]) < 30:
				q_crit_vel_thresh_vol += cell_vol
			# Other thresholding for vortex volumes
			if np.mean([q1, q2, q3, q4]) > 100:
				if np.mean([v1_mag, v2_mag, v3_mag, v4_mag]) < 30:
					q_crit_vel_thresh_11_vol += cell_vol
				if np.mean([v1_mag, v2_mag, v3_mag, v4_mag]) < 37:
					q_crit_vel_thresh_12_vol += cell_vol
				if np.mean([v1_mag, v2_mag, v3_mag, v4_mag]) < 20:
					q_crit_vel_thresh_13_vol += cell_vol
				if np.mean([v1_mag, v2_mag, v3_mag, v4_mag]) < 55:
					q_crit_vel_thresh_14_vol += cell_vol
			# Other thresholding for vortex volumes
			if np.mean([q1, q2, q3, q4]) > 500:
				if np.mean([v1_mag, v2_mag, v3_mag, v4_mag]) < 30:
					q_crit_vel_thresh_21_vol += cell_vol
				if np.mean([v1_mag, v2_mag, v3_mag, v4_mag]) < 37:
					q_crit_vel_thresh_22_vol += cell_vol
				if np.mean([v1_mag, v2_mag, v3_mag, v4_mag]) < 20:
					q_crit_vel_thresh_23_vol += cell_vol
				if np.mean([v1_mag, v2_mag, v3_mag, v4_mag]) < 55:
					q_crit_vel_thresh_24_vol += cell_vol
		# Threshold Positive Helicity
		if np.mean([hd1, hd2, hd3, hd4]) > 0:
			pos_hd_vol += cell_vol
		elif np.mean([hd1, hd2, hd3, hd4]) < 0:
			neg_hd_vol += cell_vol



	# print("WARNING: The %d out of %d total cells were not included in the volume calculation due to a missing point." % (len(missing_cell), num_cells))
	# print(str(missing_cell))
	integrated_results = [clip_vol, clip_vort_avg, clip_q_avg, clip_hd_avg, clip_relhd_avg, 
						  q_crit_thresh_vol, pos_hd_vol, neg_hd_vol, clip_ke, clip_el, q_crit_vel_thresh_vol, q_crit_vel_thresh_11_vol, q_crit_vel_thresh_12_vol, 
						  q_crit_vel_thresh_13_vol, q_crit_vel_thresh_14_vol, q_crit_vel_thresh_21_vol, q_crit_vel_thresh_22_vol, q_crit_vel_thresh_23_vol, 
						  q_crit_vel_thresh_24_vol, clip_absrelhd_avg]

	return integrated_results


def cl_smallcut_integrate_slice(ModelName, num_group, group_roi_center, group_roi_tan, group_roi_maxR, seg_location, cl_type):
	"""Cut a small center proportion from previously created slice to approximate centerline velocity (30% of original slice)
	Args:
		ModelName (str): Full path to the volume data that contains the velocity and is confined anatomically
		num_group (int): number of segments in the volume data (based on centerlines)
		group_roi_center (list): center coordinates of the slice for each group
		group_roi_tan (list): tangent vector of the centerline at the slice for each group
		group_roi_maxR (list): maximum inscribed radius from the centerline to the volume surface for each slice
		seg_location (float): Float form 0 to 1 to create slice down length of vessel segments
		cl_type (str): 'sv' or 'vmtk'
	Returns: Integrated and spatially averaged values of the following for each slice in a group
		group_avgvel
		group_maxvel
	"""
	# Initialize results arrays
	group_avgvel = [[]]*num_group
	group_maxvel = [[]]*num_group

	# Cut a slice and integrate for each group
	for j in range(0, num_group):
		# Calling for each output
		if cl_type == 'sv':
			vtu_slice_file = ModelName[0:-4] + "_group" + str(j) + "_slice" + str(seg_location*10) + "sv.vtu"
		elif cl_type == 'vmtk':
			vtu_slice_file = ModelName[0:-4] + "_group" + str(j) + "_slice" + str(seg_location*10) + "vmtk.vtu"
		vtu_slice = read_polydata(vtu_slice_file)
		
		# Confine slice to solely the maximum inscribed radius plus a little extra
		sphere = vtk.vtkSphereSource()
		sphere.SetCenter(group_roi_center[j][0], group_roi_center[j][1], group_roi_center[j][2])
		sphere.SetRadius(group_roi_maxR[j]*0.3)
		sphere.Update()

		implicitPolyDataDistance = vtk.vtkImplicitPolyDataDistance()
		implicitPolyDataDistance.SetInput(sphere.GetOutput())

		# Create additional arrays of information
		signedDistances = vtk.vtkFloatArray()
		signedDistances.SetNumberOfComponents(1)
		signedDistances.SetName("SignedDistances")
		# print"cutter num points",vtu_slice.GetNumberOfPoints()
		if vtu_slice.GetNumberOfPoints() == 0:
			print("0 sliced point for group = ",j)
			exit()
		for pointId in range(vtu_slice.GetNumberOfPoints()):
			tmp_point = vtu_slice.GetPoint(pointId)
			signedDistance = implicitPolyDataDistance.EvaluateFunction(tmp_point)
			signedDistances.InsertNextValue(signedDistance)

		# Compute cross product of velocity and slice normal (to obtain if velocity vector is perpendicular to forward direction)
		velocity_array = vtu_slice.GetPointData().GetArray("velocity")
		uxn_array = nps.numpy_to_vtk(np.dot(velocity_array, group_roi_tan[0]))
		uxn_array.SetName("u.n")
		vtu_slice.GetPointData().AddArray(uxn_array)

		# Clip the slice
		image_array = vtu_slice.GetPointData().GetArray("image")
		vtu_slice.GetPointData().SetScalars(signedDistances) 
		vtu_slice.GetPointData().AddArray(image_array)
		clipper = vtk.vtkClipDataSet()
		clipper.SetInputData(vtu_slice)
		clipper.InsideOutOn()
		clipper.SetValue(0.0)
		clipper.Update()

		# Filter to only slice region in proximity of centerline (sometimes will get slice of the RPA)
		cnnct_filter = vtk.vtkConnectivityFilter()
		cnnct_filter.SetInputData(clipper.GetOutput())
		cnnct_filter.SetExtractionModeToClosestPointRegion()
		cnnct_filter.SetClosestPoint(group_roi_center[j][0], group_roi_center[j][1], group_roi_center[j][2])
		cnnct_filter.Update()
		
		vtu_slice = cnnct_filter.GetOutput()
		print("number of points clip",vtu_slice.GetNumberOfPoints())
		if vtu_slice.GetNumberOfPoints() == 0:
			print("0 clipped point for group=", j)
			exit()

		# Integrate quantities of interest over slice
		vn = [0.0,0.0,0.0]
		prev_vn = vn
		slice_area = 0.0
		slice_flow = 0.0
		max_vel = 0.0
		velocity_array = vtu_slice.GetPointData().GetArray("velocity")

		missing_cell = []
		for cellId in range(vtu_slice.GetNumberOfCells()):
			cell = vtu_slice.GetCell(cellId)
			p0 = cell.GetPoints().GetPoint(0)
			p1 = cell.GetPoints().GetPoint(1)
			# For 4DMRI data, if slice goes through centers slices, it will only include full cells with 3 points
			try:
				p2 = cell.GetPoints().GetPoint(2) 
			except:
				# print('ERROR: Point 2 not found in cell %2d' % cellId)
				missing_cell.append(cellId)
				continue
			cell_area = vtk.vtkTriangle().TriangleArea(p0, p1, p2)
			slice_area = slice_area + cell_area
			nodeids = cell.GetPointIds()
			# Compute Flow
			v1 = np.array(velocity_array.GetTuple3(nodeids.GetId(0)))
			v2 = np.array(velocity_array.GetTuple3(nodeids.GetId(1)))
			v3 = np.array(velocity_array.GetTuple3(nodeids.GetId(2)))
			vtk.vtkTriangle().ComputeNormal(p0, p1, p2, vn)
			if prev_vn != vn:
				print("Previous cell does not have same normal as current cell: " + str(vn))
			prev_vn = vn
			slice_flow = slice_flow + sum((v1+v2+v3)/3.0*vn)*cell_area # 1 point Gauss quad rule
			v1_mag = np.linalg.norm(v1)
			v2_mag = np.linalg.norm(v2)
			v3_mag = np.linalg.norm(v3)
			if np.mean([v1_mag, v2_mag, v3_mag]) > max_vel:
				max_vel =  np.mean([v1_mag, v2_mag, v3_mag])

		group_avgvel[j] = slice_flow/slice_area
		group_maxvel[j] = max_vel

	results = pd.DataFrame()
	results['group_Num'] = range(0, len(group_maxvel))
	results['cl_avg_Vel'] = group_avgvel
	results['cl_max_Vel'] = group_maxvel

	return results