"""Shape registration between systolic and diastolic models using Deformetrica

Requires the following additional python packages:
- deformetrica

=================================================================================================================================
---Created--- | ------Owner------- | Notes---------------------------------------------------------------------------------------
=================================================================================================================================
   10-12-2020   Melody Dong          Created
=================================================================================================================================
"""

import os
from os import listdir, path
from os.path import isfile, join
import sys

import deformetrica as dfca
import itertools
import multiprocessing as mp
import pandas as pd
import pdb
import time
import vtk
from vmtk import vtkvmtk, vmtkscripts
import warnings
warnings.filterwarnings("ignore")

import post4DMRI_helper as post

do_registration = False
do_volumeremeshing = False
do_resample4dmri = False
do_qoiextraction = False
centerline_already_extracted = False
meshing_type = 'vmtk'
data_path = ''
subject_id = ''
mr_vtk_dir = ''
mr_subject_id = ''
dia_tmpl = ''
sys_tmpl = ''
max_iterations = 0
total_timepoints = 0
remesh_edge_length = 0.0
seg_locations = []
procs = 1

def main(input_file):

	# Define input variables
	read_input_file(input_file)
	# Assumes templates files are named in the following method: 'mag#.vtp' with 1 being the first timepoint
	if int(dia_tmpl[3:-4]) < int(sys_tmpl[3:-4]):
		sys_upstroke_timepoints = int(sys_tmpl[3:-4]) - int(dia_tmpl[3:-4]) + 1
		sys_downstroke_timepoints = total_timepoints - sys_upstroke_timepoints + 1
	else:
		sys_downstroke_timepoints = int(dia_tmpl[3:-4]) - int(sys_tmpl[3:-4]) + 1
		sys_upstroke_timepoints = total_timepoints - sys_downstroke_timepoints + 1
		
	# For Centerline Extraction
	cl_path = os.path.join(data_path, 'cl_sv')

	#################################################################################################################################
	# Shape Analysis - Registration
	if do_registration:
		print(post.registration('pa_sys', subject_id, os.path.join(data_path, dia_tmpl), 
						   os.path.join(data_path, sys_tmpl), max_iterations, sys_upstroke_timepoints, 1.0))
		print(post.registration('pa_dia', subject_id, os.path.join(data_path, sys_tmpl), 
						   os.path.join(data_path, dia_tmpl), max_iterations, sys_downstroke_timepoints, 0.8))


	# Convert deformed mesh into volume meshes
	if do_volumeremeshing:
		start_mesh = time.time()
		reg_vtk_err = []
		geom_files = []
		output_files = []
		time_pt = -1
		for i, file in enumerate(sorted(listdir(os.path.join(data_path, 'registration_output')))):
			isfile_check = isfile(join(os.path.join(data_path, 'registration_output'), file))
			sys_stroke_check = 'pa_sys' in file or 'pa_dia' in file
			vtk_check = '.vtk' in file
			tp_check = '__tp_' in file
			if isfile_check and sys_stroke_check and vtk_check and tp_check:
				# Identify which timestep the geometry belongs to
				if 'pa_sys' in file:
					def_timepoint = int(file[int(file.find('__tp_'))+5:-4]) -1 + int(dia_tmpl[3:-4])
					if def_timepoint == sys_upstroke_timepoints - 1 + int(dia_tmpl[3:-4]) - 1:
						continue
				if 'pa_dia' in file:
					reg_timepoint = int(file[int(file.find('__tp_'))+5:-4])
					def_timepoint = sys_upstroke_timepoints + reg_timepoint -1 + int(dia_tmpl[3:-4]) -1
				if def_timepoint >= total_timepoints:
					def_timepoint -= total_timepoints
				output_file = '%s_%02d.vtu' % (subject_id, def_timepoint)
				# Check that the geometry has not already been meshed
				if isfile(os.path.join(data_path, 'registration_output', output_file)):
					continue
				geom_files.append(file)
				output_files.append(output_file)

		# Parallel Volume remeshing of geometry
		if int(mp.cpu_count()) < procs:
			num_processors = mp.cpu_count()
		else:
			num_processors = procs
		pool = mp.Pool(num_processors)
		output = pool.starmap(post.volume_remesh, [(geom_files[i], remesh_edge_length, output_files[i], data_path, meshing_type) 
												   for i in range(0, len(geom_files))])
		err_ind = [i for i, x in enumerate(output) if x=='ERROR']
		reg_vtk_err = [geom_files[i] for i in err_ind]

		# While loop for remeshing to catch all meshes that failed previously
		# TODO: Find bug that causes meshing to fail randomly
		n = 10
		iter_count = 0
		while len(reg_vtk_err) > 0 and iter_count <= n:
			iter_count += 1
			tmp_reg_vtk_err = []
			geom_files = []
			output_files = []
			for file in reg_vtk_err:
				isfile_check = isfile(join(os.path.join(data_path, 'registration_output'), file))
				sys_stroke_check = 'pa_sys' in file or 'pa_dia' in file
				vtk_check = '.vtk' in file
				tp_check = '__tp_' in file
				if isfile_check and sys_stroke_check and vtk_check and tp_check:
					if 'pa_sys' in file:
						def_timepoint = int(file[int(file.find('__tp_'))+5:-4])
						if def_timepoint == sys_upstroke_timepoints-1:
							continue
					if 'pa_dia' in file:
						reg_timepoint = int(file[int(file.find('__tp_'))+5:-4])
						def_timepoint = sys_upstroke_timepoints + reg_timepoint -1
					print('Attempting to remesh ' + file + '\t' + str(def_timepoint))
					output_file = '%s_%02d.vtu' % (subject_id, def_timepoint)
					# Check that the geometry has not already been meshed
					if isfile(os.path.join(data_path, 'registration_output', output_file)):
						continue
					geom_files.append(file)
					output_files.append(output_file)
			# Parallel Volume resmeshing
			output = pool.starmap(post.volume_remesh, [(geom_files[i], remesh_edge_length, output_files[i], data_path, meshing_type) 
													   for i in range(0, len(geom_files))])
			err_ind = [i for i, x in enumerate(output) if x=='ERROR']
			reg_vtk_err = [geom_files[i] for i in err_ind]
		pool.close()
		end_mesh = time.time()
		print("Volume remeshing: %d" % (end_mesh - start_mesh))


	# Resample 4DMRI with MPA mesh
	if do_resample4dmri:
		start_resample = time.time()
		mpa_meshes = []
		for file in listdir(os.path.join(data_path, 'registration_output')):
			if 'vtu' in file and subject_id in file:
				mpa_meshes.append(file)
		mpa_meshes = sorted(mpa_meshes)
		mr_files = []
		mesh_files = []
		for i in range(0, total_timepoints):
			file_str = '%s_%02d.vtk' % (mr_subject_id, i)
			if not (os.path.exists(os.path.join(mr_vtk_dir, 'RESULTS'))):
				os.mkdir(os.path.join(mr_vtk_dir, 'RESULTS'))
			mr_files.append(os.path.join(mr_vtk_dir, file_str))
			mesh_files.append(os.path.join(data_path, 'registration_output', mpa_meshes[int(i)]))
		# Parallelize
		if int(mp.cpu_count()) < procs:
			num_processors = mp.cpu_count()
		else:
			num_processors = procs
		pool = mp.Pool(num_processors)
		output = pool.starmap(post.resample_4dmr, [(mesh_files[i], mr_files[i], os.path.join(mr_vtk_dir, 'RESULTS')) 
												   for i in range(0, len(mr_files))])
		pool.close()
		end_resample = time.time()
		print("Resampled 4DMRI: %d" % (end_resample - start_resample))


	# Extract quantities of interest from 4DMRI data
	if do_qoiextraction:
		start = time.time()
		# Create centerline path to ensure that it exists
		if not(os.path.exists(cl_path)):
			os.mkdir(cl_path)

		# Create lists of filenames and iteratables for analysis
		surf_model_names = []
		resampled_model_names = []
		centerlines_outputs = []
		cl_file_names = []
		for i in range(0, total_timepoints):
			surf_model_names.append(os.path.join(data_path, 'registration_output', '%s_%02d.vtu' % (subject_id, i)))
			resampled_model_names.append(os.path.join(mr_vtk_dir, 'RESULTS', mr_subject_id + '_%02d_resampled.vtu' % i))
			cl_file_names.append(os.path.join(cl_path, 'centerlines-result_%02d.vtp' % i))
			# Centerline Extraction or Reading
			if centerline_already_extracted:
				centerlines_outputs.append(post.read_polydata(cl_file_names[-1]))
			else:
				centerlines_outputs.append(post.centerline_extraction(surf_model_names[-1], cl_file_names[-1]))
		# List of slice and timepoint combinations
		slice_iters = list(itertools.product(list(range(0, total_timepoints)), seg_locations))
		
		if int(mp.cpu_count()) < procs:
			num_processors = mp.cpu_count()
		else:
			num_processors = procs
		pool = mp.Pool(num_processors)
		# Create merged centerlines using VMTK
		start_merge = time.time()
		cl_vmtk_merged_file_names = pool.starmap(post.merged_centerline, [(surf_model_names[i], cl_file_names[i], os.path.join(mr_vtk_dir, 'RESULTS')) 
																		  for i in range(0, total_timepoints)])
		print("CL Merging: %d" % (time.time() - start_merge))
		# Compute additional QoI based on vortical quantities
		start_vort = time.time()
		pool.starmap(post.compute_vortical_quantities, [(x,) for x in resampled_model_names])
		print("Vortical Quantities: %d" % (time.time() - start_vort))
		# Volume clipping of Model based on MPA/LPA/RPA branches and VMTK and SV centerlines
		start_volclip = time.time()
		results_vol = pool.starmap(post.volume_qoi, [(resampled_model_names[i], cl_vmtk_merged_file_names[i], cl_file_names[i]) 
													 for i in range(0, total_timepoints)])
		print("Volume Clipping: %d" % (time.time() - start_volclip))
		# Slicing of Model in the MPA/LPA/RPA with SV and VMTK centerlines
		start_svslice = time.time()
		results_sv = pool.starmap(post.extract_qoi, [(resampled_model_names[i], cl_file_names[i], j, 'sv') for i,j in slice_iters])
		print("SimVascular Slicing: %d" % (time.time() - start_svslice))
		start_vmtkslice = time.time()
		results_vmtk = pool.starmap(post.extract_qoi, [(resampled_model_names[i], cl_vmtk_merged_file_names[i], j, 'vmtk') for i,j in slice_iters])
		print("VMTK Slicing: %d" % (time.time() - start_vmtkslice))
		pool.close()
		pool.join()

		# Summarize results in pandas dataframe
		qoi_vol_4dmri = pd.concat(results_vol, keys=list(range(0, total_timepoints))).reset_index()
		qoi_vol_4dmri.drop(['level_1'], axis=1)
		qoi_vol_4dmri.rename(columns={'level_0': 'Time'})
		slices_sv = pd.concat(results_sv, keys=slice_iters).reset_index()
		slices_sv = slices_sv.drop(['level_2'], axis=1)
		slices_sv = slices_sv.rename(columns={'level_0': 'Time', 'level_1': 'Slice'})
		slices_vmtk = pd.concat(results_vmtk, keys=slice_iters).reset_index()
		slices_vmtk = slices_vmtk.drop(['level_2'], axis=1)
		slices_vmtk = slices_vmtk.rename(columns={'level_0': 'Time', 'level_1': 'Slice'})
		qoi_slice_4dmri = pd.concat([slices_vmtk[slices_vmtk['group_Num']==0], 
									 slices_sv[slices_sv['group_Num']==1],
									 slices_sv[slices_sv['group_Num']==2]])
		qoi_vol_4dmri.to_csv(os.path.join(mr_vtk_dir, 'RESULTS', 'qoi_volclip.csv'))
		qoi_slice_4dmri.to_csv(os.path.join(mr_vtk_dir, 'RESULTS', 'qoi_slice.csv'))

		print("Time Elapsed for QoI Extraction: %d" % (time.time()-start))


def read_input_file(input_file):
	"""Read Input file txt file and define user-defined variables into global variables"""
	global do_registration
	global do_volumeremeshing
	global do_resample4dmri
	global do_qoiextraction
	global centerline_already_extracted
	global meshing_type
	global data_path
	global mr_subject_id
	global subject_id
	global mr_vtk_dir
	global dia_tmpl
	global sys_tmpl
	global max_iterations
	global remesh_edge_length
	global seg_locations
	global total_timepoints
	global procs

	with open(input_file, 'r') as f:
		for line in f:
			# Skip commented lines
			if not(line.split()) or line.split()[0] == '#':
				continue
			curr_line = line.split('=')
			# Parse for arguments
			if 'do_registration' in curr_line[0]:
				do_registration = eval(curr_line[-1])
			elif 'do_volumeremeshing' in curr_line[0]:
				do_volumeremeshing = eval(curr_line[-1])
			elif 'do_resample4dmri' in curr_line[0]:
				do_resample4dmri = eval(curr_line[-1])
			elif 'do_qoiextraction' in curr_line[0]:
				do_qoiextraction = eval(curr_line[-1])
			elif 'centerline_already_extracted' in curr_line[0]:
				centerline_already_extracted = eval(curr_line[-1])
			elif 'meshing_type' in curr_line[0]:
				meshing_type = curr_line[-1].rstrip().split("'")[1]
			elif 'data_path' in curr_line[0]:
				data_path = curr_line[-1].rstrip().split("'")[1]
			elif 'mr_subject_id' in curr_line[0]:
				mr_subject_id = curr_line[-1].rstrip().split("'")[1]
			elif 'subject_id' in curr_line[0]:
				subject_id = curr_line[-1].rstrip().split("'")[1]
			elif 'mr_vtk_dir' in curr_line[0]:
				mr_vtk_dir = curr_line[-1].rstrip().split("'")[1]
			elif 'dia_tmpl' in curr_line[0]:
				dia_tmpl = curr_line[-1].rstrip().split("'")[1]
			elif 'sys_tmpl' in curr_line[0]:
				sys_tmpl = curr_line[-1].rstrip().split("'")[1]
			elif 'max_iterations' in curr_line[0]:
				max_iterations = int(curr_line[-1].rstrip())
			elif 'total_timepoints' in curr_line[0]:
				total_timepoints = int(curr_line[-1].rstrip())
			elif 'remesh_edge_length' in curr_line[0]:
				remesh_edge_length = float(curr_line[-1].rstrip())
			elif 'seg_locations' in curr_line[0]:
				tmp_seg_locations = curr_line[-1].split("'")[0].split(",")
				seg_locations = [float(x) for x in tmp_seg_locations]
			elif 'procs' in curr_line[0]:
				procs = int(curr_line[-1].rstrip())

"""
USER CALL
"""
# main()
# MAIN FUNCTION
if __name__ == "__main__":

	if len(sys.argv) < 2:
		print("Error: Arguments Missing")
		print("Usage:")
		print(" python {} [Input text file]".format(sys.argv[0]))
		sys.exit(0)

	input_file = sys.argv[1]
	try:
		input = open(input_file, 'r')
		input.close()
	except IOError:
		print('Input Filename cannot be opened/does not exist')
		sys.exit()

	sys.exit(main(input_file))