import argparse
import os
import sys

sys.path.insert(0, '/home/melody/Software/SimVascular-Tests/new-api-tests/applications/create-surface-caps-from-centerlines/')
sys.path.insert(1, '/home/melody/Software/SimVascular-Tests/new-api-tests/graphics/')
from centerlines import Centerlines
from surface import Surface
import graphics as gr
import create_surface_caps as surfcap
import pdb


def parse_args():
	'''Parse command-line arguments.
	'''
	parser = argparse.ArgumentParser()

	parser.add_argument("--path-name", required=True, type=str, default='',
		help="Path name to folder containing surface files.")

	parser.add_argument("--prefix-name", required=True, type=str, default='',
		help="Prefix of all surface files.")

	parser.add_argument("--timepoints", required=True, type=int, default=20,
		help="Total # of timepoints to do this for.")

	parser.add_argument("--clip-distance", type=float, default=0.1, 
		help="The distance from the end of a centerline branch to clip a surface.")

	parser.add_argument("--clip-width-scale", type=float, default=1.0, 
		help="The width multiplied by the centerline branch end radius to define the width of the box used to clip a surface.")

	parser.add_argument("--mesh-scale", type=float, default=1.0, 
		help="The factor used to scale the fe volume meshing edge size. A larger scale creates a coarser mesh. The initial edge size is determined from the largest surface triangle.")

	parser.add_argument("--remesh-scale", type=float, default=1.0, 
		help="The factor used to scale the surface remeshing edge size. A larger scale creates a coarser suface mesh. The initial edge size is determined from the largest surface triangle.")

	args = parser.parse_args()

	if len(sys.argv) == 1:
		parser.print_help()
		sys.exit(1)

	return args

def main():

	# Get command-line arguments.
	args = parse_args()

	cl_path = os.path.join('/'.join(args.path_name.split('/')[:-1]), 'cl_sv')
	if not(os.path.exists(cl_path)):
		os.mkdir(cl_path)

	for t in range(args.timepoints):

		## Create renderer and graphics window.
		win_width = 500
		win_height = 500
		renderer, renderer_window = gr.init_graphics(win_width, win_height)

		surface_file_name = os.path.join(args.path_name, args.prefix_name + '__%02d_surfacemesh.vtk' % t)

		## Read in the segmentation surface.
		surface = Surface(gr, renderer_window, renderer)
		surface.read(surface_file_name)
		gr_geom = gr.add_geometry(renderer, surface.geometry, color=[0.8, 0.8, 1.0])
		surface.vtk_actor = gr_geom 
		#gr_geom.GetProperty().SetOpacity(0.5)

		## Create a Centerlines object used to clip the surface.
		centerlines = Centerlines()
		centerlines.graphics = gr
		centerlines.surface = surface
		centerlines.window = renderer_window 
		centerlines.renderer = renderer
		centerlines.clip_distance = args.clip_distance
		centerlines.clip_width_scale = args.clip_width_scale
		centerlines.remesh_scale = args.remesh_scale
		centerlines.mesh_scale = args.mesh_scale

		# surface.compute_centerlines(data=surface)
		# surface.create_model_automatically(data=centerlines)

		print("---------- Alphanumeric Keys ----------")
		print("a - Compute model automatically for a three vessel surface with flat ends.")
		print("c - Compute centerlines.")
		print("m - Create a model from the surface and centerlines.")
		print("q - Quit")
		print("s - Select a centerline source point.")
		print("t - Select a centerline target point.")
		print("u - Undo the selection of a centerline source or target point.")

		## Create a mouse interactor for selecting centerline points.
		picking_keys = ['s', 't']
		event_table = {
			'a': (surface.create_model_automatically, centerlines),
			'c': (surface.compute_centerlines, surface),
			'm': (centerlines.create_model, surface),
			's': surface.add_centerlines_source_node,
			't': surface.add_centerlines_target_node
		}
		interactor = gr.init_picking(renderer_window, renderer, surface.geometry, picking_keys, event_table)

		## Display window.
		interactor.Start()

		os.rename(os.path.join(args.path_name, args.prefix_name + '__%02d_surfacemesh-centerlines.vtp' % t), os.path.join(cl_path, 'centerlines-result_%02d.vtp' % t))


if __name__ == '__main__':
	main()
