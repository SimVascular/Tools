# arterys_mapping
Python class to register Arterys image and corresponding npy velocity data and convert to viewable vtk files
Functionally equivalent to the Matlab scripts and functions driven by ArterysDataProcessing.m

## Background
Scripts in this repository are used to process Arterys image and velocity data to register them and align them correctly onto the same x,y,z coordinates. Registered data can be saved as vtk files (saving to matlab and tecplot in progress).

Data downloaded from Arterys includes image files and a numpy array of velocity data. 

Image files are downloadable from Arterys as .tar or .zip files and can include 4DMRI scans and standard CMRs.

## Usage

arterysdataprocessing.py is the main class with subroutines to extract dicoms from compressed Arterys files, register velocity onto the image, and save final scan into a vtk file. This file uses subroutines from imagesorting.py within the ArterysDataProcessing repo

arterys_mapping.py is the driver script to execute and call arterysdataprocessing.py. This is an example of how to call arterysdataprocessing and can be changed to prescribe the following user options:
* imgFile: [path] to compressed folder of dicom files or parent folder where IMAGE/ has already been extracted
* velFileName: [path] to velocity numpy array from Arterys
* saveFileName: file name prefix to save outputs with
* scalingFactor: Scaling factor (double) to scale all images and velocities to different units. Example- 0.1 scales mm units to cm units
* zeroOrigin: True means coordinates of the image and velocity will be aligned to (0,0,0) on the x,y,z coordinates. False means image and velocity will be aligned to the patient position from metadata found in the dicom file. Aligning to patient position will match the origin of the .vti outputted from the imagesorting.py script.
* alignSimVascular: True means coordinates of image and velocity will be aligned with the same origin as the vti that SimVascular outputs.
* [not functional] medianFiltering: True means applying median filtering to velocity data, false means no filtering
* [not functional] calculateQoI: True means computing vorticity and Q-criterion from the velocity data and spatial coordinates
* [not functional] writeToTecplot: True means data will be written as Tecplot plt files for each phase, false means no Tecplot files will be written
* imagesExtracted: True if compressed folder already extracted
* writeOutput: True if want vtk files with velocity and image written to file
