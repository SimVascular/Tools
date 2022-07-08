"""Class for processing Arterys 4D flow MRI files

Given the compressed file of dicom images and velocity numpy array from Arterys

Given an input file with defined pathnames, this script will:
	- Process velocity numpy array
	- Packages image data and velocity data together
	- Currently outputs in VTK (TODO: or Tecplot) - implementing other options in the future
	- TODO: Currently calculates velocity magnitude, vorticity, Q-criterion - implement other calculations in the future

Processes all dicom images in a given directory. To obtain folder of dicoms sorted into a folder of only 
magnitude dicoms, use imagesorting.py. You may choose to run imagesorting.py within this script.

Dependency requirements:
- pydicom
- vmtk, vtk
- numpy

=================================================================================================================================
---Created--- | ------Owner------- | Notes---------------------------------------------------------------------------------------
=================================================================================================================================
   	   1-2019   Nicole Schiavone     Created original Matlab scripts
   	   3-2021	Melody Dong			 Translated to Python
=================================================================================================================================
"""
import os
from os import listdir, path
from os.path import isfile, join
import sys

import copy
import numpy as np
import pydicom
import pdb
import vtk
from vtk.util.numpy_support import vtk_to_numpy, numpy_to_vtk
from vmtk import vtkvmtk, vmtkscripts

sys.path.insert(0, os.path.split(sys.path[0])[0])
import imagesorting



class arterysDataProcessing:

	def __init__(self):
		self.imgFile = ''
		self.magImgDir = ''
		self.velFileName = ''
		self.saveFileName = ''
		self.scalingFactor = 0.1
		self.zeroOrigin = False
		self.alignSimVascular = False
		self.medianFiltering = False
		self.calculateQoI = True
		self.writeToVTK = True
		self.writeToTecplot = False
		self.imagesExtracted = True
		self.writeOutput = True


	def Execute(self):
		# Extract images if not already extracted
		if not(self.imagesExtracted):
			self.extractSortImages(self.imgFile)
			imgDir = os.path.join('/'.join(self.imgFile.split('/')[:-1]), 'IMAGE', 'mag')
			saveDir = os.path.join('/'.join(self.imgFile.split('/')[:-1]))
		else:
			imgDir = os.path.join(self.imgFile, 'IMAGE', 'mag')
			saveDir = os.path.join(self.imgFile)

		# Read Velocity and Filter/Scale
		velocity = self.readVelocity()
		velocity = self.scaleVelocity(velocity)
		if self.medianFiltering:
			velocity = self.medianFiltering(velocity)

		# Read Imageas and map velocity onto it
		# data4dmri = []
		magImgList = [f for f in listdir(imgDir) if os.path.isdir(os.path.join(imgDir, f))]
		for i in range(len(magImgList)):
			print("\tTimepoint: %d" % i)
			# Read Image
			magImgDir = os.path.join(imgDir, 'mag' + str(i+1))
			magImg = self.readMagImg(magImgDir)

			# Map Velocity
			magVel = velocity[i]
			imgVel = self.mapVelocity(magImgDir, magVel, magImg.Image)
			# data4dmri.append(imgVel)

			# Save mapped image and velocity
			if self.writeOutput:
				saveOutputName = os.path.join(saveDir, self.saveFileName + '_%02d.vtk' % (i))
				self.write4dmriData(saveOutputName, imgVel)

		# return data4dmri



	# TODO: Calculate other quantities of interest - vorticity, Q-criterion

	# TODO: Write out 4DMRI data to Tecplot file


	def mapVelocity(self, magImgDir, originalVel, img):
		"""Map velocity array to same shape as image array onto vtk object
		Args:
			- magImgDir (str): path to folder containing dicom files
			- originalVel (npy array): numpy velocity array to map
			- img (vtkImageData): vtk object containing image data
		Returns:
			- img (vtkStructuredGrid): vtk object with velocity and image together
		"""
		# Get header information for image
		firstImgPath = listdir(magImgDir)[0]
		imgHeader = pydicom.dcmread(os.path.join(magImgDir, firstImgPath), force=True)

		# Scan Ras
		scanRas = imgHeader[0x0019, 0x1018].value
		scanRasLoc = imgHeader[0x0019, 0x1019].value

		# reshapes to x, y, z coordinates
		imgNPY = vtk_to_numpy(img.GetPointData().GetArray('ImageScalars'))
		imgNPY = imgNPY.reshape(img.GetDimensions(), order='F')

		# To convert velocity numpy array to same dimensions of img - use swapaxes
		velx = originalVel[0]
		vely = originalVel[1]
		velz = originalVel[2]
		newVelX = np.zeros(imgNPY.shape)
		newVelY = np.zeros(imgNPY.shape)
		newVelZ = np.zeros(imgNPY.shape)		
		velx = velx.swapaxes(0, 2)
		vely = vely.swapaxes(0, 2)
		velz = velz.swapaxes(0, 2)
		# Adjust size of velocity arrays to match size of image array
		minCropIdx_0 = int((imgNPY.shape[0]-velx.shape[0])/2)
		minCropIdx_1 = int((imgNPY.shape[1]-velx.shape[1])/2)
		minCropIdx_2 = int((imgNPY.shape[2]-velx.shape[2])/2)
		newVelX[minCropIdx_0:minCropIdx_0+int(velx.shape[0]), 
				minCropIdx_1:minCropIdx_1+int(velx.shape[1]), 
				minCropIdx_2:minCropIdx_2+int(velx.shape[2])] = velx
		newVelY[minCropIdx_0:minCropIdx_0+int(velx.shape[0]), 
				minCropIdx_1:minCropIdx_1+int(velx.shape[1]), 
				minCropIdx_2:minCropIdx_2+int(velx.shape[2])] = vely
		newVelZ[minCropIdx_0:minCropIdx_0+int(velx.shape[0]), 
				minCropIdx_1:minCropIdx_1+int(velx.shape[1]), 
				minCropIdx_2:minCropIdx_2+int(velx.shape[2])] = velz

		# Try converting to original velocity shape
		newVelX = newVelX.swapaxes(0, 2)
		newVelY = newVelY.swapaxes(0, 2)
		newVelZ = newVelZ.swapaxes(0, 2)

		if scanRas == 'P' or scanRas == 'A':
			cp_newVelY = copy.copy(newVelY)
			cp_newVelZ = copy.copy(newVelZ)
			newVelY = -cp_newVelZ
			newVelZ = cp_newVelY

		numPts = newVelX.shape[0]*newVelX.shape[1]*newVelX.shape[2]
		newVelX = newVelX.ravel().reshape((numPts, 1))
		newVelY = newVelY.ravel().reshape((numPts, 1))
		newVelZ = newVelZ.ravel().reshape((numPts, 1))
		newVelocity = np.append(newVelX, newVelY, 1)
		newVelocity = np.append(newVelocity, newVelZ, 1)

		# Convert velocity to vtk
		velVtk = numpy_to_vtk(newVelocity)
		velVtk.SetName('velocity')

		# Add velocity to image
		img.GetPointData().AddArray(velVtk)

		# Convert ImageData to StructuredGrid
		imgConvert = vtk.vtkImageDataToPointSet()
		imgConvert.SetInputData(img)
		imgConvert.Update()

		return imgConvert.GetOutput()


	def scaleVelocity(self, originalVel):
		"""Scale the velocity based on the scalingFactor to match image scaling
		Args:
			- velocity numpy array
		"""
		return originalVel*self.scalingFactor


	def medianFilteringVelocity(self, originalVel):
		"""Apply median filtering to velocity numpy array, use a default filter of 
		3 pixels neighboring.

		Args:
			- velocity numpy array
		"""
		# Check that user wants median filtering
		if self.medianFiltering:
			from scipy import ndimage
			return ndimage.median_filter(originalVel, 3)


	def write4dmriData(self, saveOutputName, dataVTK):
		"""After filtering and centering magnitude image and velocity, write out to a file per timepoint
		Args:
			- saveOutputName (str): full path + file name to save output to
			- dataVTK (vtk): vtk object with image + velocity + other QoIs
		"""
		vtkDataWriter = vtk.vtkDataSetWriter()
		vtkDataWriter.SetFileName(saveOutputName)
		vtkDataWriter.SetInputData(dataVTK)
		vtkDataWriter.Update()


	def readMagImg(self, magImgDir):
		"""Read images in one timepoint into vtk object. Magnitude images must be organized before
		running using imagesorting.py. Center images based on user-based preferences.

		Args:
			magImgDir: pathname of folder containing magnitude dicom images from one timepoint

		Returns:
			vtk object of image scalars for volumetric object of one timepoint
		"""
		img = listdir(magImgDir)[0]

		# Read Image
		reader = vmtkscripts.vmtkImageReader()
		reader.InputFileName = os.path.join(magImgDir, img)
		reader.Execute()

		# Apply Scaling to Image
		spacing = reader.Image.GetSpacing()
		reader.Image.SetSpacing(spacing[0]*self.scalingFactor, spacing[1]*self.scalingFactor, spacing[2]*self.scalingFactor)

		# Set Origin if user specified different origin
		if self.zeroOrigin:
			reader.Image.SetOrigin(0, 0, 0)
		if self.alignSimVascular:
			origin = reader.Image.GetOrigin()
			reader.Image.SetOrigin(-origin[0]*self.scalingFactor, -origin[1]*self.scalingFactor, origin[2]*self.scalingFactor)

		return reader


	def readVelocity(self):
		"""Read velocity .npy file and return the velocity numpy array"""
		vel_file = os.path.join(self.velFileName)
		velocity = np.load(vel_file)
		return velocity


	def readImageInfo(self):
		"""Read dicom meta-data info for images from Arterys"""
		img = listdir(self.imgDir)
		img_header =  pydicom.dcmread(os.path.join(self.imgDir, img[0]), force=True)
		return img_header


	def extractSortImages(self, imgFile):
		"""Extract compressed image files and organize by calling imagesorting.py subroutines.
		Will uncompress files to the parent directory where the file is living.
		"""
		# Check that file is a compressed zip/tar file
		if imgFile[-4:] == '.zip' or imgFile[-4:] == '.tar' or imgFile[-4:] == '.tgz':
			imagesorting.process_sort_4DMR_data(True, self.scalingFactor, imgFile, False)
		else:
			sys.exit("ERROR: Attempted to extract compressed file, but none were given. \n \
					  imgFile = %s" % imgFile)
