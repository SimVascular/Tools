"""Processing of Arterys 4D flow MRI files to combine image and velocity to a single readable file

Main script that interacts with ArterysDataProcessing class to carry out functional subroutines to:
	- extract dicom files from zipped folder of images in Arterys scan
	- read images
	- read and filter velocity
	- write file with image and velocity mapped together

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

import arterysdataprocessing as adp
import pdb


# Call arterysDataProcessing class object
def arterysMapping(img, vel, saveName, imageExtracted):
	
	scan = adp.arterysDataProcessing()
	scan.imgFile = img
	scan.velFileName = vel
	scan.saveFileName = saveName
	scan.imagesExtracted = imageExtracted
	scan.writeOutput = True
	scan.Execute()



"""
USER CALL
"""
if __name__ == "__main__":

	if len(sys.argv) != 5 or '-h' in sys.argv:
		print("Usage:")
		print(" python {} [imageFile] [velFile] [saveFileName] [imageExtracted]".format(sys.argv[0]) + " \n \
					imageFile \t:\t [path] to compressed folder of dicom files or parent folder where IMAGE/ has already been extracted\n \
					velFile \t:\t [path] to velocity numpy array from Arterys \n \
					saveFileName \t:\t file name prefix to save outputs with \n \
					imageExtracted \t:\t True or False if compressed folder already extracted")	
		sys.exit(0)

	img = sys.argv[1]
	vel = sys.argv[2]
	saveName = sys.argv[3]
	imageExtracted = eval(sys.argv[4])
	sys.exit(arterysMapping(img, vel, saveName, imageExtracted))