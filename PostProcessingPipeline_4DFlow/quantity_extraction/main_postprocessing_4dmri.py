"""Script to call postprocessing_4dmri for multiple instances. Specifically used for parallelization
without initializing new nodes every time.

To call:
Input file must contain the paths to the input filenames used by postprocessing_4dmri

=================================================================================================================================
---Created--- | ------Owner------- | Notes---------------------------------------------------------------------------------------
=================================================================================================================================
    07-22-2021  Melody Dong          Created
=================================================================================================================================
"""

import os
import sys
import pdb

import postprocessing_4dmri

def call_postprocessing(input_file):
	input_file_names = []
	with open(input_file, 'r') as f:
		for line in f:
			input_file_names.append(line)

	for file in input_file_names:
		print(file.split('/')[-5])
		try:
			if '\n' in file.split('/')[-1]:
				print(postprocessing_4dmri.main(file[:-1]))
			else:
				print(postprocessing_4dmri.main(file))
		except:
			pdb.set_trace()
			print("ERROR: COULD NOT PROCESS %s" % file.split('/')[-5])
			continue


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

	sys.exit(call_postprocessing(input_file))