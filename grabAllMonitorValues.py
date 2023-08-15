import numpy as np
import matplotlib.pyplot as plt
import argparse
from astropy.io import fits
import glob
import tqdm
import datetime
import os
import collections


ampFileStructure = ["*_1_*_U.fits", "*_1_*_L.fits", "*_2_*_U.fits", "*_2_*_L.fits"]


if __name__ == '__main__':
	

	# Setup command line arguments
	parser = argparse.ArgumentParser(description="Compute the exposure of DAMIC at SNOLAB experiment")

	parser.add_argument("-d", "--indirs", nargs="*", help="Directories to search for fits files")
	parser.add_argument("-o", "--outfile", default="test.npy", help="Output file")

	args = parser.parse_args()

	indirs = args.indirs
	outfile = args.outfile

	# Define array of things we are tracking
	# Structur is [1U, 1L, 2U, 2L]
	maxNumberOfImages = int(1e4)
	run = np.zeros((maxNumberOfImages, 4))
	date = np.zeros((maxNumberOfImages, 4))
	imageID = np.zeros((maxNumberOfImages, 4))
	calibration = np.zeros((maxNumberOfImages, 4))
	sigma = np.zeros((maxNumberOfImages, 4))


	# get all the fits files (only need one amplifier, will multiply by four)
	totalCount = 0
	for k, directory in enumerate(indirs):

		print(directory)
		for i, amp in enumerate(ampFileStructure):
			print(amp)
			files = sorted(glob.glob(os.path.join(directory, amp ) ) )
			nImagesInDataSet = len(files)
			print(nImagesInDataSet)
			# Now get the header information from all the files
			for j, fname in enumerate(files):
				with fits.open(fname) as hdu:

					# image id
					date[totalCount+j, i] = datetime.datetime.strptime(hdu[0].header["DATESTART"], "%Y-%m-%dT%H:%M:%S").timestamp()
					imageID[totalCount+j, i] = int( fname.split("_")[-2] )
					calibration[totalCount+j, i] = float( hdu[0].header["CAL"] )
					sigma[totalCount+j, i] = float( hdu[0].header["NOISE"] )
					run[totalCount+j, i] = k


			if(i == len(ampFileStructure) - 1):
				totalCount+= nImagesInDataSet

	np.savez(outfile, date=date[:int(totalCount), :], imageid=imageID[:int(totalCount), :], k=calibration[:int(totalCount), :], sigma=sigma[:int(totalCount), :], run=run[:int(totalCount), :], allow_pickle=True)


	# Open the files, read the header and compute the exposure
