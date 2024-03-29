import argparse
import sys
import os
import tabulate
import numpy as np
import lmfit
import scipy
import astropy.io
import palettable
import matplotlib.pyplot as plt
from cycler import cycler
import tqdm
import trackMasking as tk
import glob
from pathlib import Path
import datetime

import readFits
import PoissonGausFit
import DamicImage

# from plotAllAmplifierPixDistShifter import plotAllAmplifierPixelDist
# from plotAllAmplifierColumnValues import plotAllAmplifierColumnValues
# from plotImageMetrics import plotImageMetrics
# from createJoinedFitsFile import createJoinedFitsFile

from shifterPlots import *

# Specific name of files to be processed. Need to update when changing science runs
processDir = "/data/damic/snolab/upgrade/processed/science/run6b"
filenameStruct = {"L1":[], "U1":[], "L2":[],"U2":[],}

imageMask = {"U1":np.load("/home/apiers/snolab/mask/run5a_1U_course_mask.npy"), "L1":np.load("/home/apiers/snolab/mask/run5a_1L_course_mask.npy"), "U2":np.load("/home/apiers/snolab/mask/run5a_2U_course_mask.npy"), "L2":np.load("/home/apiers/snolab/mask/run5a_2L_course_mask.npy") }


if __name__ == '__main__':

	# Code to prepare plots for the DAMIC at SNOLAB shifters

	# Get all the average fits file in the science data folder
	avgImageFiles = glob.glob(os.path.join(processDir, "avg/avg_*.fits"), recursive=True)

	searchDate = datetime.date.today() - datetime.timedelta(0)
	#searchDate = datetime.date(2023, 1, 13)
	# Filter by ones that were processed today, use these images for the shifter
	avgImageFilesToday = [f for f in avgImageFiles if datetime.date.fromtimestamp(os.path.getctime(f)) == searchDate]

	# Parse filenames into a more regular structure to be used by plotting functions
	vImageID = []
	for file in avgImageFilesToday:

		directory, filename = os.path.split(file)

		# Get info from the filename
		splitFilename = os.path.splitext(filename)[0].split("_")

		ccdNumber, imageID, amplifier = splitFilename[-3:]
		vImageID.append(int(imageID))

		# print("%s, %s, %s, %s"%(runID, ccdNumber, imageID, amplifier))

		filenameStruct[amplifier+ccdNumber].append( (imageID, file))


	# Get output directory to save files, make new dir if necessary
	inputdir = os.path.dirname(avgImageFilesToday[0])
	outputdir = (Path(inputdir).parents[0]).joinpath("monitor")
	# Create the directory if it doesn't exist
	if not outputdir.exists():
		os.mkdir(outputdir, mode=0o770)



	# Plot pixel distribution
	print("Creating pixel distribution plot...")
	figPix, axsPix = plotAllAmplifierPixelDist(filenameStruct, imageMask)
	figPix.savefig(outputdir.joinpath("{:s}_pixel_dist_imageid_{}_{}.pdf".format(searchDate.strftime("%Y-%m-%d"), np.min(vImageID), np.max(vImageID))), bbox_inches="tight")

	# Plot CCD median column values
	print("Creating median column value plot...")
	figCol, axsCol = plotAllAmplifierColumnValues(filenameStruct, imageMask)
	figCol.savefig(outputdir.joinpath("{:s}_med_col_val_imageid_{}_{}.pdf".format(searchDate.strftime("%Y-%m-%d"), np.min(vImageID), np.max(vImageID))), bbox_inches="tight")

	# Plot image metrics
	print("Creating image metrics plot...")
	figMet, axsMet = plotImageMetrics(filenameStruct, imageMask)
	figMet.savefig(outputdir.joinpath("{:s}_image_metrics_imageid_{}_{}.pdf".format(searchDate.strftime("%Y-%m-%d"), np.min(vImageID), np.max(vImageID))), bbox_inches="tight")


	# Saving combined fits images
	print("Combining fits files...")
	combinedFitsFilenames = outputdir.joinpath("{:s}_combined_fits_imageid_{}_{}".format(searchDate.strftime("%Y-%m-%d"), np.min(vImageID), np.max(vImageID)))
	createJoinedFitsFile(filenameStruct, combinedFitsFilenames)
