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

import readFits
import PoissonGausFit
import DamicImage

skoffset = 1
ltaAmplifier = {"L":2, "U":4}
bw = 0.05

imageMask = {"U1":np.flip(np.load("/home/apiers/snolab-monitor/mask/run1a/U1_mask_RUNID133-139.npy")), "L1":np.load("/home/apiers/snolab-monitor/mask/run1a/L1_mask_RUNID_1a_53-61.npy"), "U2":np.flip(np.load("/home/apiers/snolab-monitor/mask/run1a/U2_mask_RUNID133-139.npy")), "L2":np.load("/home/apiers/snolab-monitor/mask/run1a/L2_mask_RUNID133-139.npy") }

if __name__ == '__main__':
	

	# Parse relevent arguments
	parser = argparse.ArgumentParser(description="Command line interface for damic-m image preprocessing."	)

	# Add arguments
	parser.add_argument("-f", "--filenames", nargs="*",	help="Processes all files passed to the code",)
	parser.add_argument(
		"-l", 
		"--lta",
		action="store_true",
		help="For LTA data")
	parser.add_argument("-q", "--quick", action="store_true", help="Quick analysis from header information")


	commandArgs = parser.parse_args()

	# Set command line arguments
	filenames = commandArgs.filenames
	uselta = commandArgs.lta
	quickAnalysis = commandArgs.quick

	print(np.sum(imageMask["U1"]))
	print(np.sum(imageMask["L1"]))
	print(np.sum(imageMask["U2"]))
	print(np.sum(imageMask["L2"]))


	runID = {"U1":[], "L1":[], "U2":[], "L2":[] }
	imageEndTime = {"U1":[], "L1":[], "U2":[], "L2":[] }
	calibrationConstant = {"U1":[], "L1":[], "U2":[], "L2":[] }
	noise = {"U1":[], "L1":[], "U2":[], "L2":[] }
	darkCurrent = {"U1":[], "L1":[], "U2":[], "L2":[] }
	noiseErr = {"U1":[], "L1":[], "U2":[], "L2":[] }
	darkCurrentErr = {"U1":[], "L1":[], "U2":[], "L2":[] }

	for file in tqdm.tqdm(filenames):

		# Parse runid, CCD number
		filebase, ext = os.path.splitext(file)

		runid = int(filebase.split("_")[-2])
		ccdNumber = int(filebase.split("_")[-3])
		ltaRunID = int(filebase.split("_")[-2])
		amplifier = filebase.split("_")[-1]

		key = amplifier+str(ccdNumber)

		if quickAnalysis:
			# Read header information for the necessary info

			hdu = astropy.io.fits.open(file)

			header = hdu[0].header

			runID[key].append(runid);
			calibrationConstant[key].append(header["CAL"])
			noise[key].append(header["NOISE"])
			darkCurrent[key].append(header["LAMBDA"])



		else:
			# Perform fit of pixel distribution

			# Read pixel distribution
			hdu = astropy.io.fits.open(file)
			header = hdu[0].header
			data = hdu[0].data

			# Apply a mask for clusters. Things >= 5 e-
			clusterMask = np.logical_not(tk.mask(data, 20, 50, 10))
			extendedColMask = np.reshape( np.tile( imageMask[key], data.shape[0]), data.shape, order="C")
			fullMask = np.logical_or(extendedColMask, clusterMask)


			img = DamicImage.DamicImage(data[np.logical_not(fullMask)], reverse=False, bw=bw, minRange=10)

			# Perform fit of pixel distribution
			minresult = PoissonGausFit.computeGausPoissDist(img, offset=-0.01, aduConversion=-1, darkCurrent=-0.05, npoisson=10, sigma=-0.18)

			fitparams = PoissonGausFit.parseFitMinimum(minresult)
			runID[key].append(runid);
			calibrationConstant[key].append(header["CAL"])
			noise[key].append(fitparams["sigma"][0])
			darkCurrent[key].append(fitparams["lambda"][0])
			noiseErr[key].append(fitparams["sigma"][1])
			darkCurrentErr[key].append(fitparams["lambda"][1])


			# if key == "U1":
			# 	print(lmfit.fit_report(minresult))
			# 	print(np.sum(fullMask)/data.size)


	fig, axs = plt.subplots(3, 1, figsize=(16, 8), sharex=True)
	color = palettable.colorbrewer.qualitative.Paired_6.mpl_colors

	for ax in axs:
		ax.set_prop_cycle(cycler(color=color))
		ax.tick_params(labelsize=14)
		ax.grid(True, which="both")

	for key in runID:

		if quickAnalysis:
			# Plotting the calibration constant
			axs[0].plot(runID[key], calibrationConstant[key], "o", label=key)
			axs[0].set_ylabel("k (ADU / e-)", fontsize=14)

			#  Plotting the noise
			axs[1].plot(runID[key], np.array(noise[key]) / np.array(calibrationConstant[key]), "o", label=key)
			axs[1].set_ylabel("Noise (e-)", fontsize=14)


			#  Plotting the DC value
			axs[2].plot(runID[key], darkCurrent[key], "o", label=key)
			axs[2].set_ylabel(r"$\lambda$ (e- / pix)", fontsize=14)
			axs[2].set_xlabel("RunID", fontsize=16)
		else:

			# Preprocess error bars in the case the fit didn't converge
			amplifierNoiseErr = np.array(noiseErr[key])
			amplifierDCErr = np.array(darkCurrentErr[key])
			amplifierNoiseErr[ amplifierNoiseErr == None ] = 0
			amplifierDCErr[ amplifierDCErr == None ] = 0


			# Plotting the calibration constant
			axs[0].plot(runID[key], calibrationConstant[key], "o", label=key)
			axs[0].set_ylabel("k (ADU / e-)", fontsize=14)

			#  Plotting the noise
			axs[1].errorbar(runID[key], np.array(noise[key]), yerr=amplifierNoiseErr, fmt="o", label=key)
			axs[1].set_ylabel("Noise (e-)", fontsize=14)

			#  Plotting the DC value
			axs[2].errorbar(runID[key], darkCurrent[key], yerr=amplifierDCErr, fmt="o", label=key)
			axs[2].set_ylabel(r"$\lambda$ (e- / pix)", fontsize=14)
			axs[2].set_xlabel("RunID", fontsize=16)


	for ax in axs:
		ax.legend(fontsize=16)

	plt.show()
