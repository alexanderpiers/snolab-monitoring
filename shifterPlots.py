import sys
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import numpy as np
import palettable
import argparse
import scipy.stats
from matplotlib.gridspec import GridSpec
import astropy.io
import glob
import os
import datetime
from cycler import cycler


import PoissonGausFit as pgf
import readFits
import DamicImage
import lmfit
import trackMasking as tk

nElectronMask = 20
imageMask = {"U1":np.load("/home/apiers/snolab-monitor/mask/run0d/U1_mask_RUNID133-139.npy"), "L1":np.load("/home/apiers/snolab-monitor/mask/run0d/L1_mask_RUNID133-139.npy"), "U2":np.load("/home/apiers/snolab-monitor/mask/run0d/U2_mask_RUNID133-139.npy"), "L2":np.load("/home/apiers/snolab-monitor/mask/run0d/L2_mask_RUNID133-139.npy") }
processDir = "/data/damic/snolab/upgrade/processed/science"

def plotAllAmplifierPixelDist(filestructure, columnMask, gain=-1, bw=0.05):

	fig, axs = plt.subplots(2, 2, figsize=(16, 10))
	imageID = 0
	axflat = axs.flatten()

	# For each amplifier plot pixel distribution
	for i, (key, images) in enumerate(filestructure.items()):

		imageID = np.array(images)[:, 1].astype(int)

		# Read all the images in the file struct object
		allImages = np.array([])
		for fileinfo in images:
			filename = fileinfo[-1]

			hdu = astropy.io.fits.open(filename)
			data = hdu[0].data

			# Stack all data. Stacking order is not important because we are just doing the pixel distribution in the end
			if allImages.size == 0:
				allImages = data
			else:
				allImages = np.vstack((allImages, data))

		# Perform fit to raw data
		img = DamicImage.DamicImage(allImages, bw=bw, reverse=False, minRange=np.abs(20*gain))
		minres = pgf.computeGausPoissDist(img, aduConversion=gain, npoisson=10, darkCurrent=-0.02, offset=-img.med, sigma=-0.2)
		paramsRaw = minres.params

		# Perform fit to masked data
		clusterMask = tk.mask(allImages, nElectronMask, 60, 10)
		extendedColMask = np.reshape( np.tile( columnMask[key], allImages.shape[0]), allImages.shape, order="C")
		fullMask = np.logical_or(extendedColMask, np.logical_not(clusterMask))
		imgMask = DamicImage.DamicImage(allImages[np.logical_not(fullMask)], bw=0.05, reverse=False, minRange=np.abs(20*gain))
		minresMask = pgf.computeGausPoissDist(imgMask, aduConversion=gain, npoisson=10, darkCurrent=-0.02, offset=-imgMask.med, sigma=-0.2)
		paramsMask = minresMask.params

		# Plot results
	    # Plot pixel distribution
		axflat[i].errorbar(img.centers, img.hpix, yerr=np.sqrt(img.hpix), fmt="ok", markersize=3, alpha=0.7)
		axflat[i].errorbar(imgMask.centers, imgMask.hpix, yerr=np.sqrt(imgMask.hpix), fmt="ob", markersize=3, alpha=0.7)

		# Plot fit results
		x = np.linspace(img.centers[0], img.centers[-1], 20000)
		axflat[i].plot(x, pgf.fGausPoisson(x, *pgf.paramsToList(paramsRaw)), "--k", linewidth=3, label=r"Raw. $\sigma$=%.2f e-, $\lambda$=%.2g e- / pix"%(paramsRaw["sigma"].value / paramsRaw["ADU"].value, paramsRaw["lamb"].value))
		axflat[i].plot(x, pgf.fGausPoisson(x, *pgf.paramsToList(paramsMask)), "--b", linewidth=3, label=r"Mask. $\sigma$=%.2f e-, $\lambda$=%.2g e- / pix"%(paramsMask["sigma"].value / paramsMask["ADU"].value, paramsMask["lamb"].value))

		# Set up axis and legend
		axflat[i].set_ylabel(r"Counts (%0.2f e-)${}^{-1}$"%(bw), fontsize=14)
		axflat[i].set_xlabel("Pixel Value (e-)", fontsize=14)
		axflat[i].legend(fontsize=14)
		axflat[i].set_yscale("log")
		axflat[i].set_ylim(0.05, 1.5*np.max(img.hpix))
		axflat[i].set_xlim(-1, 5)
		axflat[i].set_title(key, fontsize=14)

	fig.suptitle("Pixel Distribution of ImageID: {}-{}".format(np.min(imageID), np.max(imageID)), fontsize=16)


	return fig, axs

def plotAllAmplifierColumnValues(filestructure, columnMask, gain=-1, bw=0.05):

	fig, axs = plt.subplots(2, 2, figsize=(16, 10))
	imageID = 0
	axflat = axs.flatten()

	# For each amplifier plot pixel distribution
	for i, (key, images) in enumerate(filestructure.items()):

		imageID = np.array(images)[:, 1].astype(int)

		# Read all the images in the file struct object
		allImages = np.array([])
		for fileinfo in images:
			filename = fileinfo[-1]

			hdu = astropy.io.fits.open(filename)
			data = hdu[0].data

			# Stack all data. Stacking order is not important because we are just doing the pixel distribution in the end
			if allImages.size == 0:
				allImages = data
			else:
				allImages = np.vstack((allImages, data))

		# Compute median value of each column

		medColValue = np.median(allImages, axis=0)
		axflat[i].plot(medColValue, "ok")
		axflat[i].fill_between(np.arange(allImages.shape[1]), 0, 1, where=columnMask[key], alpha=0.4, color="r", transform=axflat[i].get_xaxis_transform())

		# Plot fit results
		
		# Set up axis and legend
		axflat[i].set_ylabel(r"Median Column Value (e-)", fontsize=14)
		axflat[i].set_xlabel("Column Number", fontsize=14)
		# axflat[i].legend(fontsize=14)
		# axflat[i].set_yscale("log")
		# axflat[i].set_ylim(0.05, 1.5*np.max(img.hpix))
		axflat[i].set_xlim(0, allImages.shape[1])
		axflat[i].set_title(key, fontsize=14)

	fig.suptitle("Median Column Values of ImageID: {}-{}".format(np.min(imageID), np.max(imageID)), fontsize=16)


	return fig, axs

def plotImageMetrics(filestructure, columnMask, bw=0.05):

	fig, axs = plt.subplots(3, 1, figsize=(16, 10), sharex=True)
	imageID = 0
	
	color = palettable.colorbrewer.qualitative.Paired_6.mpl_colors

	for ax in axs:
		ax.set_prop_cycle(cycler(color=color))
		ax.tick_params(labelsize=14)
		ax.grid(True, which="both")	

	# For each amplifier plot pixel distribution
	for i, (key, images) in enumerate(filestructure.items()):

		imageID = np.array(images)[:, 1].astype(int)

		# Read all the images in the file struct object
		allImages = np.array([])
		darkCurrent = []
		darkCurrentErr = []
		noise = []
		noiseErr = []
		calibration = []

		for fileinfo in images:
			filename = fileinfo[-1]

			# Read pixel distribution
			hdu = astropy.io.fits.open(filename)
			header = hdu[0].header
			data = hdu[0].data

			# Apply a mask for clusters. Things >= 5 e-
			clusterMask = np.logical_not(tk.mask(data, 20, 50, 10))
			extendedColMask = np.reshape( np.tile( columnMask[key], data.shape[0]), data.shape, order="C")
			fullMask = np.logical_or(extendedColMask, clusterMask)


			img = DamicImage.DamicImage(data[np.logical_not(fullMask)], reverse=False, bw=bw, minRange=10)

			# Perform fit of pixel distribution
			minresult = pgf.computeGausPoissDist(img, offset=-0.01, aduConversion=-1, darkCurrent=-0.05, npoisson=10, sigma=-0.18)

			fitparams = pgf.parseFitMinimum(minresult)
			calibration.append(header["CAL"])
			noise.append(fitparams["sigma"][0])
			darkCurrent.append(fitparams["lambda"][0])
			noiseErr.append(fitparams["sigma"][1])
			darkCurrentErr.append(fitparams["lambda"][1])

		# Preprocess error bars in the case the fit didn't converge
		amplifierNoiseErr = np.array(noiseErr)
		amplifierDCErr = np.array(darkCurrentErr)
		amplifierNoiseErr[ amplifierNoiseErr == None ] = 0
		amplifierDCErr[ amplifierDCErr == None ] = 0

		# Plot results
		axs[0].plot(imageID, calibration, "o", label=key)
		axs[0].set_ylabel("k (ADU / e-)", fontsize=14)

		#  Plotting the noise
		axs[1].errorbar(imageID, np.array(noise), yerr=amplifierNoiseErr, fmt="o", label=key)
		axs[1].set_ylabel("Noise (e-)", fontsize=14)

		#  Plotting the DC value
		axs[2].errorbar(imageID, darkCurrent, yerr=amplifierDCErr, fmt="o", label=key)
		axs[2].set_ylabel(r"$\lambda$ (e- / pix)", fontsize=14)
		axs[2].set_xlabel("RunID", fontsize=16)
		


	for ax in axs:
		ax.legend(fontsize=16)
	fig.suptitle("ImageID: {}-{}".format(np.min(imageID), np.max(imageID)), fontsize=16)


	return fig, axs

def createJoinedFitsFile(filestructure, basefilename):

	for i, (key, images) in enumerate(filestructure.items()):

		imageIDs = np.array(images)[:, 1].astype(int)

		imageIDRange = np.arange(np.min(imageIDs), np.max(imageIDs)+1)

		# Read all the images in the file struct object
		allImages = np.array([])
		header = ""

		for imageID in imageIDRange:
			idx = np.where(imageIDs == imageID)[0][0]

			filename = images[idx][2]

			hdu = astropy.io.fits.open(filename)
			data = hdu[0].data
			header = hdu[0].header

			# Stack all data. Stacking order is not important because we are just doing the pixel distribution in the end
			if allImages.size == 0:
				allImages = data
			else:
				allImages = np.vstack((allImages, data))

		# save the joint fits file
		header["COM"] = "Several fits image joined together."
		outfilename = basefilename.parent.joinpath(basefilename.name + "_" + key + ".fits")
		astropy.io.fits.writeto(outfilename, allImages, header=header, overwrite=True)




if __name__ == '__main__':
	
	filenameStruct = {"U1":[], "L1":[], "U2":[], "L2":[]}

	avgImageFiles = glob.glob(os.path.join(processDir, "**/avg_*.fits"), recursive=True)

	# Filter by ones that were processed today, use these images for the shifter
	avgImageFilesToday = [f for f in avgImageFiles if datetime.date.fromtimestamp(os.path.getctime(f)) == (datetime.date.today() - datetime.timedelta(days=1)) ]


	# Parse filenames into a more regular structure to be used by plotting functions
	for file in avgImageFilesToday:

		directory, filename = os.path.split(file)

		# Get info from the filename
		splitFilename = os.path.splitext(filename)[0].split("_")

		runID, ccdNumber, imageID, amplifier = splitFilename[-4:]

		print("%s, %s, %s, %s"%(runID, ccdNumber, imageID, amplifier))

		filenameStruct[amplifier+ccdNumber].append( (runID, imageID, file))

	fig, axs = plotAllAmplifierPixelDist(filenameStruct, imageMask)

	plt.show()