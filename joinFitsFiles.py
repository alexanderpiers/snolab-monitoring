import numpy as np
import os
import glob
from astropy.io import fits
import argparse

nrows = 210

def createJoinedFitsFiles(filenames, outdir, chunksize=30):

	# sort files first, using image ID to sort
	filenames.sort()
	
	imageID = [ int(x.split("_")[-2]) for x in filenames]
	noutfiles = len(imageID) // chunksize + 1
	startID = imageID[0]

	for i in range(noutfiles):
	

		print("Joining image ids {}-{}".format(startID + i*chunksize, np.min([(i+1)*chunksize, len(imageID)])))
		# get some info on the size and fits header for writing
		hdu = fits.open(filenames[i*chunksize])
		data = hdu[0].data
		variance = hdu[1].data
		header = hdu[0].header

		nrows = data.shape[0]
		# do the joining
		dataJoined = np.zeros((data.shape[0] * chunksize, data.shape[1]))
		varianceJoined = np.zeros((data.shape[0] * chunksize, data.shape[1]))

		for j in range(chunksize):
			try:
				hdu = fits.open(filenames[i*chunksize + j])
				data = hdu[0].data
				variance = hdu[1].data

				dataJoined[j*nrows:(j+1)*nrows, :] = data
				varianceJoined[j*nrows:(j+1)*nrows, :] = variance
			except IndexError:
				pass
				

		# now write the file
		infname = os.path.split(filenames[0])[1]
		outfilename = "joined_" + "_".join(infname.split("_")[:-2]) + "_{:03}_".format(i+1) + infname.split("_")[-1]
		print(outfilename)

		hduList = fits.HDUList()
		hduList.append(fits.ImageHDU(data=dataJoined, header=header))
		hduList.append(fits.ImageHDU(data=varianceJoined, header=header))
		# Save file as a new fits
		hduList.writeto(os.path.join(outdir, outfilename), overwrite=True)



if __name__ == "__main__":

	print("combining fits files")

	parser = argparse.ArgumentParser(description="Join multiple fits files")
	
	parser.add_argument("-f", "--filenames", nargs="*", help="files to process")
	parser.add_argument("-o", "--outdir", help="output directory")
	parser.add_argument("-c", "--chunk", default=40, type=int, help="number of fits to combine together")

	commandArguments = parser.parse_args()

	filenames = commandArguments.filenames
	outdir = commandArguments.outdir
	chunksize = commandArguments.chunk


	createJoinedFitsFiles(filenames, outdir, chunksize=chunksize)
