#!/usr/bin/env/python
##============================================================

##============================================================

import numpy
import scipy
import os
import sys; sys.path.append("../misc/")
import pyfits
import Image
import time
import getopt
from scipy.signal import fftconvolve as convolve
from ImageMethods import linear_rescale, stretch_params,\
												 png_pix, pngcropwhite, writearray,\
												 fits_pix, fits_hdr, makeFITS
from OSMethods import file_seek
from common import lineno
from sys import argv

##============================================================

## Global variables
from common import cutoff

##============================================================
##============================================================

def main(argv):
	"""
	NAME
	  DeconvolveToTargetPSF.py

	PURPOSE
		1. Uses a PSF image to construct a target 2D Gaussian PSF
		2. Finds convolution kernel which maps between the observed and target PSFs
		3. Deconvolves an entire image with this kernel

	COMMENTS
		- Uses standard deviation to characterise Gaussian rather than FWHM;
		f=sqrt(2ln2)*s.
	  
	USAGE
		python DeconvolveToTargetPSF.py [flags] function args
	
	FUNCTION
		1. psf
		2. kern
		3. imdec
	
	ARGS
		1. None
		2. Kernel-size
			 [9,9]
		3. Image-file  kernel-file
		(See examples)

	FLAGS
	  -h            Print this message

	INPUTS
		Depends on function selected		

	OPTIONAL INPUTS

	OUTPUTS
	  stdout        Useful information

	EXAMPLES
		NEED EXAMPLES
	                    
	BUGS
		- There is some disagreement between my convolution method and SciPy's.
		- Can't save directly to PNG in GPSF() -- save as FITS then convert.
		
		- Use flags rather than relying on argv
	
	TO DO
		- Make AsecToPix method
		
	HISTORY
	  2012-06-12 started Sandford (NYU)
	"""
	
	## -------------------------------------------------------------------
	## Options and help-string
	
	try:
		opts, args = getopt.getopt(argv, "hvpk:d:",["help","psf","kern","imdec"])
	except getopt.GetoptError, err:
		## Print help information and exit:
		print str(err)
		print main.__doc__
		return
			   
	for o,a in opts:
		print o,a
		if o in ("-h", "--help"):
			print main.__doc__
			return
		elif o is "-v":
			vb=True
		else:
			assert False, "Unhandled option"
	  
		"""
		elif o in ("-p", "--psf"):
			GPSF()	### Should have arguments!
		elif o in ("-k", "--kern"):
			all_kernels()	###
		elif o in ("-d", "--imdec"):
		###
		"""	  
	## -------------------------------------------------------------------
	
	t_init = time.time()
	
	if args[1]=="psf":
		all_psfs(args)
	elif args[1]=="kern":
		## Find kernel dimensions from CLA -- e.g. "[5,5]"
		kdim=[int(x) for x in args[2][1:-1].split(",")]
		all_kernels(kdim)
	elif args[1]=="imdec":
		deconvolve_image(args[2],args[3])	
	
	return

##============================================================


##============================================================
##============================================================
## PSF Model
##============================================================
##============================================================

## Copy all the psf models in directory to Gaussians
def all_psfs(topdir="~/Dropbox/Hogg_2012/Data/PSFs/", filestring="snap_*.fits"):	
	for fitsfile in file_seek(topdir, filestring):
		GPSF(fitsfile)
	return

##============================================================

## Create a Gaussian PSF like the original image

def GPSF(filename):
	## Get 2D Gaussian data
	Gdat = Gauss_2D(*moments(fits_pix(filename)))
	##		Horrible method: write fits then convert to png
	if os.path.exists("temp.fits"): os.remove("temp.fits")
	hdu = pyfits.PrimaryHDU(Gdat)
	hdu.writeto("temp.fits")
	if zoom==True:
		outfile = os.path.splitext(filename)[0]+"_Gauss_z8.png"
		os.system("ds9 temp.fits -colorbar no -minmax -zoom 8 -saveimage png "+outfile+" -exit")
	else:
		outfile = os.path.splitext(filename)[0]+"_Gauss.png"
		os.system("ds9 temp.fits -colorbar no -minmax -saveimage png "+outfile+" -exit")
	os.remove("temp.fits")
	return Gdat
	
##============================================================
## From data, extract parameters needed to make a naive Gaussian fit
	## Argument data is 2D array

def moments(data):
	dim = data.shape
	L = max(dim)
	if dim[0]!=dim[1]: print "DeconvolveToTargetPSF.py: moments: image not square. Continuing."
	total = data.sum()
	X,Y = numpy.indices(dim)
	m_x = (X*data).sum()/total
	m_y = (Y*data).sum()/total
	ccol = data[:, int(m_y)]
	crow = data[int(m_x), :]
	width_x = numpy.sqrt(abs((numpy.arange(dim[1])-m_y)**2*ccol).sum()/ccol.sum())
	width_y = numpy.sqrt(abs((numpy.arange(dim[0])-m_x)**2*crow).sum()/crow.sum())
	return [L, width_x, width_y, m_x, m_y, total]
  
##============================================================
## Generate 2D Gaussian array

def Gauss_2D(size, w_x, w_y, x0, y0, scale=1.0):
	## Create mesh
	x=numpy.arange(0.0,size,1.0,float)
	y=x[:,numpy.newaxis]
	## Compute Gaussian
	gaussian=numpy.exp(-0.5*(((x-x0)/w_x)**2+((y-y0)/w_y)**2))
	## Cut out pixels below threshold
	gaussian[gaussian<7.0e-6]=0.0
	## Normalise
	gaussian/=gaussian.sum()
	## Make same scale as original image	
	return scale*gaussian
	

##============================================================
## Translate a width in arcseconds to pixels

def AsecToPix(w, infile):
	
	
	return

##============================================================	


##============================================================
##============================================================
## Mapping kernel
##============================================================
##============================================================


## Once all GPSFs are calulated, work out kernels to and from
def all_kernels(kdim, writefile=True,makeimage=False, vb=False):
		
	## Get files -- path hard-coded
	allpngs = numpy.sort(file_seek("../Data/GPSF/","snap*d0_CFHTLS_03_g*.png"))
	
	## Collect pairs of files
	for i in range (0,len(allpngs),2):
	
		## Extract image info
		image1 = png_pix(allpngs[i])
		image2 = png_pix(allpngs[i+1])
		
		## Cut out extraneous white pixels
		image1 = pngcropwhite(image1)
		image2 = pngcropwhite(image2)
		
		## Deconvolve
		A=get_kernel(image1,image2, kdim)
		B=get_kernel(image2,image1, kdim)	

##------------------------------------------------------------
			
		## Write kernels to file and make images
		
		## Gaussian to PSF
		outfile=os.path.splitext(allpngs[i])[0][:-4]+"_GausstoPSF_"+str(A.shape[0])+".krn"
		f=open(outfile,"w")
		f.write("## Convolution kernel taking 2D Gaussian to PSFEx image\n\n\n")
		writearray(f,A,True)
		if vb: print "DeconvolveToTargetPSF.py: kernel written to",outfile
		
		## PSF to Gaussian		
		outfile=os.path.splitext(allpngs[i])[0][:-4]+"_PSFtoGauss_"+str(B.shape[0])+".krn"
		f=open(outfile,"w")
		f.write("## Convolution kernel taking PSFEx image to 2D Gaussian\n\n\n")
		writearray(f,B,True)
		if vb: print "DeconvolveToTargetPSF.py: kernel written to",outfile
		
		print "\n"

##------------------------------------------------------------	
	
	return

##============================================================	

## Find convolution kernel whch maps image1 to image2 (both arrays)
## image1 * kernel = image2
	## Kernel must be odd and should be square
	## Some information from image2 gets shaved off
def get_kernel(image1, image2, kernel_dim, vb=False):
	
##------------------------------------------------------------		

	## Enforce valid kernel shape
	if kernel_dim[1]!=kernel_dim[0]:
		kernel_dim[1] = kernel_dim[0]
		if vb: print "DeconvolveToTargetPSF.py: get_kernel: oblong kernel. Using",kernel_dim,"instead."
	if kernel_dim[0]%2==0:
		kernel_dim = [x-1 for x in kernel_dim]	## Unnecessarily Pythonic
		if vb: print "DeconvolveToTargetPSF.py: get_kernel: even kernel. Using",kernel_dim,"instead."
	
##------------------------------------------------------------		
		
	## Rearrange image1 pixels into a suitable matrix 
	orig = stacker(image1, *kernel_dim)
	## Shave image2 vector (for "honest" solution)
	## Can do this because it's overconstrained
	desired_dim = numpy.array(image1.shape)-numpy.array(kernel_dim)+1
	proc = shave(desired_dim, image2).flatten()
	
	t0=time.time()
	## Compute kernel
	kernel = LSsolve(orig,proc)
	if vb: print "DeconvolveToTargetPSF.py: get_kernel: \
								kernel computation took",round(time.time()-t0,2),"seconds."
	
	## Normalise and shape
	kernel[kernel<cutoff] = 0.0
	kernel /= kernel.sum()
	kernel = numpy.reshape(kernel, kernel_dim)
	
	return kernel
	
##============================================================
## Stacks an array of pixel values into a matrix for dotting with flattened kernel
	## Arguments are image array and kernel shape (flipped)
	
def stacker(img_arr, krn_h,krn_w):
	
	## Dimensions
	h, w = img_arr.shape
	numcols = krn_h*krn_w
	numrows = (img_arr.shape[0]-krn_h+1)*(img_arr.shape[1]-krn_w+1)	
	
	## Prepare result array
	stacked = numpy.zeros([numrows, numcols])
	
	## Loop over rows (i) in the new array
	i=0;j=0;k=0
	while k<numrows:
		## Each row in new array is a kernel-shaped slice from old
		try:
			stacked[k,:] = img_arr[j:j+krn_w,i:i+krn_h].flatten()
		except ValueError:
			print "DeconvolveToTargetPSF.py: stacker: kernel larger than PSF image! Abort."
			return stacker(img_arr, *img_arr.shape)
		i+=1
		k+=1
		## Move down a row in original and start from column 0
		if i+krn_w-1==w:
			i=0
			j+=1
			continue
	scipy.misc.imsave("stac.png",stacked)	
	return stacked
	
##============================================================
## Delete rows and columns from edges of a 2D array
## so that it has shape = eudim.
def shave(eudim, arr):
	if ((numpy.array(arr.shape)-eudim)%2).all()!=0:
		print "shave: odd shave -> over-cut."
	while arr.shape[0]>eudim[0]:
		arr=numpy.delete(arr, 0, 0)	
		arr=numpy.delete(arr, -1, 0)	
	while arr.shape[1]>eudim[1]:
		arr=numpy.delete(arr, 0, 1)	
		arr=numpy.delete(arr, -1, 1)
	return arr
	
## Add rows and columns to edges of a 2D array
## so that it has shape = eudim.
def pad(eudim, arr):
	dimdiff = (eudim-numpy.array(arr.shape))
	if (dimdiff%2).all()!=0:
		print "pad: odd pad -> over-pad."
	newarr = numpy.zeros(eudim)
	newarr[dimdiff[0]/2:-dimdiff[0]/2,dimdiff[1]/2:-dimdiff[1]/2] = arr
	return newarr
	

##============================================================
## Marix equation solvers

## Minimising square error via matrix multiplication
	## Seems to be "unstable", don't use
def puresolve(X,y):
	from numpy import dot
	t0=time.time()
	XT = numpy.transpose(X)
	ans = dot(dot(numpy.linalg.inv(dot(XT,X)),XT),y)
	print "puresolve time",round(time.time()-t0,5)
	return ans

## Solve system using QR decomposition	
def QRsolve(X,y):
	from numpy import dot
	t0=time.time()
	Q,R = numpy.linalg.qr(X)
	ans = dot(dot(numpy.linalg.inv(R),Q.T),y)
	print "QRsolve time",round(time.time()-t0,5)
	return ans

## Solve system using least squares		
def LSsolve(X,y):
	#t0=time.time()
	kvec,resi,rank,sing = numpy.linalg.lstsq(X,y)
	#print "LSsolve time",round(time.time()-t0,3),"seconds"
	return kvec

## Not working	
def sympysolve(X,y):
	t0=time.time()
	X = sympy.Matrix(X)
	#ans= X.LDLsolve(y)
	ans= X.cholesky_solve(y)
	print "time",round(time.time()-t0,5)
	return ans
 


##============================================================
##============================================================
## Total Image Deconvolution
##============================================================
##============================================================
"""
Using the kernels for transforming Gaussian PSF to observed PSF,
deconvolve an entire image to a constant PSF
"""

##============================================================

def deconvolve_image(imagefile, kernel, vb=False):
	
##------------------------------------------------------------
	## Read in image
	
	## Determine image file type and get pixels
	try:
		imgext = os.path.splitext(imagefile)[1]
		if imgext==".fits":
			imgarr = fits_pix(imagefile)
		elif imgext==".png":
			assert False, "Cannot handle PNG format."
			
		## For some reason the image is flipped at this point, so un-flip
		imgarr = imgarr[::-1,:]
	
	## Or, if the imagefile is already an array	
	except AttributeError:
		imgarr = imagefile

##------------------------------------------------------------
	## Read in kernel
	
	## Distinguish between array-kernel and file-kernel
	if type(kernel) is str:	
		## Ensure the right kernel file has been selected
		if "PSFtoGauss" in kernel:
			kernel = numpy.loadtxt(kernel)
		else:
			print "DeconvolveToTargetPSF.py: deconvolve_image: wrong kernel file. Abort."
			return
	elif type(kernel) is numpy.array:
		pass
	
	## Kernel dimensions
	kernel_h,kernel_w = kernel.shape	
	
	#scipy.misc.imsave("../doc/Images/kernel.png", kernel)				#####
	#scipy.misc.imsave("../doc/Images/image.png",imgarr)					#####
##------------------------------------------------------------
	
	## Have kernel and original image; now we need to solve the matrix equation
	## Kx = b, where K is kernel, x is target image, and b is our observed image.
	
	## FIRST, dimensions of data image should be greater than those of scene
	## to ensure a stable solution.
	## Pad so that the scene has same dimensions of original data.
	refimg_dim = imgarr.shape
	refimg_size = imgarr.size
	obsimg_dim = numpy.array(imgarr.shape) + numpy.array(kernel.shape) - 1
	obsimg_size = obsimg_dim[0]*obsimg_dim[1]	
	
	
	## SECOND, construct matrix K from kernel array.	
	## kernmat should be obsimg_size by refimg_size
			
	## What is first row of kernmat?
	krow = numpy.zeros(refimg_size)
	for i in range(kernel_h):
		j = i*refimg_dim[0]	## Stride number
		krow[j:j+kernel_w] = kernel[i,:]
	## First column of kernmat (declare too large)
	kcol = numpy.zeros(obsimg_size); kcol[0]=krow[0]
	
	## Create dense Toeplitz matrix (not square)
	## This contains all the relevant rows; we'll chop it up later
	kernmat = scipy.linalg.toeplitz(kcol,krow)	
	
	## Some rows of the kernel matrix will be no good -- they represent pixel where the convolution
	## oversteps the image border
		## There must be a better way of doing this
	badrow = [numpy.arange(i*refimg_dim[0]-kernel_w+2, i*refimg_dim[0]+1) for i in range(1,refimg_dim[1])]
	## Delete spurious rows
	kernmat = numpy.delete(kernmat, badrow, axis=0)	
	## Also delete end rows from declaring kcol too large
		##convimgsiz is how big the kernel matrix should be if we are true to convolution
	convimgsiz = (refimg_dim[0]-kernel.shape[0]+1)*(refimg_dim[1]-kernel.shape[1]+1)
	kernmat = kernmat[:convimgsiz,:]
	
	#scipy.misc.imsave("../doc/Images/kernmat.png", kernmat)				#####
	
	## Now we need to pad out the kernel matrix to reflect padding of the data image	
	## do inner bits
	insertionrow = numpy.arange(refimg_dim[0]-kernel.shape[0]+2, convimgsiz, refimg_dim[0]-kernel.shape[0]+1).repeat(2*(obsimg_dim[0]-imgarr.shape[0]))
	kernmat = numpy.insert(kernmat, insertionrow, numpy.zeros(refimg_size), axis=0)
	## do ends
	thickpad = numpy.zeros([(obsimg_dim[0]+1)*(obsimg_dim[0]-imgarr.shape[0]),refimg_size])
	kernmat = numpy.append(kernmat, thickpad, axis=0)
	kernmat = numpy.append(thickpad, kernmat, axis=0)
	
	#scipy.misc.imsave("../doc/Images/kernmat_pad.png", kernmat)				#####
	
	## Don't need these any more
	kcol=None; krow=None; badrow=None; insertionrow=None; thickpad=None
	scipy.misc.imsave("kernmat.png",kernmat)
	
	## Turn dense kernmat into a sparse matrix
	#kernmat = scipy.sparse.dia_matrix(kernmat)

	
	## THIRD, construct vector b from observed image
	## This is padded with zeros so to make the matrix equation well-defined
	imgvec = pad(obsimg_dim, imgarr).flatten()
	
	#scipy.misc.imsave("../doc/Images/image_pad.png", imgvec.reshape(obsimg_dim))				#####
		
	## FOURTH, solve the equation Kx = b
	## The system of equations is UNDERdetermined
	t0=time.time()
	
	## Decide whether to use initial guess or not
	if False:
		print "Using intital guess."
		## Make initial guess; must be conformable
		refimg0 = shave(refimg_dim, imgarr).flatten()
		## Solve matrix equation K*(dx)=res
		residual = imgvec-numpy.dot(kernmat,refimg0)
		#correction,istop,itn = scipy.sparse.linalg.lsqr(kernmat, residuals, calc_var=False)[:3]
		correction,resi = numpy.linalg.lstsq(kernmat, residual)[:2]
		## Add initial solution to correction
		refimgvec = refimg0 + correction
		## Tidy up
		correction = None; residual = None; refimg0 = None
		
	## Or just straight solve
	else:
		## Account for noise
		#####
		## Sparse method or regular
		#refimgvec,istop,itn,r1norm,r2norm = scipy.sparse.linalg.lsqr(kernmat, imgvec, calc_var=False)[:5]
		refimgvec,resi = numpy.linalg.lstsq(kernmat, imgvec)[:2]
	
	if vb:
		print "image_deconvolve: least-squares converged to reference image in",\
					round(time.time()-t0,2),"seconds."

##------------------------------------------------------------					
	## Check: reverse procedure
	check = False
	if check is True:
		resi = (numpy.dot(kernmat,refimgvec)-imgvec).reshape(obsimg_dim)
		scipy.misc.imsave("check_resi.png", resi)
		scipy.misc.imsave("check_recon.png", numpy.dot(kernmat,refimgvec).reshape(obsimg_dim))
					
##------------------------------------------------------------
	## Reshape into an image
	refimg = refimgvec.reshape(refimg_dim)
	
##------------------------------------------------------------	

	## Outfile name
	try:
		imagedir,imagename = os.path.split(imagefile)
		info = imagename[imagename.find("CFHT"):imagename.find("sci")+3]
		outfile = imagedir+"/gDeconvolved_"+info
		header = fits_hdr(imagefile)
	except (TypeError, AttributeError):
		## If the input was in array format
		outfile = "deconvolved"
		header = None
		
##------------------------------------------------------------	
	
	## Rescale (linear stretch)
	if True:
		## Scaling parameters
		vmin,vmax = stretch_params(refimg)
		## Stretch
		refimg = linear_rescale(refimg, vmin, vmax)
		## Modify filename
		outfile = outfile+"_rescaled"
	
	## Delete pre-existing images
	if os.path.isfile(outfile+".fits"): os.remove(outfile+".fits")
	if os.path.isfile(outfile+".png"):  os.remove(outfile+".png")
	
	## Save image as FITS and as PNG
	#makeFITS(outfile+".fits", refimg, header)
	#scipy.misc.imsave(outfile+".png", refimg)
	if vb: print "DeconvolveToTargetPSF.py: deconvolve_image: image saved to",outfile+".fits"
	
	return refimg
	
	
##============================================================	
##============================================================
if __name__=="__main__":
	main(argv)