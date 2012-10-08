"""
Works as a standalone module, or integrated into HumVI

This means there's some extra baggage.

If an argument variable doesn't have the correct type, e.g. expects list but gets None, the relevant steps of the code are skipped.


vb needs to be a parameter
"""

##============================================================

## Modules
import numpy
import scipy
import os
from ImageMethods import png_pix, fits_pix
from sys import argv
try:
	import io_pjm
except ImportError:
	print "wherry.py: cannot import io_pjm.py -- cannot take channel input."

##============================================================

## Global variables
from common import cutoff

##============================================================

def main(args):	
	return nw_rgb_make(*args[1:])

##============================================================
def nw_rgb_make(Rim,Gim,Bim,\
								name=None,scalefactors=None,nonlinearity=3.0,\
								yrebin=1.0,xrebin=1.0,origin=[0.,0.,0.], \
                saturatetowhite=True,overlay=False,underlay=False,invert=False):
  
	vb = False
  
	## Read in data if given file names
	if type(Rim) is str:
		imagestack = readin(Rim,Gim,Bim)
	## If given arrays, make stack
	elif type(Rim) is numpy.ndarray:
		imagestack = numpy.dstack(Rim,Gim,Bim)
	## If given channel instances from fits2colorjpg
	elif isinstance(Rim, io_pjm.channel):
		## Need to flip the images
		R = Rim.image[::-1,:]; G = Gim.image[::-1,:]; B = Bim.image[::-1,:]
		imagestack = numpy.dstack((R,G,B))
		del(R,G,B)
	else:
		print "wherry: nw_rgb_make: unknown input data type. Abort."
		return
    
  ## Outfile name
	outfile = "./WS.jpg"
	if name!=None:
		outfile = name
	
	## Subtract median if this hasn't been done by fits2colorjpg
	if isinstance(Rim, io_pjm.channel)==False:
		for k in range (3):
			imagestack[:,:,k] -= numpy.median(imagestack[:,:,k])
	
	## Rescale if this hasn't already been done
	if type(scalefactors) is list and scalefactors!=[1.0,1.0,1.0]:
		if vb: print "wherry: nw_rgb_make: rescaling by", scalefactors
		imagestack = nw_scale_rgb(imagestack, scalefactors)
	else:
		if vb: print "wherry: nw_rgb_make: no rescale."

	## Rebin
	if rebin[0] == rebin[1] == 1.0:
		if vb: print "wherry: nw_rebin_image: no rebinning."
	else:
		imagestack = nw_rebin_image(imagestack, [yrebin,xrebin])
	  
	## Kill noise -- is this allowed?
	imagestack[imagestack<cutoff] = 0.0

	## Arsinh stretch
	try:
		nonlinearity = float(nonlinearity)
		if vb: print "wherry: nw_rgb_make: arsinh stretch with nonlinearity",nonlinearity
		imagestack = nw_arsinh(imagestack, nonlinearity)
	except (ValueError, TypeError):
		if vb: print "wherry: nw_rgb_make: no stretch."

	## Box-cut
	if saturatetowhite==False:
		if vb: print "wherry: nw_rgb_make: box-cut."
		imagestack = nw_cut_to_box(imagestack, origin)

	## Overlay / underlay
		#### Not sure what this does CHECK	####	
	if overlay is not False:
		imagestack = numpy.clip(imagestack,1,overlay)
		#imagestack[imagestack>overlay]=0.0
		#imagestack[imagestack<1]=0.0
	if underlay is not False:
		imagestack = numpy.clip(imagestack,underlay,0)
		#imagestack[imagestack<underlay]=0.0
		#imagestack[imagestack>0]=0
		
	## Invert
		## Does this make sense?
	if invert==True:
		if vb: print "wherry: nw_rgb_make: inverting image."
		for i in [0,1,2]:
		 	imagestack[:,:,i] = imagestack[:,:,i].max()+imagestack[:,:,i].min()-imagestack[:,:,i]

	## Write to JPEG
	
	"""
	## If we don't trust SciPy to convert image to PNG correctly
	imagestack -= imagestack.min()
	imagestack /= imagestack.max()	
	#imagestack = numpy.clip(imagestack, 0,1)
	imagestack = numpy.clip(numpy.floor(imagestack * 256.).astype(int), 0, 255)#.astype(bytes)
	"""
	
	scipy.misc.imsave(outfile, imagestack)
	if vb: print "wherry: nw_rgb_make: saved RGB JPEG to\n",outfile,"\n"

	return None
  
  
##============================================================
##============================================================

## Readin

def readin(Rim,Gim,Bim):
 
  ## Test filetype (assuming they are all the same)
  if os.path.splitext(Rim)[1]==".png":
  	print "wherry: readin: PNG DOESN'T WORK."
  	pixel_reader = png_pix
  elif os.path.splitext(Rim)[1]==".fits":
  	pixel_reader = fits_pix  	
  else:
  	print "wherry: readin: invalid file extension. Abort."
  	return
  
  ## FITS / PNG: read in rgb images to 3-array
  imagestack = numpy.dstack(pixel_reader(fi) for fi in [Rim,Gim,Bim])
  
  ## Flip image
  imagestack = imagestack[::-1,...]
  
  ## Test whether images are the same size
  dim = imagestack[:,:,0].shape
  if dim!=imagestack[:,:,1].shape!=imagestack[:,:,2].shape:
  	print "wherry: nw_rgb_make: Images have different dimensions. Abort."
  	return
  	
  return imagestack
  	
  
##============================================================
## Rebin.

## Can take a single image or a stack
def nw_rebin_image(arr, rebin=[1.0,1.0]):
	
	## Old and new dimensions
	in_h,in_w = arr.shape[0],arr.shape[1]
	out_h,out_w = int(in_h*rebin[0]),int(in_w*rebin[1])
	
##------------------------------------------------------------
	
	## For each axis...
	for ax in [0,1]:
	
		## Oversample if we increase the size
		if rebin[ax] > 1.0:
			## If rebinfactor is an integer
			if rebin[ax]%1.0==0.0:
				arr = numpy.repeat(arr,rebin[ax],axis=ax)
			## If rebinfactor is non-integer
			else:
				print "wherry: nw_rebin_image: need something more fancy!"
				return None
				## Just find closest multiple etc...
				
		## Interpolate if we decrease the size
		## Works only if rebin is a factor -- maybe this is okay for now
		elif rebin[ax] < 1.0:
			inv = int(1/rebin[ax])
			if ax==0:
				arr = arr.reshape((out_h,inv,out_w,1)).mean(1)
			elif ax==1:
				arr = arr.reshape((out_h,1,out_w,inv)).mean(3)
	
	## Ensure correct shape
	if arr.shape[:2] != (out_h,out_w):
		print "wherry: nw_rebin_image: unexpected dimensions,",arr.shape[:2],(out_h,out_w)
	
	return arr
	
		
##============================================================

## Rescale

## Expects array stack
def nw_scale_rgb(arrstk, scalefactors=[1.,1.,1.]):
	for i in range (3):
		arrstk[:,:,i]*=scalefactors[i]
	return arrstk

##============================================================

## Arsinh stretch (note "arcsinh" is wrong).

## Expects array stack
def nw_arsinh(arrstk, nonlin=3.0):
	
	## No stretch
	if nonlin==0.0:
		print "wherry: nw_arsinh: no change made."
		return arrstk
	
	## Add up all pixels to get broadband intensity
	pixtot = arrstk.sum(axis=-1)
	pixtot[pixtot==0.0] = 1.0
	## Arsinh factor
	fac = numpy.arcsinh(pixtot*nonlin)/(pixtot*nonlin)
	
	## Multiply each array in stack	by factor
	for i in range (3):
		arrstk[:,:,i]*=fac	
	
	return arrstk

##============================================================

## Box-cut

## Pixels saturate to a specific colour rather than white.
## Expects array stack
def nw_cut_to_box(arrstk,origin=[0.0,0.0,0.0]):
	
	#pos_dist = 1.0 - numpy.array(origin)
	dim = arrstk.shape[:2]
	
	## Highest pixel-value at given position
	maxpix = numpy.max(arrstk, axis=-1)
	maxpix[maxpix<1.0]=1.0
		
	## Rescale and offest
	for k in range (3):
		arrstk[:,:,k] /= maxpix
		arrstk[:,:,k] += origin[k]
	
	return arrstk

##============================================================

if __name__=="__main__":
	main(argv)
