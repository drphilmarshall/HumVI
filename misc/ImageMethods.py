#!/usr/bin/env/python
##============================================================

import os
import glob
import numpy
import scipy
import pylab
import Image
import ImageEnhance
import pyfits
import getopt
from sys import argv
from OSMethods import file_seek

##============================================================

def main():
	"""
	NAME
	  ImageMethods.py

	PURPOSE
	  A variety of image manipulation methods -- filetype conversions,
	  cropping function, image anhancement, and stretching.

	COMMENTS
	  

	USAGE
	  ImageMethods.py [flags] function args
	  
	  * args depend on function

	FLAGS
	  -h            Print this message

	INPUTS
	  Directory and/or file name depending on function

	OPTIONAL INPUTS

	OUTPUTS
	  stdout        Useful information
	  outfile       When converting/modifying images
	  array					When modifying image data

	EXAMPLES
	  ImageMethods.py F2P ../Images smileyimage?
	  ImageMethods.py cropwhite ./whiteborderimg.png
	                    
	BUGS
		No bugs per se, but a lot of usability improvements / tidying up
		to do.

	HISTORY
	  2012-06-12 started Sandford (NYU)
	"""
	
	## -------------------------------------------------------------------
	## Options and help-string
	
	try:
		opts, args = getopt.getopt(argv, "h",["help"])
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
	        
	## -------------------------------------------------------------------	
	
	if args[1]=="F2P": FITStoPNG(args)
	elif args[1]=="P2P": PSFtoPNG_all(args)
	elif args[1]=="cropwhite": pngcropwhite(args[2])
	elif args[1]=="brighten": bright_all()
	elif args[1]=="stretch": stretch(args[2])	
	elif args[1]=="stretchall": stretch_all(0)	
	elif args[1]=="Stretchall": stretch_all(1)
	else: print "ImageScript.py: select valid option:\
							\n\tF2P\n\tP2P\n\tbrighten\n\tstretch\n\tcropwhite"
	
	return
	
##============================================================

##============================================================
##============================================================
## IMAGE CONVERSION
##============================================================
##============================================================

## Converts FITS images in subdirectories of argv[2] to PNG format
  ## argv[3] is string (e.g. "CFHTLS_")
def FITStoPNG(args):
  
  i=0 
  for f_img in file_seek(args[2], args[3]+"*.fits"):
  
		## Prepare destination directory
		outdir=os.path.dirname(f_img)#+"/Enhanced_Original/" ## Needs to be commented when brightening
		if os.path.isdir(outdir)==False: os.mkdir(outdir)
		
		## Outfile name
		p_img = outdir+"/"+os.path.basename(f_img)[:-5]+".png"
		
		## File conversion
		os.system("ds9 "+f_img+" -minmax -colorbar no -saveimage png "+p_img+" -exit")
		i+=1
	
	## Shave whitespace?
	##
	##
		
  print "ImageScript.py: FITStoPNG: Converted",i,"files to png format."
  return

##============================================================

## Converts all .psf to .png in subdirectories of argv[2]
def PSFtoPNG_all(args):
	for psffile in file_seek(args[2], "*.psf"):
		PSFtoPNG(psffile)		
	return

##============================================================

## Convert .psf to a number of .png
	
def PSFtoPNG(imag_in):
	
	## Prepare directory for output
	out_dir=os.path.dirname(imag_in)+"/PSF_PNG_Breakdown/"
	if os.path.isdir(out_dir)==False: os.mkdir(out_dir)
	
	## Load image
	IMG = pyfits.open(imag_in)
	
	## Extract image data
	psfdata = IMG[1].data[0][0]
	Nim = IMG[1].header["psfaxis3"]
	Xdim = IMG[1].header["psfaxis1"]
	Ydim = IMG[1].header["psfaxis2"]
	
	## Reshape to separate images
	psfdata = numpy.reshape(psfdata, [Nim,Xdim,Ydim])#, dtype="unit8"
	
	## Load data into new image to save individually
	for i in range (Nim):
		## Outfile name
		imag_out = out_dir+os.path.basename(imag_in)[:-4]+"_fit"+str(i)+".png"
		#image = Image.new("L",[Xdim,Ydim],None)
		#image.putdata(psfdata[i,:,:],1.0)
		#image = Image.frombuffer('L',(Xdim,Ydim),psfdata,'raw','L',0,1)
		#image.save(imag_out)
		scipy.misc.imsave(imag_out, psfdata[i,:,:])
	
	print "ImageScript.py: PSFtoPNG: made",Nim,"PNG files from",os.path.basename(imag_in)
	
	return
		
##============================================================	
	
## Convert from .ext to .jpg
	## Call from Python
def EXTtoJPG(infile):
  f,e=os.path.splitext(infile)
  outfile=f+".jpg"
  if infile!=outfile:
    try:
      Image.open(infile).save(outfile)
    except IOError:
      print "ImageScript.py: EXTtoJPG: cannot convert", infile
  return

##============================================================

## Crop white pixels from image border
	## Expects a pixel-array with 255-pixels around central image
	
def pngcropwhite(imgarr, level=254.9):
	
	## Check whether imgarr is filename or array
	if isinstance(imgarr, numpy.ndarray):
		returntype="array"
	elif isinstance(imgarr, str):
		filename = imgarr
		imgarr = png_pix(imgarr)
		returntype="pngfile"
	
	## Original dimensions and centre
	img_h,img_w = imgarr.shape
	crow,ccol = int(0.5*img_h), int(0.5*img_w)
	
	## Find dimensions of the central image by stepping inwards from periphery
	for i in range (img_h):
		if imgarr[i,ccol]<level:
			top=i
			break
	for i in range (img_h-1,0,-1):
		if imgarr[i,ccol]<level:
			bottom=i+1
			break
	for i in range (img_w):
		if imgarr[crow,i]<level:
			left=i
			break
	for i in range (img_w-1,0,-1):
		if imgarr[crow,i]<level:
			right=i+1
			break
	
	## Retain slice of interest
	imgarr = imgarr[top:bottom, left:right]
	
	## Output depends on input -- do we want an array or a file?
	if returntype is "array":
		return imgarr[top:bottom, left:right]
	elif returntype is "pngfile":
		scipy.misc.imsave(filename, imgarr)
		return
		

##============================================================
## Image to array conversion

## Read in .fits and return array of (decimal) pixel values
	## Works for PSFEx snap - other fits have different heirarchy
def fits_pix(filename):	
	im = pyfits.open(filename)
	return numpy.array(im[0].data)

def fits_hdr(filename):
	im = pyfits.open(filename)
	return im[0].header

## Convert png file to an array of (0-255) pixel values
def png_pix(filename):
	img = Image.open(filename)
	dim = [img.size[1],img.size[0]]	## Because of course they need to be reversed
	img = img.convert("L")
	dat = numpy.array(img.getdata())
	dat = numpy.reshape(dat, dim)
	return dat


##============================================================

## Make FITS image from array of pixels
def makeFITS(outfile, array, head=None):
	hdu = pyfits.PrimaryHDU(array[::-1])
	if head is not None:
		hdu.header = head
	hdu.writeto(outfile)
	return

	
##============================================================
## Pixel file to PNG
## Greyscale

def pixfiletopng(inname):
	DAT = numpy.loadtxt(inname)
	outname = os.path.splitext(inname)[0]+".png"
	scipy.misc.imsave(outname,DAT)
	return

##============================================================
##============================================================
## IMAGE ENHANCEMENT    
##============================================================
##============================================================

## Automate image stretching
	## sel=0: Stretch original images only
	## sel=1: Stretch original images, and apply the same stretch to diagnostic images
	
def stretch_all(sel=0):

	for string in ["CFHTLS_*_sci.fits"]:#
		for fits in file_seek("../Data/",string):
			if sel==0: stretch(fits)
			else: stretch_plus(fits)
			
	return
	
##============================================================	

## Stretches a file and then its children

def stretch_plus(fits_original):
	
	## Stretching parameters	
	vmin,vmax=stretch(fits_original)
	
	## Apply to diagnostic files
	homedir=os.path.dirname(fits_original)
	for string in ["snap_*_sci.fits","samp_*_sci.fits","resi_*_sci.fits"]:#"CFHTLS_*_sci.fits"
		for fits in file_seek(homedir,string):
			stretch(fits, vmin,vmax)
			
	return
	
##============================================================	
	
def stretch(fits_original, vmin=None,vmax=None):
	
	## Prepare destination directory
	outdir=os.path.dirname(fits_original)+"/"#Enhanced_Original/"
	if os.path.isdir(outdir)==False: os.mkdir(outdir)
	
	## Outfile name
	outfile = outdir+os.path.basename(fits_original)[:-5]+"_rescaled.fits"

	## Load image
	IMG = pyfits.open(fits_original)	
	## Array of pixel values
	im_data = IMG[0].data	
	
	## If we haven't provided stretch parameters, calculate them
	if vmin==None or vmax==None: vmin,vmax = stretch_params(im_data)
	
	## Rescale
	im_data = linear_rescale(im_data, vmin,vmax)

	## Update FITS data
	IMG[0].data = im_data
	## Delete pre-existing files
	if os.path.isfile(outfile): os.remove(outfile)
	## Save
	IMG.writeto(outfile)
	IMG.close()
	    
	return vmin,vmax
	


##============================================================
## Linear scaling of input array
def linear_rescale(inputArray, scale_min=None, scale_max=None):
	
	imageData=numpy.array(inputArray, copy=True)
	
	if scale_min == None:
		scale_min = imageData.min()
	if scale_max == None:
		scale_max = imageData.max()

	imageData = imageData.clip(min=scale_min, max=scale_max)
	imageData = (imageData-scale_min)/(scale_max-scale_min)
	imageData[imageData<0.] = 0.0
	imageData[imageData>1.] = 1.0
	
	return imageData


##============================================================

## Parameters for re-scale/stretch
def stretch_params(array):
	## Find median and 97.5% value
	a=numpy.median(array)
	b=numpy.sort(array.flatten())[0.975*len(array.flatten())]
	## Stretch parameters
	sigma = 0.5*(b-a)
	vmin = a-5.0*sigma
	vmax = a+10.0*sigma
	return vmin,vmax


##============================================================ 
## Standard enhancements

## Contrast
def cont(imagefile, factor=10):	
	im=Image.open(imagefile)
	enhancer = ImageEnhance.Contrast(im)
	enhancer.enhance(factor).show("Contrast %f" % factor)
	return
## Brightness
def brig(imagefile,factor=10,brightfile=None):
	im=Image.open(imagefile)
	enhancer = ImageEnhance.Brightness(im)
	if brightfile==None: brightfile=imagefile[:-4]+"_br"+str(factor)+".png"
	else: brightfile=os.path.dirname(imagefile)+"/Enhanced_Original/"+os.path.basename(imagefile)[:-4]+"_br"+str(factor)+".png"
	enhancer.enhance(factor).save(brightfile)
	return

##============================================================

## Brighten ALL originals
	
def bright_all(factor=8):
	for original in file_seek("../Data/","CFHTLS_*sci.png"):
		brig(original,factor,1)
	return


##============================================================
##============================================================
## Miscellanious
##============================================================
##============================================================

## Write 2D array to a file in easy-parsed fashion
def writearray(f,arr,close=True,dp=5):
	for i in range (arr.shape[0]):
		for j in range (arr.shape[1]):
			f.write(str(round(arr[i,j],dp))+"\t")
		f.write("\n")
	if close == True: f.close()
	return

##============================================================
if __name__=="__main__":
  main()
