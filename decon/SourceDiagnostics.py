#!/usr/bin/env/python
##============================================================


import os
import numpy
import time
from sys import argv
from OSMethods import file_seek

##============================================================

## Seek out all fits image files in directory and make psfs
	## argv[1] is the directory to look at
	## argv[2] is option for PSFEx (maybe later)
	
def main(args):
	"""
	NAME
	  SourceDiagnostics.py

	PURPOSE
		Identifies sources in a FITS image, and calculates the image.

	COMMENTS
		Uses SExtractor and PSFEx software.
	  
	USAGE
		python SourceDiagnostics.py [flags] directory/file

	FLAGS
	  -h            Print this message
	  -vb						Verbose

	INPUTS
		Requires the following configuration files:
		- default.conv
		- default.nnw
		- prepsfex.param
		- prepsfex.sex
		- default.psfex
		
		Input directory or filename

	OPTIONAL INPUTS

	OUTPUTS
	  stdout        Useful information	  
	  image.cat file	Diagnostics about each source
	  snap_image.fits	PSF model
	  image.psf				More sophisticated PSF model

	EXAMPLES
		py SourceDiagnostics.py CFHTLS_03_g_sci.fits
		py SourceDiagnostics.py CFHTLS_03/
		
	                    
	BUGS / ISSUES
		- Hard coded paths to configuration files

	HISTORY
	  2012-06 started Sandford (NYU)
	"""
	
	## -------------------------------------------------------------------
	## Options and help-string
	
	try:
		opts, args = getopt.getopt(argv, "hv",["help"])
	except getopt.GetoptError, err:
		## Print help information and exit:
		print str(err)
		print main.__doc__
		return
	
	vb=False
			   
	for o,a in opts:
		print o,a
		if o in ("-h", "--help"):
			print main.__doc__
			return
		elif o in ("-v"):
			vb=True
		else:
			assert False, "Unhandled option"
			
	##--------------------------------------------------------------------
	

	t0=time.time()
	
	## If the input is a single file
	if os.path.isfile(args[0]):	fits_list = [args[0]]
	## If the input is a directory of files
	elif os.path.isdir(args[0]):	fits_list = file_seek(args[0], "CFHTLS_*sci.fits")
	
	## Run SExtractor and PSFEx on the image
	for fitsfile in fits_list:
		## Run SExtractor
		SEx(fitsfile, vb)		
		## Run PSFEx on fresh .cat
		PSFEx(fitsfile[:-4]+"cat", "snap", vb)
	
	
	if vb:
		print "\nSourceDiagnostics.py: analysed",len(fits_list),"FITS images."
		print "\tTook",round(time.time()-t0,2),"seconds."
		
	return
	
##============================================================

## DEFUNCT
## Runs SExtractor or PSFEx depending on input file
	## argv[1] is a FILE
def choose(argv):
	
	if argv[1][-4:]=="fits": SEx(argv)
	elif argv[1][-3:]=="cat": PSFEx(argv)
	else: print "SourceDiagnostics.py: provide suitable input file."
	
	return

##============================================================

def SEx(imgname, vb):
	
	## Make absolute path
	imgname = os.path.abspath(imgname)
	
	## Check the input
	if imgname[-5:]!=".fits":
		print "SourceDiagnostics.py: SEx: invalid input image."
		return
		
	## PSFEx needs output in FITS_LDAC format, and different data in parameter file
	confile = " -c ~/Dropbox/Hogg_2012/HumVI/decon/AstromaticConfig/prepsfex.sex "
	confile += " -FILTER_NAME ~/Dropbox/Hogg_2012/HumVI/decon/AstromaticConfig/default.conv "
	params  = " -PARAMETERS_NAME ~/Dropbox/Hogg_2012/HumVI/decon/AstromaticConfig/prepsfex.param "
	cattype = " -CATALOG_TYPE FITS_LDAC "
	catname = imgname[:-5]+".cat"
	catdecl = " -CATALOG_NAME "+catname
	
	## Put it all together
	commands = "sex "+imgname+confile+params+cattype+catdecl
	if vb: commands += " -VERBOSE_TYPE NORMAL "
		
	## Execute commands
	os.system(commands)
	
	if vb: print "SourceDiagnostics.py: SEx: saved source catalogue as "+catname
	
	return catname

##============================================================

## Modify arguments of PSFEx

def PSFEx(catfile, mode, vb):

	## Check the input
	if catfile[-4:]!=".cat":
		print "SourceDiagnostics.py: PSFEx: invalid input catalogue. Abort."
		return
	
	## Input catalogue
	outdir,catname = os.path.split(catfile)
	## Destination directory
	outdir = outdir+"/Diagnostics/"
	if os.path.isdir(outdir)==False: os.mkdir(outdir)

##------------------------------------------------------------
	imgname = None
	config = " -c ~/Dropbox/Hogg_2012/HumVI/decon/AstromaticConfig/default.psfex "
	
	## If we want to return only a PSF snapshot
	if mode=="snap":
	
		## Only one check image: snap_n1d0*.fits
		imgname = outdir+"snap_n1d0_"+catname[:-4]+".fits"
		imgdest = " -CHECKIMAGE_TYPE SNAPSHOTS  -CHECKIMAGE_NAME "+outdir+"snap_n1d0"
		params  = " -PSFVAR_NSNAP 1  -PSFVAR_DEGREES 0  "
		## No check plots (png)
		pltdest = " -CHECKPLOT_TYPE NONE "
		## No XML file
		xmldest = " -WRITE_XML N "

	elif mode=="none":
		## No check image
		imgdest = " -CHECKIMAGE_TYPE NONE "
		## No check plots (png)
		pltdest = " -CHECKPLOT_TYPE NONE "
		## No XML file (xml)
		xmldest = " -WRITE_XML N "
		## Default 1-degree fit	
		params  = " -PSFVAR_DEGREES 1   BASIS_TYPE PIXEL_AUTO "
			
	## Compute everything -- takes longer
			## Need to re-check this stuff
	elif mode=="all":
		## Destination for check images (fits)
		imgdest = outdir+"chi.fits,"+outdir+"proto.fits,"+outdir+"samp.fits,"+\
							outdir+"resi.fits,"+outdir+"snap_n1d0.fits,"+\
							outdir+"moffat.fits,"+outdir+"submoffat.fits,"+outdir+"subsym.fits"
		imgdest = " -CHECKIMAGE_NAME "+imgdest
		## Destination for check plots (png)
		pltdest = outdir+"fwhm,"+outdir+"ellipticity,"+outdir+"counts,"+\
							outdir+"countfrac,"+outdir+"chi2,"+outdir+"resi"
		pltdest = " -CHECKPLOT_NAME "+pltdest
		## Destination for XML file (xml)
		xmldest = " -XML_NAME "+outdir+catname+".xml "
		params  = " "
		
	else:
		print "SourceDiagnostics.py: invalid PSFEx option specified. Abort."
		return
	
##------------------------------------------------------------
	
	## Stick it all together
	commands = "psfex "+catfile+imgdest+pltdest+xmldest+params+config
	if vb: commands += " -VERBOSE_TYPE NORMAL "
	
	## Execute commands
	os.system(commands)
	
	if vb:
		if mode is "snap":
			print "SourceDiagnostics.py: PSF model image saved to",imgname
	
	return imgname

##============================================================

## DEFUNCT
## Generate snapshots of varying quality and coverage
	## Called from python line (use OSScript.file_seek first)
def snapall(catfile):
	for nsnap in [1,9]:
		for deg in [0,1,2,3]:
			PSFEx(["dummy",catfile,"snap",str(nsnap),str(deg)])
	return

	
##============================================================

if __name__=="__main__":
	main(argv[1:])
