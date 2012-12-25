#!/usr/bin/env/python
##============================================================

from decon import SourceDiagnostics as SD
from decon import DeconvolveToTargetPSF as DT
import ImageMethods as IM
from OSMethods import file_seek
import time
import getopt
import os
import scipy
from sys import argv

	
def main(argv):
	"""
	NAME
	  deconvolve.py

	PURPOSE
		Deconvolution pipeline. Procedure is essentially
		1. Identify sources in image.
		2. Estimate PSF for image.
		3. Generate target PSF.
		4. Find convolution kernel fot this target.
		5. Deconvolve original image with this kernel.

	COMMENTS
	  
	USAGE
		python deconvolve.py [flags] filename/directory

	FLAGS
		-h            Print this message
		-vb						Verbose
	  
		-p	--psfw		f		PSF size in arcseconds
		-k	--kern		i		Deconvolution kernel size in pixels
		-s	--instr		s		Select files containing the following string
				--sim					Deconvolve all inputs to the same resolution

	INPUTS

	OPTIONAL INPUTS

	OUTPUTS
	  stdout        Useful information
	  
	  image catalogue
	  .psf
	  ... MORE
	  psf snapshot
	  

	EXAMPLES
		py deconvolve.py ../Data/test_03g/
		
		py deconvolve.py --sim ../Data/test03/CFHTLS_03_i_sci.fits ../Data/test03/CFHTLS_03_r_sci.fits  ../Data/test03/CFHTLS_03_g_sci.fits
	                
	                    
	BUGS / ISSUES		
		- Need an elegant way of combining the single and simultaneous
		deconvolutions		
		- Need a better way of getting the PSF widths -- PSFEx must tell us!

	HISTORY
	  2012-09 started Sandford (NYU)
	"""
	
	t_init = time.time()
	
	##-------------------------------------------------------------------
	## Options and help-string
	
	try:
		opts, args = getopt.getopt(argv[1:], "hvp:k:s:",\
																		["help","psfw","kern","fstr","sim"])
	except getopt.GetoptError, err:
		## Print help information and exit:
		print str(err)
		print main.__doc__
		return
	
	## Default flag settings
	vb=False
	simul=False
	psf_width=None
	file_string="CFHTLS*sci.fits"
	krn_size=15
		
	## Search for flags		   
	for o,a in opts:
		if o in ("-h", "--help"):
			print main.__doc__
			return
		elif o in ("-v"):
			vb=True
		elif o in ("-p", "--psfw"):
			psf_width = float(a)
		elif o in ("-k", "--kern"):
			krn_size = int(a)
		elif o in ("-s", "--fstr"):
			file_string = a
		elif o in ("--sim "):
			simul=True
		else:
			assert False, "Unhandled option."
			
	##--------------------------------------------------------------------
	
	## Is input a collection of files or directory?
	## filelist is a list of raw image files.
	filelist=[]
	for arg in args:
		if os.path.isfile(arg):  filelist+=[arg]
		elif os.path.isdir(arg): filelist+=file_seek(arg,file_string)
	
	if simul is False:
		single_deconvolve(filelist, krn_size, psf_width, vb)
	else:
		simultaneous_deconvolve(filelist, krn_size, vb)
	
	##--------------------------------------------------------------------
	
	if vb: print "deconvolve.py: deconvolution of",len(filelist),\
								"file(s) took",round(time.time()-t_init,3),"seconds."
	
	return

##============================================================
	
def single_deconvolve(filelist, krn_size, psf_width, vb):

	## Deconvolve all FITS images in list
	for infile in filelist:
	
		## Run SExtractor to generate source-catalogue
		catfile = SD.SEx(infile, vb)
		## Run PSFEx to generate model image of PSF
		psf_obs = SD.PSFEx(catfile, "snap", vb)
		
		##--------------------------------------------------------------------
		
		## Intermediate step: convert PSF image to an array
		psf_obs = IM.pngcropwhite(IM.fits_pix(psf_obs))
		
		##--------------------------------------------------------------------		
		
		## Make target psf array
		gaussparams = DT.moments(psf_obs)
		if type(psf_width) is float:			
			gaussparams[1]=gaussparams[2]=psf_width#AsecToPix(psf_width,infile)
		psf_ref = DT.Gauss_2D(*gaussparams)
		#scipy.misc.imsave("./ref_"+str(psf_width)+".png",psf_ref)
		
		##--------------------------------------------------------------------		
		
		## Calculate deconvolution-kernel array
			## Kernel size = psf width -- assume this is ~ correlation length-scale
		krn_size = 3*int(max(gaussparams[1:3]))	### ONLY APPROPRIATE when using simdec
		kernarr = DT.get_kernel(psf_obs, psf_ref,[krn_size,krn_size],vb)	
		## Deconvolve image
		DT.deconvolve_image(infile, kernarr, vb)
	
	return
	
##============================================================

### WORKING ON THIS ### NEEDS UPDATING
def simultaneous_deconvolve(filelist, krn_size, vb):
	
	s2_max = 0.0
	
	## Find largest PSF in ensemble
	for infile in filelist:
	
		## Run SExtractor to generate source-catalogue
		catfile = SD.SEx(infile, vb)
		## Run PSFEx to generate model image of PSF
		psf_obs = SD.PSFEx(catfile, "snap", vb)
		
		### NEED to find FWHM from images without going through the
		### following rigmarole:
		
		## Convert PSF image to an array
		psf_obs = IM.pngcropwhite(IM.fits_pix(psf_obs))
		
		## Here is how we decide what PSF to use
		## Product of the two widths = s2
		m=DT.moments(psf_obs)
		s2 = m[1]*m[2]
		## Which image has largest s2? Call this psf_OBS.
		if s2 > s2_max:
			s2_max = s2
			psf_OBS = psf_obs

##------------------------------------------------------------
	
	## Make target psf array
	psf_ref = DT.Gauss_2D(*DT.moments(psf_OBS))
	## Calculate kernel array for obs->ref
	kernarr = DT.get_kernel(psf_OBS, psf_ref, [krn_size,krn_size],vb)

##------------------------------------------------------------	

	## Deconvolve to the target	PSF
	for infile in filelist:
		DT.deconvolve_image(infile, kernarr, vb)	
	
	return



##============================================================
if __name__=="__main__":
	main(argv)
