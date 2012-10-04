#!/usr/bin/env python
## =====================================================================

## Globally useful modules:

import numpy
import sys,getopt,Image
import humvi

## =====================================================================

def compose(argv):
	"""
	NAME
	  compose.py

	PURPOSE
	  Make a color composite JPG image from 3 FITS images, using the Lupton
	  algorithm.

	COMMENTS
	  Reads filter name and zero point from header and tries to set scales 
	  automatically - unless scales are set on command line. Note 
	  telescope/survey has to be recognised for this to work...

	USAGE
	  compose.py [flags] red.fits green.fits blue.fits

	FLAGS
	  -h            Print this message
	  -v            Verbose operation

	INPUTS
	  red.fits etc  Names of FITS files containing image data

	OPTIONAL INPUTS
	  -s --scales   3,2,4     Comma-separated scales for R,G,B channels [None]
	  -p --parameters  5,0.05 Non-linearity parameters Q,alpha
	  -o --output   outfile   Name of output filename [guessed]

	OUTPUTS
	  stdout        Useful information
	  outfile       Output plot in jpg or png format


	EXAMPLES

	  compose.py  CFHTLS_27_i_sci.fits CFHTLS_27_r_sci.fits CFHTLS_27_g_sci.fits
	                    
	BUGS

	HISTORY
	  2012-05-11    started Marshall (Oxford)
	  2012-07-??    integrated with wherry.py Sandford (NYU)
	"""

	## -------------------------------------------------------------------

	try:
	    opts, args = getopt.getopt(argv, "hvs:p:n:o:l",\
	    									["help","verbose","png","scales","pars","output","lupton"])
	except getopt.GetoptError, err:
	    ## print help information and exit:
	    print str(err) ## will print something like "option -a not recognized"
	    print compose.__doc__
	    return

	vb = False

	png = False
	outfile = 'output.png'

	pars = '10,0.04'
	scales = 'Auto'

	LuptonStretch=False

	for o,a in opts:
	    if o in ("-h", "--help"):
	        print compose.__doc__
	        return
	    elif o in ("-v", "--verbose"):
	        vb = True
	    elif o in ("-p","--parameters"):
	        pars = a
	    elif o in ("-s","--scales"):
	        scales = a
	    elif o in ("-o","--output"):
	        outfile = a
	    elif o in ("-l","--lupton"):
	        LuptonStretch=True
	        if vb:	print "Lupton method selected."
	    else:
	        assert False, "Unhandled option"

	## Check for datafiles in array args:

	if len(args) == 3:
	  rfile = args[0]
	  gfile = args[1]
	  bfile = args[2]
	  if vb:
	    print "Making color composite image of data in following files:",rfile,gfile,bfile
	    print "Output will be written to",outfile

	else:
		print compose.__doc__
		return
	  
	## Parse nonlinearity parameters:  
	Qs,alphas = pars.split(',')
	Q = float(Qs)
	alpha = float(alphas)

	if scales != 'Auto':
	  x,y,z = scales.split(',')
	  rscale = float(x)
	  gscale = float(y)
	  bscale = float(z)
	  
	## -------------------------------------------------------------------
	## Read in images, set and apply scales etc:

	red = humvi.channel(rfile)
	green = humvi.channel(gfile)
	blue = humvi.channel(bfile)
		
	## Check shapes are equal:
	humvi.check_image_shapes(red.image,green.image,blue.image)

	## Subtract backgrounds (median):
	red.subtract_background()
	green.subtract_background()
	blue.subtract_background()

	## Set scales:
	if scales == 'Auto':
	  red.set_scale()
	  green.set_scale()
	  blue.set_scale()
	  rscale,gscale,bscale = humvi.normalize_scales(red.scale,green.scale,blue.scale)

	red.set_scale(manually=rscale)
	green.set_scale(manually=gscale)
	blue.set_scale(manually=bscale)
	  
	if vb:	print 'Scales set to:',red.scale,green.scale,blue.scale

	# Scale images - only do once:
	red.apply_scale()
	green.apply_scale()
	blue.apply_scale()

	## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	## Lupton scaling
	
	if LuptonStretch:

	  if vb:	print "Nonlinearity parameters Q,alpha:",Q,alpha

	  I = humvi.lupton_intensity(red.image,green.image,blue.image)
	  
	  r = humvi.lupton_stretch(red.image,I,Q,alpha)
	  g = humvi.lupton_stretch(green.image,I,Q,alpha)
	  b = humvi.lupton_stretch(blue.image,I,Q,alpha)

	  ## Package into a python Image, and write out to file:

	  image = humvi.pack_up(r,g,b) 
	  image.save(outfile)
	  if vb: print "Image saved to:",outfile

		return

	## -------------------------------------------------------------------
	## Wherry scaling -- default
	
	else:
	  
		## Arguments are (R,G,B, outfilename,scale,nonlinearity,\
		## yrebin,xrebin,origin,saturatetowhite,overlay,underlay,invert)
				## Note that scaling has been done already so scale should be None
		humvi.nw_rgb_make(red,green,blue, outfile, None, Q*alpha)

		return	## Redundant but keep for clarity
	
	
# ======================================================================

if __name__ == '__main__':
  compose(sys.argv[1:])

# ======================================================================
