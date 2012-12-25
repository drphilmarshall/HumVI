#!/usr/bin/env python
## =====================================================================

## Globally useful modules:

import numpy
import sys,getopt,Image
import humvi
import os

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
        -b --subtract-background 
        -l            Follow Lupton et al [default]
        -w            Follow Wherry et al

      INPUTS
        red.fits etc  Names of FITS files containing image data

      OPTIONAL INPUTS
        -s --scales   3,2,4     Comma-separated scales for R,G,B channels [None]
        -p --parameters  5,0.05 Non-linearity parameters Q,alpha
        -o --output   outfile   Name of output filename [guessed]
        -x --saturate-to style  Saturation method, color or white [white]
        -z --offset 0.1         Offset image level (+ve = gray background)    
        -m --mask -1.0          Mask images below level [don't]

      OUTPUTS
        stdout        Useful information
        outfile       Output plot in jpg or png format


      EXAMPLES

        compose.py  -v -s 0.8,1.0,1.0 -z 0.0 -p 1.0,0.03 -m -1.0 \
          -o examples/CFHTLS_27_gri.png \
          examples/CFHTLS_27_i_sci.fits \
          examples/CFHTLS_27_r_sci.fits \
          examples/CFHTLS_27_g_sci.fits

        py compose.py ../Data/test03/gDeconvolved_CFHTLS_03_i_sci_rescaled.fits ../Data/test03/gDeconvolved_CFHTLS_03_r_sci_rescaled.fits ../Data/test03/gDeconvolved_CFHTLS_03_g_sci_rescaled.fits


      BUGS
      
      - Should cope with just two or even one input image file.

      HISTORY
        2012-05-11    started Marshall (Oxford)
        2012-07-??    integrated with wherry.py Sandford (NYU)
        2012-12-19    defaults set for CFHTLS images Marshall (Adler)
      """

      ## -------------------------------------------------------------------

      try:
          opts, args = getopt.getopt(argv, "hvs:p:n:o:x:z:blwm:",\
            ["help","verbose","png","scales","pars","output","saturate-to","offset","subtract-background","lupton","wherry","mask"])
      except getopt.GetoptError, err:
          ## print help information and exit:
          print str(err) ## will print something like "option -a not recognized"
          print compose.__doc__
          return

      vb = False

      png = False
      outfile = os.path.dirname(args[0])+"/composed"

      pars = '1,0.03'
      scales = 'Auto'
      backsub = False

      LuptonStretch = True
      saturation = 'white'
      offset = 0.0
      mask = False
      masklevel = -1.0

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
          elif o in ("-x","--saturate-to"):
              saturation = a
          elif o in ("-z","--offset"):
              offset = float(a)
          elif o in ("-m","--mask"):
              mask = True
              masklevel = float(a)
          elif o in ("-o","--output"):
              outfile = a
          elif o in ("-b","--subtract-background"):
              backsub = True
          elif o in ("-l","--lupton") or LuptonStretch is True:
              LuptonStretch = True
							outfile += "_lupton"
							if vb:	print "Lupton's method selected."
          elif o in ("-w","--wherry") or LuptonStretch is False:
              LuptonStretch = False
							outfile += "_wherry"
							if vb:	print "Wherry's method selected."
          else:
              assert False, "Unhandled option"

			## Add info to outfile name
			outfile += "_"+pars+"_"+scales+".png"
			
			
      ## Check for datafiles in array args:

      if len(args) == 3:
        rfile = args[0]
        gfile = args[1]
        bfile = args[2]
        if vb:
          print "Making color composite image of data in following files:",rfile,gfile,bfile
          print "Output will be written to",outfile
          if mask: print "Masking stretched pixel values less than",masklevel

      else:
            print compose.__doc__
            return

      ## Parse nonlinearity parameters:
      Qs,alphas = pars.split(',')
      Q = float(Qs)
      alpha = float(alphas)

      ## -------------------------------------------------------------------
      ## Read in images, set and apply scales etc:

      red = humvi.channel(rfile)
      green = humvi.channel(gfile)
      blue = humvi.channel(bfile)

      ## Check shapes are equal:
      humvi.check_image_shapes(red.image,green.image,blue.image)

      # Subtract backgrounds (median, optional):
      if backsub:
        red.subtract_background()
        green.subtract_background()
        blue.subtract_background()

      # Set scales:
      if scales == 'Auto':
        red.set_scale()
        green.set_scale()
        blue.set_scale()
        rscale,gscale,bscale = humvi.normalize_scales(red.scale,green.scale,blue.scale)
      else:
        x,y,z = scales.split(',')
        rscale = float(x)
        gscale = float(y)
        bscale = float(z)
        rscale,gscale,bscale = humvi.normalize_scales(rscale,gscale,bscale)
      red.set_scale(manually=rscale)
      green.set_scale(manually=gscale)
      blue.set_scale(manually=bscale)

      if vb: print 'Scales normalized to:',red.scale,green.scale,blue.scale

      # Scale images - only do once:
      red.apply_scale()
      green.apply_scale()
      blue.apply_scale()

      ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ## Lupton code (default):

      if LuptonStretch:

            if vb: 
              print "Stretch parameters Q,alpha:",Q,alpha
              print "At low surface brightness levels, the channel images are further rescaled by alpha"
              print "Nonlinearity sets in at about 1/Q*alpha in the scaled intensity image:",1.0/(Q*alpha)
               
            # Compute total intensity image and the arcsinh of it:
            I = humvi.lupton_intensity(red.image,green.image,blue.image,type='sum')
            stretch = humvi.lupton_stretch(I,Q,alpha)

            # Apply stretch to channel images:
            r = stretch * red.image
            g = stretch * green.image
            b = stretch * blue.image

            if mask:
              # Mask problem areas - exact zeros or very negative patches should 
              # be set to zero.
              r,g,b = humvi.pjm_mask(r,g,b,masklevel)

            # Offset the stretched images to make zero level appear dark gray.
            # Negative offset makes background more black...
            r,g,b = humvi.pjm_offset(r,g,b,offset)
            
            # This should really be a snap to box option...
            if saturation == 'color':            
              # Saturate to colour at some level - might as well be 1, since
              # Q redefines scale?:
              threshold = 1.0
              r,g,b = humvi.lupton_saturate(r,g,b,threshold)
            # Otheriwse saturate to white.
            
            # Package into a python Image, and write out to file:
            image = humvi.pack_up(r,g,b)
            image.save(outfile)
            if vb: print "Image saved to:",outfile
            
      ## -------------------------------------------------------------------
      ## Wherry code:

      else:

            ## Arguments are (R,G,B, outfilename,scale,nonlinearity,\
            ## yrebin,xrebin,origin,saturatetowhite,overlay,underlay,invert)
                        ## Note that scaling has been done already so scale should be None
            humvi.nw_rgb_make(red,green,blue, outfile, None, Q*alpha)

## ======================================================================

			if vb: print "Image saved to:",outfile
			return
			
## ======================================================================

if __name__ == '__main__':
  compose(sys.argv[1:])

# ======================================================================
