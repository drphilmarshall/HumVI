#!/usr/bin/env python
# =====================================================================

# Globally useful modules:

import numpy
import sys,getopt,Image
import humvi
import os

# =====================================================================

def compose(argv):
    """
    NAME
        compose.py

    PURPOSE
      Make a color composite PNG image from 3 FITS images, using the Lupton
      algorithm. Part of the HumVI package.

    COMMENTS
      Reads filter name and zero point from header and tries to set scales
      automatically - unless scales are set on command line. Note
      telescope/survey has to be recognised for this to work...

    USAGE
      compose.py [flags] red.fits green.fits blue.fits

    FLAGS
      -h        Print this message
      -v        Verbose operation
      -b --subtract-background

    INPUTS
      red.fits etc  Names of FITS files containing image data

    OPTIONAL INPUTS
      -s --scales   3,2,4     Comma-separated scales for R,G,B channels [None]
      -p --parameters  5,0.05 Non-linearity parameters Q,alpha
      -o --output   outfile   Name of output filename [guessed]
      -x --saturate-to style  Saturation method, color or white [white]
      -z --offset 0.1       Offset image level (+ve = gray background)
      -m --mask -1.0        Mask images below level [don't]

    OUTPUTS
      stdout      Useful information
      outfile     Output plot in png format


    EXAMPLES
    
        Note the different alpha required, to cope with survey depth!
    
        CFHTLS:
          compose.py  -v -s 0.4,0.6,1.7 -z 0.0 -p 1.7,0.09 -m -1.0 \
          -o examples/CFHTLS_27_gri.png \
          examples/CFHTLS_27_i_sci.fits \
          examples/CFHTLS_27_r_sci.fits \
          examples/CFHTLS_27_g_sci.fits

        PS1 (needs optimizing on a larger image):
          compose.py  -v -s 0.4,0.6,1.7 -z 0.0 -p 1.7,300.0 -m -1.0 \
          -o examples/H1413+117_10x10arcsec_riz.png \
          examples/H1413+117_10x10arcsec_55377.34051_z_sci.fits \
          examples/H1413+117_10x10arcsec_55665.51546_i_sci.fits \
          examples/H1413+117_10x10arcsec_55664.39704_r_sci.fits

        DES (very experimental - vary alpha first!):
          compose.py -v -s 2.0,1.5,2.5 -z 0.0 -p 2.0,0.04 -m -1.0 \
          -o bullet_gri.cropped.png \
          bullet_i.10.cropped.fits \       
          bullet_r.12.cropped.fits \
          bullet_g.11.cropped.fits

    BUGS

    - Renormalized scales will not be appropriate if one image is in very
      different units to another, or if images are in counts, not counts per
      second or AB maggies.

    - Code currently assumes one file per channel, whereas we might want to
      use N>3 images in combination.  The image to channel transformation 
      would be performed after scaling to flux units but before
      stretching. Masking would need to be done at this point too.

    - Should cope with just two, or even one, input image file(s).

    HISTORY
      2012-05-11    started Marshall (Oxford)
      2012-07-??    integrated with wherry.py Sandford (NYU)
      2012-12-19    defaults set for CFHTLS images Marshall (Oxford)
      2013-02-21    defaults set for PS1 images Marshall (Oxford)
      2013-02-26    experimenting with DES images Hirsch (UCL)
    """

    # -------------------------------------------------------------------

    try:
        opts, args = getopt.getopt(argv, "hvp:s:n:o:x:z:blwm:",\
        ["help","verbose","scales","pars","output","saturate-to","offset","subtract-background","lupton","wherry","mask"])
    except getopt.GetoptError, err:
        # print help information and exit:
        print str(err) # will print something like "option -a not recognized"
        print compose.__doc__
        return

    vb = False

    outfile = "color.png"

    # Defaults optimized for CFHTLS...
    pars = '1.7,0.09'
    scales = '0.4,0.6,1.7'
    
    # More general sensible choices:
    backsub = False
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
        else:
            assert False, "Unhandled option"


    # Check for datafiles in array args:

    print len(args)

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

    # Parse nonlinearity parameters:
    Qs,alphas = pars.split(',')
    Q = float(Qs)
    alpha = float(alphas)
    
    # Parse channel colour scales:
    x,y,z = scales.split(',')
    rscale,gscale,bscale = float(x),float(y),float(z)

    # -------------------------------------------------------------------
    # Read in images, calibrated into flux units:

    band3 = humvi.channel(rfile)
    band2 = humvi.channel(gfile)
    band1 = humvi.channel(bfile)

    # Check shapes are equal:
    humvi.check_image_shapes(band1.image,band2.image,band3.image)

    # Subtract backgrounds (median, optional):
    if backsub:
      band1.subtract_background()
      band2.subtract_background()
      band3.subtract_background()

    # -------------------------------------------------------------------
    
    # BUG: as it stands, this code assumes one file one channel, whereas
    # in practice we might like to be able to make composites based on 
    # N bands. Need to overload + operator for channels? Calib etc will
    # need altering as well as image.

    red = band3
    green = band2
    blue = band1

    # -------------------------------------------------------------------
    # Set scales determining color balance in composite:

    rscale,gscale,bscale = humvi.normalize_scales(rscale,gscale,bscale)
    red.set_scale(manually=rscale)
    green.set_scale(manually=gscale)
    blue.set_scale(manually=bscale)
    if vb: print 'Scales normalized to:',red.scale,green.scale,blue.scale

    # Scale images - only do once:
    red.apply_scale()
    green.apply_scale()
    blue.apply_scale()

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Stretch images to cope with high dynamic range:

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

        # BUG: this should have been done after scaling but before conversion
        # to channels, as its the individual images that have problems...

        r,g,b = humvi.pjm_mask(r,g,b,masklevel)

    # Offset the stretched images to make zero level appear dark gray.
    # Negative offset makes background more black...
    r,g,b = humvi.pjm_offset(r,g,b,offset)

    if saturation == 'color':
        # Saturate to colour at some level - might as well be 1, since
        # Q redefines scale?:
        threshold = 1.0
        r,g,b = humvi.lupton_saturate(r,g,b,threshold)
    # Otherwise, saturate to white.

    # Package into a python Image, and write out to file:
    image = humvi.pack_up(r,g,b)
    image.save(outfile)


# ======================================================================

    if vb: print "Image saved to:",outfile
    return

# ======================================================================

if __name__ == '__main__':
    compose(sys.argv[1:])

# ======================================================================
