#!/usr/bin/env python
# =====================================================================

# Globally useful modules:

import numpy
import sys,getopt
from PIL import Image
import humvi
import os

# =====================================================================

def compose(rfile, gfile, bfile, scales=(1.0,1.0,1.0), Q=1.0, alpha=1.0, \
            masklevel=None, saturation='color', offset=0.0, backsub=False, \
            vb=False, outfile='color.png'):
    """
    Compose RGB color image.
    """

    # -------------------------------------------------------------------

    if vb:
        print "HumVI: Making color composite image of data in following files:",rfile,gfile,bfile
        print "HumVI: Output will be written to",outfile
        if masklevel is not None: print "HumVI: Masking stretched pixel values less than",masklevel

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

    rscale,gscale,bscale = humvi.normalize_scales(scales)
    red.set_scale(manually=rscale)
    green.set_scale(manually=gscale)
    blue.set_scale(manually=bscale)
    if vb: print 'HumVI: Scales normalized to:',red.scale,green.scale,blue.scale

    # Scale images - only do once:
    red.apply_scale()
    green.apply_scale()
    blue.apply_scale()

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Stretch images to cope with high dynamic range:

    if vb:
        print "HumVI: Stretch parameters Q,alpha:",Q,alpha
        print "HumVI: At low surface brightness levels, the channel images are further rescaled by alpha"
        print "HumVI: Nonlinearity sets in at about 1/Q*alpha in the scaled intensity image:",1.0/(Q*alpha)

    # Compute total intensity image and the arcsinh of it:
    I = humvi.lupton_intensity(red.image,green.image,blue.image,type='sum')
    stretch = humvi.lupton_stretch(I,Q,alpha)

    # Apply stretch to channel images:
    r = stretch * red.image
    g = stretch * green.image
    b = stretch * blue.image

    if masklevel is not None:
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

    if vb: print "HumVI: Image saved to:",outfile

    return

# ======================================================================
