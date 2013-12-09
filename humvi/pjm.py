# ======================================================================

"""
Functions for implementing Phil's tweaks to the Lupton algorithm for
making color composite images.
"""

# ======================================================================
# Globally useful modules:

import numpy,Image

# ======================================================================
# Add small offset to image, to make background look dark gray not black:

def pjm_offset(r,g,b,offset):

    rr = r + offset
    gg = g + offset
    bb = b + offset

    return rr,gg,bb

# ----------------------------------------------------------------------
# Detect problem areas in any of the channel images, and mask out:

def pjm_mask(r,g,b,threshold):

    tiny = 1e-10
    mask = r*0.0 + 1.0
    
    for image in (r,g,b):

        image[numpy.isnan(image)] = 0.0
        image[numpy.isinf(image)] = 0.0

        mask[image < threshold] = 0.0
        mask[(image > -tiny) & (image < tiny)] = 0.0

    return r*mask,g*mask,b*mask

# ======================================================================
