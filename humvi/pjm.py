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
    mask[r < threshold] = 0.0
    mask[(r > -tiny) & (r < tiny)] = 0.0
    mask[g < threshold] = 0.0
    mask[(g > -tiny) & (g < tiny)] = 0.0
    mask[b < threshold] = 0.0
    mask[(b > -tiny) & (b < tiny)] = 0.0

    return r*mask,g*mask,b*mask

# ======================================================================
