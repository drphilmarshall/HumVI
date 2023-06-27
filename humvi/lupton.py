# ======================================================================

"""
Functions for implementing the Wherry et al extension to the
Lupton et al algorithm for making color composite images.

Phil Marshall, Winter 2013

"""

# ======================================================================
# Globally useful modules:

import numpy
from PIL import Image

# ======================================================================

def lupton_intensity(r,g,b,type='sum'):

    if type == 'sum':
        return (r+g+b) + 1e-10

    elif type == 'rms':
        return numpy.sqrt(r*r+g*g+b*b) + 1e-10

# ----------------------------------------------------------------------

def lupton_stretch(I, Q, alpha):

    return numpy.arcsinh(alpha*Q*I) / (Q*I)

# ----------------------------------------------------------------------
# Clip high values to box:

def lupton_saturate(r, g, b, threshold):

    x = numpy.dstack((r, g, b))

    # Highest pixel-value at given position
    maxpix = numpy.max(x, axis=-1)
    maxpix[maxpix<1.0] = 1.0

    rr = r/maxpix
    gg = g/maxpix
    bb = b/maxpix

    return rr, gg, bb

# ======================================================================
# Testing:

if __name__ == '__main__':

    print("No tests defined")

# ======================================================================
