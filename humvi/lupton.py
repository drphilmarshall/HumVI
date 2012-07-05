# ======================================================================

''' 
Functions for implementing the Lupton algorithm for making color 
composite images.
'''

# ======================================================================
# Globally useful modules:

import numpy,Image

# ======================================================================

def lupton_intensity(r,g,b):

  return (r+g+b)/3.0 + 1e-20

# ----------------------------------------------------------------------

def lupton_stretch(image,mean,Q,alpha):
 
  return image*numpy.arcsinh(alpha*Q*mean) / (Q*mean)

# ======================================================================
# Testing:

if __name__ == '__main__':
    
  print "No tests defined"

# ======================================================================
