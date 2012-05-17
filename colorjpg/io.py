# ======================================================================

''' 
Functions for reading and writing sets of fits files, including header
information.
'''

# ======================================================================
# Globally useful modules:

import numpy,pyfits,Image

vb = 0

# ======================================================================

class channel:

   def __init__(self,fitsfile):
   
      # Read in image and header:
      hdulist = pyfits.open(fitsfile)
      self.image = hdulist[0].data      
      self.hdr = hdulist[0].header
      hdulist.close()

      return
      
   def set_scale(self,manually=False):
   
      if manually:
        self.scale = manually
        
      else:  
        # Get zero point and exptime:
        self.zpt = extract_zeropoint(self.hdr)
        self.exptime = self.hdr['EXPTIME']

        # Extract filter and infer wavelength, telescope:
        self.filter = self.hdr['filter']
        self.wavelength = filter2wavelength(self.filter)

        # Suggest scale for this image based on exptime and zpt:
        # self.scale = numpy.sqrt(self.exptime)*10**(self.zpt + 2.5*numpy.log10(self.wavelength) - 26.0)
        self.scale = 10**(self.zpt + 2.5*numpy.log10(self.wavelength) - 26.0)
        # This will be some crazy large number, in general.
        
      return
   
   def apply_scale(self):
   
      self.image *= self.scale
      
      return
      
# ======================================================================

def normalize_scales(s1,s2,s3):
   mean = (s1 + s2 + s3)/3.0
   return s1/mean, s2/mean, s3/mean

# ----------------------------------------------------------------------

def extract_zeropoint(hdr):

   # Where did the data come from?
   origin = hdr['ORIGIN']

   # Look for zero points:
   if origin == 'CFHT':
     zpt = hdr['PHOT_C']
   
   else:
     raise "Unrecognised origin: "+origin
   
   return zpt
   
# ----------------------------------------------------------------------
   
def filter2wavelength(fname):

   # CFHT MegaCam (from http://www.cfht.hawaii.edu/Instruments/Imaging/Megacam/specsinformation.html)
   if fname == 'u.MP9301':
     L = 3740
   elif fname == 'g.MP9401':
     L = 4870
   elif fname == 'r.MP9601':
     L = 6250
   elif fname == 'i.MP9701':
     L = 7700
   elif fname == 'z.MP9801':
     L = 9000

   # SDSS:
   
   # DES:
   
   # etc
   return L

# ----------------------------------------------------------------------
   
def check_image_shapes(r,g,b):

   if (numpy.shape(r) != numpy.shape(g)) or \
      (numpy.shape(r) != numpy.shape(b)):
      raise "Image arrays are of different shapes, exiting"

   return
  
# ----------------------------------------------------------------------
# Make an 8 bit integer image cube from three channels:
   
def pack_up(r,g,b):

   NX,NY = numpy.shape(r)

   x = numpy.zeros([NX,NY,3])
   x[:,:,0] = numpy.flipud(r)
   x[:,:,1] = numpy.flipud(g)
   x[:,:,2] = numpy.flipud(b)

   x = numpy.clip(x,0.0,1.0)
   x = x*255

   return Image.fromarray(x.astype(numpy.uint8))

# ======================================================================


