# ======================================================================

"""
Functions for reading and writing sets of fits files, including header
information.
"""

# ======================================================================
# Globally useful modules:

import numpy,pyfits,os
from PIL import Image

vb = 0

# ======================================================================

# BUG: This code implies one file one channel, whereas we want to make
# composites based on N images... Class should be image, not channel.
# RGB channels should be constructed *after* scaling but *before* stretching

class channel:

    def __init__(self,fitsfile):

        self.input = fitsfile
        # Read in image and header:
        hdulist = pyfits.open(self.input)
        # self.hdr = hdulist[0].header
        # self.image = hdulist[0].data
        # Picking -1 header assumes we have 1 extension or PS1 (2 ext, image is last)
        self.image = hdulist[-1].data
        self.hdr = hdulist[-1].header
        self.calibrate()
        hdulist.close()

        return

    def calibrate(self):
        
        # Which telescope took these data?
        self.get_origin()
        
        # Get zero point, exptime:
        self.get_zeropoint()
        # EXPTIME is 1.0 for images in counts/s - but headers do not always
        # get this right... get_exptime gets this right.
        self.get_exptime()
        # Airmass? Gain? Should be included in zeropoint.

        # print "Image statistics for "+self.input
        # print "  ",self.origin,self.exptime,self.zpt

        # # Report 5 sigma depth:
        # image = self.image.copy()
        # mean = numpy.average(image)
        # stdev = numpy.std(image)
        # nsigma = 1
        # clip = 3
        # # Apply clipping:
        # while nsigma > 0.01:
        #     index = numpy.where(abs((self.image - mean)/stdev) < clip)[0]
        #     image = image[index]
        #     newmean = numpy.average(image)
        #     newstdev = numpy.std(image)
        #     nsigma = abs(mean - newmean)/newstdev
        #     mean = newmean
        #     stdev = newstdev
        # print "  Before calibration, mean, rms = ",mean,stdev
        # depth = -2.5*numpy.log10(5.0*stdev) + self.zpt
        # print "  Approximate 5-sigma limiting magnitude: ",depth
        
        # Compute calibration factor for image pixel values to 
        # convert them into flux units. The 30 is arbitrary, and 
        # simply determines the absolute value of alpha required 
        # for a nice image. 
        self.calib = (10.0**(0.4*(30.0 - self.zpt))) / self.exptime
        self.image *= self.calib

        # # Report 5 sigma depth:
        # image = self.image.copy()
        # mean = numpy.average(image)
        # stdev = numpy.std(image)
        # nsigma = 1
        # clip = 3
        # # Apply clipping:
        # while nsigma > 0.01:
        #     index = numpy.where(abs((self.image - mean)/stdev) < clip)[0]
        #     image = image[index]
        #     newmean = numpy.average(image)
        #     newstdev = numpy.std(image)
        #     nsigma = abs(mean - newmean)/newstdev
        #     mean = newmean
        #     stdev = newstdev
        # print "  After calibration, mean, rms = ",mean,stdev
        
        return
        
    def get_origin(self):
        if 'TELESCOP' in self.hdr:  
            if self.hdr['TELESCOP'] == 'CFHT 3.6m':
                self.origin = 'CFHT'
            elif self.hdr['TELESCOP'] == 'ESO-VLT-U0':
                # Assume that all data from ESO-VLT-U0 is from KIDS.
                self.origin = "KIDS"
            else:
                self.origin = self.hdr['TELESCOP']
        elif 'ORIGIN' in self.hdr:  
            if self.hdr['ORIGIN'] == 'CFHT':
                self.origin = 'CFHT'
            elif self.hdr['ORIGIN'] == 'DES':
                self.origin = 'DES'
            else:
                self.origin = self.hdr['ORIGIN']
        elif 'PSCAMERA' in self.hdr:
            self.origin = 'PS1'
        elif 'FID_ZP' in self.hdr:
            self.origin = 'DES'
        elif 'PROV' in self.hdr:
            self.origin = 'VICS82'
        else:
            self.origin = 'UNKNOWN'
        return

    def get_zeropoint(self):
        if self.origin == 'CFHT':
            if 'MZP_AB' in self.hdr:
                self.zpt = self.hdr['MZP_AB']
            elif 'MAGZP' in self.hdr:
                self.zpt = self.hdr['MAGZP']
            # elif 'PHOT_C' in self.hdr:
            #     self.zpt = self.hdr['PHOT_C']
            else:
                self.zpt = 30.0    
        elif self.origin == 'PS1':
            self.zpt = self.hdr['HIERARCH FPA.ZP']
        elif self.origin == 'DES':
            self.zpt = self.hdr['MZP_AB']
        elif self.origin == 'VICS82':
            if 'MZP_AB' in self.hdr:
                self.zpt = self.hdr['MZP_AB']
            else:
                self.zpt = 30.0
        elif self.origin == 'KIDS':
            # KiDS coadds are calibrated to ZPT=0.
            self.zpt = 0.0
        else: # UNKNOWN
            self.zpt = 30.0
        return

    def get_exptime(self):
        # Here we assume that both CFHT and PS1 provide images with 
        # pixel values in counts per second... or that the zero point
        # takes into account the exptime.
        if self.origin == 'CFHT':
            self.exptime = 1.0    
        elif self.origin == 'PS1':
            # self.exptime = self.hdr['EXPTIME']
            self.exptime = 1.0    
        elif self.origin == 'DES':
            # self.exptime = self.hdr['EXPTIME']
            self.exptime = 1.0    
        elif self.origin == 'VICS82':
            # self.exptime = self.hdr['EXPTIME']
            self.exptime = 1.0    
        elif self.origin == 'KIDS':
            # self.exptime = self.hdr['EXPTIME']
            self.exptime = 1.0    
        else: #UNKNOWN
            # Use 1.0 as default to ensure the program doesn't crash.
            self.exptime = 1.0
        return

    def set_scale(self,manually=False):
        if manually:
             self.scale = manually
        else:
             self.scale = 1.0
        return

    def apply_scale(self):
        self.image *= self.scale
        return

    def subtract_background(self):
        self.image -= numpy.median(self.image)
        return
        
    def writefits(self):
        self.output = str.split(self.input,'.')[0]+'_calibrated.fits'
        if os.path.exists(self.output): os.remove(self.output) 
        hdu = pyfits.PrimaryHDU()
        hdu.header = self.hdr
        hdu.data = self.image
        hdu.verify()
        hdu.writeto(self.output)
        return

# ======================================================================

def normalize_scales(s1,s2,s3):
    mean = (s1 + s2 + s3)/3.0
    return s1/mean, s2/mean, s3/mean

# ----------------------------------------------------------------------

def filter2wavelength(fname):

    # CFHT MegaCam (from http://www.cfht.hawaii.edu/Instruments/Imaging/Megacam/specsinformation.html)
    if fname == 'u.MP9301':
        L = 3740
    elif fname == 'g.MP9401':
        L = 4870
    elif fname == 'r.MP9601':
        L = 6250
    elif fname == 'i.MP9701' or 'i.MP9702':
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


