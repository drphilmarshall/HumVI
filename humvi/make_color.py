"""
Code to generate color composite images from data cubes.

This module contains code for generating color composite images using
the Lupton et al. algorithm.

Functions:
make_image -- the main function implementing the Lupton et al. algorithm.
write_color -- driver function for writing an image in a FITS file to a jpeg.
write_all -- driver function for making images of all cameras in a file.

The other functions aren't really meant for public consumption.

Note that this module requires, in addition to numpy and scipy, the
Python Imaging Library (http://www.pythonware.com/products/pil/)

Author: Patrik Jonsson
Shared freely at https://bitbucket.org/lutorm/python/src

"""

import Image, ImageColor, pyfits, os.path, types, pdb
from numpy import *
from scipy import ndimage


def extract_bands(imagedata, band):
    """
    Extracts the specified slices of an image cube.

    band can either contain a tuple of the slices or a comma-separated
    string of the slice numbers. The image cube is assumed to have the
    slice dimension last.

    """

    if isinstance(band, str):
        b = band.split(',',2)
        band =(int(b[0]),int(b[1]),int(b[2]))
    #pdb.set_trace()
    image = dstack((imagedata[:,:,band[0]],
            imagedata[:,:,band[1]],
            imagedata[:,:,band[2]]))
    return image

def make_energyunits(image, wavelengths=array((1,1,1)), **kwargs):
    """
    Convert the image data to energy units.

    The standard Sunrise images are in lambda units, ie luminosity per
    wavelength. To obtain a more wavelength-independent scale, energy
    or nu*F_nu units are preferable. This function converts a stack of
    3 images by multiplying by the wavelengths given.n

    """

    # convert image to nu*fnu
    return image*array(wavelengths)

def make_image(image,return_jpeg=True,**kwargs):
    """
    Converts a data array to a color image using the Lupton scheme.

    This function takes a 3D array as input and returns a color image,
    either as a jpeg data stream or as a 3D array that can be shown
    with imshow. The image array must be ordered (x,y,band).

    Keyword arguments:

    return_jpeg -- if True, the function returns the jpeg data for the
    image as a string that can be written to a file, otherwise it
    returns a x*y*3 RGB array with range 0-1 that can be fed to
    matplotlib.imshow.

    band -- the slices of the input data used for the r,g,b
    channels. Should be either a tuple or a comma-separated string.

    scale -- specifies the scale of the input data. If a tuple or
    comma-separated string, specifies the image values in the three
    colors where saturation happens. Can also be set to "auto" or
    "autolum", in which case the values will be determined
    automatically based on the "autopercentile" keyword (see
    below). With "auto", each color is scaled independently to its max
    value, with "autolum" the scale is the same for all 3 channels so
    the overall color is preserved.

    autopercentile -- when autoscaling, sets the percentile of the
    pixels that will be saturated, i.e. for a value of 0.1, the scale
    value will be such that 10% of the pixels will be above it. While
    in principle the scale should be set by the brightest pixel, in
    practice this often doesn't work well because a few pixels are
    often much brighter than the rest.

    The stretch is set by the parameters alpha (default .3) and Q
    (default 9). Alpha sets the low-intensity scaling and Q sets the
    saturation behavior.

    Usage notes:

    While the default values of Q and alpha are mostly good, finding a
    good scaling value can be tricky. For the normal Sunrise outputs,
    which are surface brightness in lambda units, the scaling not only
    depends on how bright the object is in the selected bands, but
    also varies drastically with the wavelength of the bands due to
    the wavelength dependence of the units. On the other hand, because
    the images are in surface brightness, they are reflective of the
    real surface brightnesses of galaxies, so should be approximately
    independent of the simulation details as long as the galaxies are
    reasonably realistic. For the default urz images
    that are commonly shown in the Sunrise papers, the scaling values
    used are (1.5, 1, 1.5).

    A typical example of using this function to show an image with
    matplotlib would be (im is the image array):

    imshow(make_color.make_image(im, return_jpeg=False, band=(6,4,2), \
         scale="autolum",autopercentile=0.1))

    """

    if kwargs.has_key("band"):
        band=kwargs["band"]
    else:
        band = (6,4,2)

    image = extract_bands(image, band)

    if kwargs.has_key("alpha"):
        alpha=float(kwargs["alpha"])
    else:
        alpha=.3

    if kwargs.has_key("Q"):
        Q=float(kwargs["Q"])
    else:
        Q=9

    if kwargs.has_key("autopercentile"):
        autopercentile=kwargs["autopercentile"]
    else:
        autopercentile=0.0005

    m=0

    sz = image.shape

    slices = sz[0]
    xs = sz [1]
    ys = sz [2]

    # deal with nan
    image = where(image != image, 0, image)

    if kwargs.has_key("scale"):
        if kwargs["scale"]=="auto":
            # set scale
            scale=find_autoscale(image, autopercentile, lum=False)
            print "autoscaling:"
            #raise RuntimeError, `scale `
        elif kwargs["scale"]=="autolum":
            scale=find_autoscale(image, autopercentile, lum=True)
            print "autoscaling:"
        elif isinstance(kwargs["scale"], str):
            sc = kwargs["scale"].split(',',2)
            scale=(float(sc[0]),float(sc[1]),float(sc[2]))
        else:
          scale=kwargs["scale"]
    else:
        scale=(1.,1.,1.)

    print "scale",scale

    image*=array(scale)

    r=image[:,:,0]
    g=image[:,:,1]
    b=image[:,:,2]

    i = (r+g+b)/3+1e-20
    r[:,:] = r*arcsinh (alpha*Q*(i-m))/(Q*i)
    g[:,:] = g*arcsinh (alpha*Q*(i-m))/(Q*i)
    b[:,:] = b*arcsinh (alpha*Q*(i-m))/(Q*i)

    image=clip(image,0.0,1.0)

    if not return_jpeg:
        # imshow wants 0-1 range numbers, so that's what we return
        return image

    # the jpeg encoding wants 0-255 uchars
    image = image*255;
    jpgimg = Image.fromarray(image.astype(uint8))
    jpg = jpgimg.tostring("jpeg","RGB", 95)

    return jpg

def write_color(file, hdu, output_file, overwrite=False,
                imagefunc=make_image, **kwargs):
    """
    Writes a color jpeg image by combining three slices of a FITS image.

    One often wants to just generate a jpeg image from an existing
    Sunrise output file. This function wraps make_image to make this
    easy.

    Inputs:
    file -- The name of the input FITS file.
    hdu -- The name of the FITS HDU in the input file that contains
    the image data.
    output_file -- The name of the output jpeg file.

    Keyword arguments:
    overwrite -- if True, a pre-existing output file is silently
    overwritten, otherwise such a case will be skipped with a warning.

    In addition, keyword arguments are passed through to make_color.

"""

    if not overwrite and os.path.exists(output_file):
        print "\tSkipping %s -- already present"%output_file
    else:
        f=open(output_file,'wb')
        f.write(imagefunc(transpose(pyfits.open(file)[hdu].data, axes=(1,2,0)), **kwargs))
        f.close()

def write_all(file, suffix, output_file, overwrite=False,
                imagefunc=make_image, **kwargs):
    """
    Writes color jpeg images for all cameras in a Sunrise a FITS file.

    One often wants to just generate jpeg images for all cameras in an
    existing Sunrise output file. This function wraps make_image to
    make this easy. The output file name is a python format string,
    where %d will be replaced with the camera number.

    Inputs:
    file -- The name of the input FITS file.
    suffix -- The suffix of the HDU to read. The HDU name is created
    by prepending CAMERA-<n> to this string.
    output_file -- The name of the output jpeg file as a python format
    string, where %d will be replaced with the camera number.

    Keyword arguments:
    overwrite -- if True, a pre-existing output file is silently
    overwritten, otherwise such a case will be skipped with a warning.

    In addition, keyword arguments are passed through to make_color.

    """

    print "Writing color images for file %s"%file
    inf=pyfits.open(file)
    ncam=inf['MCRX'].header['N_CAMERA']
    for i in range(ncam):
      hdu='CAMERA%d%s'%(i,suffix)
      ofnam=output_file%i
      if not overwrite and os.path.exists(ofnam):
        print "\tSkipping %s -- already present"%ofnam
        continue
      print "\t%s to %s"%(hdu,ofnam)
      of=open(ofnam,'wb')
      of.write(imagefunc(transpose(inf[hdu].data,(1,2,0)), **kwargs))
      of.close()
    inf.close()

identity = lambda x: x

def merge_all(files, suffixes, output_file, overwrite=False, \
          screenfile = None, screen=None, **kwargs):
    """Writes color jpgs for all camera HDUs with the specified suffix
    to the output file, where %d is replaced by the camera number. If
    screen is true, all images except the first will be attenuated
    according to the output from calculate_screen. The screen data is
    loaded from the first file."""
    print "Writing color images for file %s"%files

    if kwargs.has_key("band"):
      b = kwargs["band"].split(',',2)
      bands =(int(b[0]),int(b[1]),int(b[2]))
    else:
      bands = (2,1,0)
    # we now feed make_image the bands in RGB order
    kwargs['band']="0,1,2"

    inf=[pyfits.open(f) for f in files]
    if screen==None:
      screen=[False]*len(inf)
    ncam=inf[0]['MCRX'].header['N_CAMERA']
    for i in range(ncam):
      hdus=['CAMERA%d%s'%(i,suffix) for suffix in suffixes]
      ofnam=output_file%i
      if not overwrite and os.path.exists(ofnam):
        print "\tSkipping %s -- already present"%ofnam
        continue
      print "\t%s to %s"%(hdus,ofnam)
      #try:
      # crete data first so we don't create output file if input
      # file is truncated
      if screenfile != None:
        attn = calculate_screen(pyfits.open(screenfile), i, **kwargs)
        attenuate = lambda x: x*attn
      else:
        attenuate = identity

      attenuation_map = { True: attenuate, False: identity }
      attenuators = [ attenuation_map[x] for x in screen ]
      data=reduce(lambda x,y:x+y,
              [a(f[h].data[bands,:,:])
               for f,h,a in zip(inf,hdus, attenuators)])
      of=open(ofnam,'wb')
      of.write(make_image(data, **kwargs))
      of.close()
      #except AttributeError:
      #    print "\tError on %s -- skipping" % ofnam
    [f.close() for f in inf]


def calculate_screen(file, camera, smoothing=1, **kwargs):
    """Returns the attenuation factor images for derived from the
    surface density of metals in HDU CAMERAi-AUX. File is either a
    filename or a pyfits file object."""

    print "\tCalculating dust screen"

    if type(file)==types.StringType:
      file=pyfits.open(file)

    zslice = file['CAMERA%d-AUX'%camera].header['MM_SLICE']-1

    dust_crossection = 5.86e-6
    dust_fraction = 0.4
    # we want optical depth, which for Draine dust is
    # 5.86e-6*surface_dens at 550nm

    tauV = dust_crossection*dust_fraction* \
         file['CAMERA%d-AUX'%camera].data[zslice,:,:]
    tauV = ndimage.gaussian_filter(tauV, smoothing)
    #tauV[:,:]=10
    print "\tMax tau_v=%f"%tauV.max()
    (xs,ys) = tauV.shape

    # tauU/tauV = 1.6, tauI/tauV = 0.6

    attn = zeros((3,xs,ys))
    attn[0,:,:] = exp(-0.6*tauV)
    attn[1,:,:] = exp(-tauV)
    attn[2,:,:] = exp(-1.6*tauV)

    return attn

def make_hsv(value, hue, saturation=1, hue_smoothing=1, **kwargs):
    """Makes a color image based on two images of H and V in the HSV
    system."""

    if kwargs.has_key("alpha"):
      alpha=float(kwargs["alpha"])
    else:
      alpha=.3

    if kwargs.has_key("Q"):
      Q=float(kwargs["Q"])
    else:
      Q=9

    m=0

    (xs, ys) = value.shape

    # deal with nan
    # would be better to replace with median of surrounding pixels or something
    value = where(value != value, 0, value)
    hue = ndimage.gaussian_filter(hue, hue_smoothing)
    hue = where(hue != hue, 0, hue)

    if kwargs.has_key("vscale"):
      if kwargs["vscale"]=="auto":
        vscale=30./value.max()
        print "Autoscaling vscale=%e"%vscale
      else:
        vscale = float(kwargs["vscale"])
    else:
      vscale=1

    if kwargs.has_key("hscale"):
      if kwargs["hscale"]=="auto":
        hscale=30./hue.max()
        print "Autoscaling hscale=%e"%hscale
      else:
        hscale = float(kwargs["hscale"])
    else:
      hscale=1

    value = value*vscale + 1e-20
    value = value*arcsinh (alpha*Q*(value-m))/(Q*value) * 100
    value = clip(value,0,100)

    hue = hue*hscale + 1e-20
    hue = hue*arcsinh (alpha*Q*(hue-m))/(Q*hue) * 360
    hue = clip(hue,0,360)

    # now make the image by setting hsv values pixel by pixel
    img = Image.new("RGB",(xs,ys))
    pix = img.load()
    for x in range(xs):
      for y in range(ys):
        hsv = r"hsl(%d,%d%%,%d%%)"%(hue[x,y],saturation*100,value[x,y])
        # intentinal transpose here
        pix[y,x]=ImageColor.getrgb(hsv)

    jpg = img.tostring("jpeg","RGB")
    return jpg


def write_hsv(file,hdu,valueslice, hueslice,
          output_file, saturation=1, overwrite=False, **kwargs):
    if not overwrite and os.path.exists(output_file):
      print "\tSkipping %s -- already present"%output_file
    else:
      f=open(output_file,'wb')
      data=pyfits.open(file)[hdu].data
      f.write(make_hsv(data[valueslice,:,:],
                 data[hueslice,:,:],
                 saturation, **kwargs))
      f.close()

def write_all_hsv(file, suffix, valueslice, hueslice,
            output_file, saturation=1, overwrite=False, **kwargs):
    print "Writing hsv images for file %s"%file
    inf=pyfits.open(file)
    ncam=inf['MCRX'].header['N_CAMERA']
    for i in range(ncam):
      hdu='CAMERA%d%s'%(i,suffix)
      ofnam=output_file%i
      if not overwrite and os.path.exists(ofnam):
        print "\tSkipping %s -- already present"%ofnam
        continue
      print "\t%s to %s"%(hdu,ofnam)
      of=open(ofnam,'wb')
      of.write(make_hsv(inf[hdu].data[valueslice,:,:],
                  inf[hdu].data[hueslice,:,:],
                  saturation, **kwargs))
      of.close()
    inf.close()


def find_autoscale(image, autopercentile, lum=True):
    """Return the scale factor that saturates a certain percentage of
    pixels specified by autopercentile. If lum is true, a common
    scaling for the three channels is returned, otherwise each channel
    is scaled independently."""

    if lum:
      scim=array(image).ravel()
      scim.sort()
      autoind = -int(scim.shape[0]*autopercentile)
      scale = 1./scim[autoind]
      return (scale,scale,scale)
    else:
      # return three independent scalings
      scim = [array(sorted(array(image[:,:,i]).ravel())) for i in (0,1,2)]
      autoind = -(int(scim[0].shape[0]*autopercentile+1))
      #pdb.set_trace()
      scale = tuple([float(1./s[autoind]) for s in scim])
      return scale

def make_linear_image(image,return_jpeg=True,**kwargs):

    if kwargs.has_key("band"):
      band=kwargs["band"]
    else:
      band = (2,1,0)
    image = extract_bands(image, band)
    #pdb.set_trace()

    if kwargs.has_key("autopercentile"):
      autopercentile=kwargs["autopercentile"]
    else:
      autopercentile=0.0005

    # deal with nan
    image = where(image != image, 0, image)

    image = make_energyunits(image, **kwargs)

    if kwargs.has_key("scale"):
      if kwargs["scale"]=="auto":
        # set scale
        scale=find_autoscale(image, autopercentile, lum=False)
        print "autoscaling:"
        #raise RuntimeError, `scale `
      elif kwargs["scale"]=="autolum":
        scale=find_autoscale(image, autopercentile, lum=True)
        print "autoscaling:"
      elif isinstance(kwargs["scale"], str):
        sc = kwargs["scale"].split(',',2)
        scale=(float(sc[0]),float(sc[1]),float(sc[2]))
      else:
        scale=kwargs["scale"]
    else:
      scale=(1,1.,1)

    print "scale",scale

    # simple linear scale
    imgdata=clip(image*scale*255,0,255)

    if not return_jpeg:
      return imgdata

    jpgimg = Image.fromarray(imgdata.astype(uint8))
    jpg = jpgimg.tostring("jpeg","RGB", 95)

    return jpg


# NOT FINISHED. Need to load up images, extract slices, etc, and this
# should be done in generic functions for all the image making
# functions too.
def auto_aperture(frames, output_file, exposure_frames=20, psf=None,
            autopercentile=0.0005, imagefunc=make_linear_image,
            **kwargs):
    """Generates a set of images, simulating an autoexposure. Frames
    is a list of (file, hdu) tuples containing the movie frames, in
    order. Output_file is a format string for the output file name,
    which gets sent file and hdu."""
    scales=[]
    bands=kwargs["band"]
    kwargs["band"]=(0,1,2)
    wavelengths=kwargs["wavelengths"]
    kwargs["wavelengths"]=array((1,1,1))

    for f,h in frames:
      ofnam=output_file%(f,h)
      print "Writing color image for file %s, image %s to %s\n"%(f,h,ofnam)
      inf=pyfits.open(f)
      # transpose the slice dim before sending to extract_bands
      data=extract_bands(transpose(inf[h].data, axes=(1,2,0)),
                   band=bands)

      # make energy units
      data=make_energyunits(data, wavelengths=wavelengths)

      # get scale data and calculate boxcar smoothed scale
      scales.append(array(find_autoscale(data, autopercentile, lum=True)))
      n=len(scales)
      if n>exposure_frames:
        del scales[0]

      s=reduce(lambda x,y: x+y, scales, array((0,0,0)))/n

      of=open(ofnam,'wb')
      of.write(imagefunc(data, scale=s, **kwargs))
      of.close()
      inf.close()

def make_framelist(files, suffix):
    """Returns a list of tuples (file, hdu) containing all the images
    in the files, in order."""

    frames=[]
    for f in files:
      inf=pyfits.open(f)
      ncam=inf['MCRX'].header['N_CAMERA']
      print "File %s, %d cameras"%(f,ncam)
      curf=[(f, 'CAMERA%d%s'%(i,suffix)) for i in range(ncam)]
      frames+=curf
      inf.close()
    return frames

