#!/usr/bin/env python
# ======================================================================

# Globally useful modules:

import numpy,sys,getopt,Image

import colorjpg

# ======================================================================

def fits2colorjpg(argv):
  """
  NAME
    fits2colorjpg.py

  PURPOSE
    Make a color composite JPG image from 3 FITS images, using the Lupton
    algorithm.

  COMMENTS
    Reads filter name and zero point from header and tries to set scales 
    automatically - unless scales are set on command line. Note 
    telescope/survey has to be recognised for this to work...

  USAGE
    fits2colorjpg.py [flags] red.fits green.fits blue.fits

  FLAGS
    -h            Print this message
    -v            Verbose operation
    --png         Make PNG instead of JPG 

  INPUTS
    red.fits etc  Names of FITS files containing image data

  OPTIONAL INPUTS
    -s --scales   3,2,4    Comma-separated scales for R,G,B channels [None]
    -n --nonlinearity  1.5 Non-linearity parameter
    -o --output   outfile  Name of output filename [guessed]

  OUTPUTS
    stdout        Useful information
    outfile       Output plot in jpg or png format


  EXAMPLES

    fits2colorjpg.py  CFHTLS_27_i_sci.fits CFHTLS_27_r_sci.fits CFHTLS_27_g_sci.fits
                      
  BUGS

  HISTORY
    2012-05-11 started Marshall (Oxford)
  """

  # --------------------------------------------------------------------

  try:
      opts, args = getopt.getopt(argv, "hv:s:p:n:o:",["help","verbose","png","scales","pars","output"])
  except getopt.GetoptError, err:
      # print help information and exit:
      print str(err) # will print something like "option -a not recognized"
      print fits2colorjpg.__doc__
      return

  vb = False
  
  png = False
  outfile = 'output.jpg'
  
  pars = '8,0.02'
  scales = 'Auto'
  
  for o,a in opts:
      if o in ("-h", "--help"):
          print fits2colorjpg.__doc__
          return
      elif o in ("-v", "--verbose"):
          vb = True
      elif o in ("-p","--parameters"):
          pars = a
      elif o in ("-s","--scales"):
          scales = a
      elif o in ("--png"):
          png = True
      elif o in ("-o","--output"):
          outfile = a
      else:
          assert False, "unhandled option"
  
  # Check for datafiles in array args:

  if len(args) == 3:
    rfile = args[0]
    gfile = args[1]
    bfile = args[2]
    if vb:
      print "Making color composite image of data in following files:",rfile,gfile,bfile
      if png: 
        "Output will be in PNG format"
      else:
        "Output will be in JPG format"
  else:
    print fits2colorjpg.__doc__
    return
    
  # Parse nonlinearity parameters:  
  Qs,alphas = pars.split(',')
  Q = float(Qs)
  alpha = float(alphas)
  
  if scales != 'Auto':
    x,y,z = scales.split(',')
    rscale = float(x)
    gscale = float(y)
    bscale = float(z)
    
  # --------------------------------------------------------------------
  # Read in images, set and apply scales etc:
  
  red = colorjpg.channel(rfile)
  green = colorjpg.channel(gfile)
  blue = colorjpg.channel(bfile)
  
  # Check shapes are equal:
  colorjpg.check_image_shapes(red.image,green.image,blue.image)
  
  # Optionally, subtract backgrounds, taken to be 3-sigma clipped median:
  
  # red.subtract_background()
  
  # Set scales:
  if scales == 'Auto':
    red.set_scale()
    green.set_scale()
    blue.set_scale()
  else:
    red.set_scale(manually=rscale)
    green.set_scale(manually=gscale)
    blue.set_scale(manually=bscale)
    
  if vb: 
    print 'Scales set to:',red.scale,green.scale,blue.scale
  
  # Scale images - only do once:
  red.apply_scale()
  green.apply_scale()
  blue.apply_scale()
  
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Lupton arcsinh scaling:

  if vb: print "Nonlinearity parameters Q,alpha:",Q,alpha

  I = colorjpg.lupton_intensity(red.image,green.image,blue.image)
  
  r = colorjpg.lupton_stretch(red.image,I,Q,alpha)
  g = colorjpg.lupton_stretch(green.image,I,Q,alpha)
  b = colorjpg.lupton_stretch(blue.image,I,Q,alpha)

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Package into a python Image, and write out to file:

  image = colorjpg.pack_up(r,g,b)
 
  image.save(outfile)
  if vb: print "Image saved to:",outfile

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  return

# ======================================================================

if __name__ == '__main__':
  fits2colorjpg(sys.argv[1:])

# ======================================================================
