HumVI
=====

Human Viewable (color composite) Image creation.

Purpose
-------
Creates a composite colour image from sets of input FITS files, following the Lupton et al (2004) composition algorithm (with extensions by Hogg & Wherry.)

Authors 
-------
* Phil Marshall (KIPAC) <dr.phil.marshall@gmail.com>
* Cato Sandford (NYU)
* Anupreeta More (IPMU)
* Hugo Buddelmeijer (Kapteyn)

Usage
-----
The main executable script is `compose.py`. It takes in 3 FITS files as input, and returns
a color composite, color-saturated png image with an arcsinh stretch. Make sure this script is on your PATH, and than the `humvi` directory is on your PYTHONPATH.

**Example:** I have an image in each of three bandpasses, called i.fits, r.fits and g.fits, in the current directory. 
I wish to combine the three new images into an RGB color composite.

	compose.py  -s 0.4,1.0,1.7  -p 1.0,0.02  -o gri.png  i.fits r.fits g.fits

Run `compose.py -h` to read about these options, but basically you can choose the 
contrast via the parameters `-p Q,alpha` and then the R,G,B color balance with the 
scales `-s R,G,B`. Good strategy is first to set `alpha` with `Q=1` so you can just see 
the noise, then adjust Q if necessary to brighten up the features, and finally choose 
the scales so that the noise looks like an equal mixture of red, green and blue. This 
latter step may not be possible if the different channel images have very different 
sensitivity.

As well as Lupton's arcsinh stretch, a major advantage of the HumVI algorithm is that 
you can set the same scales and parameters for all your images, so that you can compare 
them across your survey. HumVI reads the zero points out of the FITS headers and uses 
them to put all the images on the same flux scale (ie so that all the pixel values have 
the same units), so do make sure your images are well calibrated and have informative, 
accurate headers.

**Notes:** 
In the attic there is an attempt (deconvolve.py) at a reworked version of the 
Magain, Courbin & Sohy (1998) deconvolution algorithm, that is non-operational. The problem of how to bring
images from 3 different filters to a common resolution remains open. For now, don't go in the attic!

Dependencies
------------

The composition script requires:
* the Python Image Library (PIL) available from http://www.pythonware.com/products/pil/
* numpy
