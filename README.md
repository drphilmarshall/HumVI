HumVI
=====

Human Viewable (color composite) Image creation.

Purpose
-------
Creates a composite colour image from sets of input FITS files, following the Magain, Courbin & Sohy (1998) deconvolution algorithm and the Lupton et al (2004) composition algorithm (with extensions by Hogg & Wherry.)

Authors 
-------
Phil Marshall (Oxford) <dr.phil.marshall@gmail.com>
Cato Sandford (NYU)
Anupreeta More (IPMU)

Date
----
October, 2012

Usage
-----
There are two top level scripts.
deconvolve.py takes in some number of FITS files (or a directory containing FITS files) and deconvolves each image to a Gaussian PSF.
compose.py takes in 3 FITS files as input, and returns
a color composite, color-saturated png image with an arcsinh stretch. 

Example
-------
I have an image in three bandpasses, called i.fits, r.fits and g.fits, in the current directory. I wish to deconvolve each such that the PSF is a circularly-symmetrical Gaussian with \sigma=2 pixels.

	python deconvolve.py -p 2 i.fits r.fits g.fits

Now I wish to combine the three new images.

	python compose.py gDeconvolved_i.fits gDeconvolved_r.fits gDeconvolved_g.fits

View the result:

	display output.png

Dependencies
------------

The composition script requires:
* the Python Image Library (PIL) available from http://www.pythonware.com/products/pil/

The dconvolution script requires:

Both require:
* numpy
