HumVI
=====

Human Viewable (color composite) Image creation.

Purpose
-------
Creates a composite colour image from sets of input FITS files, following the Lupton et al (2004) composition algorithm (with extensions by Hogg & Wherry.)

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
compose.py takes in 3 FITS files as input, and returns
a color composite, color-saturated png image with an arcsinh stretch. 

In the attic there is an attempt (deconvolve.py) at a reworked version of the 
Magain, Courbin & Sohy (1998) deconvolution algorithm, that is non-operational. The problem of how to bring
images from 3 different filters to a common resolution remains open. For now, don't go in the attic! 

Example
-------
I have an image in each of three bandpasses, called i.fits, r.fits and g.fits, in the current directory. 
I wish to combine the three new images into an RGB color composite.

	compose.py i.fits r.fits g.fits

Dependencies
------------

The composition script requires:
* the Python Image Library (PIL) available from http://www.pythonware.com/products/pil/
* numpy
