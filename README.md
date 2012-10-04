HumVI
=====

Human Viewable (color composite) Image creation, from sets of input FITS files following the
Lupton et al (2004) algorithm with extensions by Hogg & Wherry.

> Authors: Cato Sandford and Phil Marshall <dr.phil.marshall@gmail.com>
>
> Date: October, 2012

Based on the NW jpg creation routines in IDLutils, and Lupton et al 2004.

Usage:

The top level script compose.py takes in 3 fits files as input, and returns
a color composite, color-saturated png image with an arcsinh stretch. 
