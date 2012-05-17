colorjpg
========

Color composite image creation from sets of input FITS files following the
Lupton et al (2004) algorithm.

> Author: Phil Marshall <dr.phil.marshall@gmail.com> 
>
> Date: 16 May, 2012

Based on: Patrick Jonsson's make_color.py, the jpg creation routines in
IDLutils, and most of all Lupton et al 2004.

Usage:

The top level script fits2colorjpg.py takes 3 fits files as input, and returns
a color composite jpg image. The scales can be set manually, or estimated
automatically to reflect the nu*f_nu of the objects' SEDs. 
