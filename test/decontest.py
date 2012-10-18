#!/usr/bin/env/python
##============================================================

"""
Some tests for deconvolution

They aren't extreme
"""


import sys
sys.path.append("../decon/")
sys.path.append("../misc/")
import DeconvolveToTargetPSF as DT
import numpy,scipy

def decontest():
	
	## Original image parameters
	imwidth = 50
	o_pars = [imwidth, 5., 5., imwidth/2, imwidth/2, 1.]
	## Noise model paramters
	mu = 0.01; sigma = 0.005
	## Reference PSF paramters	
	rPSF_pars = [imwidth, 2., 2., imwidth/2, imwidth/2, 1.]
	## Deconvolution
	krn_size = 15
	
	##----------------------------------------------------------
	
	## Original image, a point at (0,0) with a Gaussian PSF
	origarr = DT.Gauss_2D(*o_pars)
	scipy.misc.imsave("original.png",origarr)
	## Noise
	origarr += 0.3*origarr.max()*scipy.random.standard_normal(origarr.shape)+0.2*origarr.max()
	scipy.misc.imsave("original+noise.png",origarr)
	
	## Sensible range
	origarr.clip(0,1)
	
	##----------------------------------------------------------
	
	## Reference PSF
	psf_ref = DT.Gauss_2D(*rPSF_pars)
	scipy.misc.imsave("decon_psf_ref.png",psf_ref)

	## Calculate deconvolution-kernel array
	kernarr = DT.get_kernel(origarr, psf_ref,[krn_size,krn_size], False)
	## Deconvolve image
	decarr = DT.deconvolve_image(origarr, kernarr)
	
	##----------------------------------------------------------
	
	## Residuals
	psf_ref = DT.shave(decarr.shape, psf_ref)
	resarr = psf_ref - decarr
	resarr -= resarr.min()	## Shift zero
	resarr /= resarr.max()	## Put in range [0,1]
	### Need to stretch ###
	scipy.misc.imsave("decon_residuals.png", resarr)
	
	return

	
	
##============================================================
##============================================================

## Extract kernel for convolved->clean Gaussian.

## WORKS

def kerntest():
		
	origarr = DT.Gauss_2D(*[50, 5., 5., 25, 25, 1.])
	nextarr = scipy.signal.fftconvolve(origarr,numpy.array([[0.25,0.5],[0.0,0.25]]), "same")
	kernarr = DT.get_kernel(origarr, nextarr, [15,15], False)
	
	scipy.misc.imsave("kern_kern.png",kernarr)
	return None

##-------------------------------------------------------------

## More complicated scenario

def kerntest2():
		
	origarr = numpy.zeros([100,100]); origarr[50,50]=1.0
	scipy.misc.imsave("kern2_orig.png",origarr)
	
	nextarr = DT.Gauss_2D(*[100, 3., 3., 50, 50, 1.])
	scipy.misc.imsave("kern2_next.png",nextarr)
	
	kernarr = DT.get_kernel(nextarr, origarr, [30,30], False)	## OR SWITCH
	scipy.misc.imsave("kern2_kern.png",kernarr)
	
	dec = scipy.signal.fftconvolve(nextarr, kernarr, "valid")
	scipy.misc.imsave("kern2_dec.png",dec)
	

	return None

##============================================================

## Deconvolve one thing
## Do and undo.

## WORKS
def dectest():	
	
	origarr = DT.Gauss_2D(*[50, 5., 5., 25, 25, 1.])
	scipy.misc.imsave("dec_orig.png",origarr)
	
	## Alternatives
	#kernel = numpy.ones([1,1])
	#kernel = numpy.array([[0,-1,0],[-1,4,-1],[0,-1,0]])
	kernel = numpy.power((numpy.arange(0,100).reshape([10,10])), 2)
	
	nextarr = scipy.signal.fftconvolve(origarr,kernel, "valid")
	scipy.misc.imsave("dec_next.png",nextarr)
	
	kernarr = DT.get_kernel(DT.shave(nextarr.shape, origarr), nextarr, [5,5], False)
	scipy.misc.imsave("dec_kern.png",kernarr)
	
	decoarr = scipy.signal.fftconvolve(origarr, kernarr, "valid")	
	scipy.misc.imsave("dec_odec.png",decoarr)
	
	return
	
##============================================================


def imdectest():
	
	## Make composite image
	origarr = DT.Gauss_2D(*[50, 7., 7., 25, 25, 1.])
	bigorig = numpy.tile(origarr, [2,2])
	scipy.misc.imsave("imdec_obig.png",bigorig)
	
	## Target
	psf_ref = DT.Gauss_2D(*[50, 2., 2., 25, 25, 1.])
	scipy.misc.imsave("imdec_psfref.png",psf_ref)
	
	## Kernel
	kernarr = DT.get_kernel(psf_ref,origarr, [5,5], False)###SWITCH
	print kernarr
	#kernarr = numpy.array([[0,-1,0],[-1,4,-1],[0,-1,0]])
	scipy.misc.imsave("imdec_kern.png",kernarr)
	
	## Deconvolution
	decarr = DT.deconvolve_image(psf_ref, kernarr)	###SWITCH
	#decarr = scipy.signal.fftconvolve(bigorig, kernarr,"same")
	scipy.misc.imsave("imdec_odec.png",decarr)
	
	return	
	
##============================================================
if __name__=="__main__":
	kerntest2()
