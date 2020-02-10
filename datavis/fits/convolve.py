#!/usr/bin/env python

""" convolve.py -- Convolve sourceimage to a lower resolution image. Outputs the lower resolution image as a fits file.

Usage: convolve [-h] [-v] [-o SAVELOC] [--overwrite] (pixel | arcsec) <fitsfile> <init_res> <final_rez>

Arguments:
    fitsfile (string)
    	Path to image to be convolved.
    init_res (float)
    	Resolution of the current image (FWHM, in arcsec).
    final_rez (float)
    	Resolution of the desired image (FWHM, in arcsec).

Options:
    -h, --help                          Show this screen
    -v, --verbose                       Show extra information [default: False]     
    -o SAVELOC, --out SAVELOC           The output mask is saved here [default: ./convolved.fits]
    --overwrite                         Overwrite output file [default: False]

Examples:
"""

import docopt
import astropy.io.fits as fits
import numpy as np
import sys, os
import scipy.signal as signal
import astropy.wcs as wcs

def convolve(saveloc, fitsfile, initres, finalres,
            inputres_pixel=False, inputres_arcsec=True,
            overwrite=False,
            verbose=False):
    # Check initial resolution is lower than final
    if initres > finalres:
        sys.exit('sys.exit Error: init_res must be a smaller number than final_res')

    # Load in source data and header information
    data,header = fits.getdata(fitsfile,header = True)
    # nan to zero
    data[np.isnan(data)] = 0
    # number of x and y pixels
    x = header['NAXIS2']
    y = header['NAXIS1']
    
    highres   = initres
    lowres   = finalres

    if inputres_arcsec:
        if inputres_pixel:
            sys.exit('Please input either pixel or arcsec as the input resolutions, not both.')
        # get the pixelsize in arcsec, where the value is in degree*3600(arcsec/degree)
        #pixelsize = header['CDELT2']*3600
        w          = wcs.WCS(header)
        pixelsizes = wcs.utils.proj_plane_pixel_scales(w)*60*60
        pixelsize  = pixelsizes[0]
        # FWHM of current image
        FWHM_highres_pix = highres/pixelsize
        # FWHM of desired resolution 
        FWHM_lowres_pix = lowres/pixelsize
        # FWHM calulated for the gaussian of convolution kernel
        FWHM_kernel_pix = np.sqrt(FWHM_lowres_pix**2 - FWHM_highres_pix**2)
    elif inputres_pixel:
        FWHM_highres_pix = highres
        FWHM_lowres_pix  = lowres
        # FWHM calulated for the gaussian of convolution kernel
        FWHM_kernel_pix = np.sqrt(FWHM_lowres_pix**2 - FWHM_highres_pix**2)
    else:
        sys.exit('Please indicate whether the input resolutions are in pixel or arcsec units')
    # gaussian consant to convert sigma to FWHM 2*sqrt(2ln(2))
    constant = 2*np.sqrt(2*np.log(2))
    # sigma for the gaussian of convolution kernel
    sigma = FWHM_kernel_pix/constant
    # making the 2-D image of the convolution kernel by making 2, 1-D gaussian
    # and normalized by the gaussian normalization factor
    gauss1 = signal.general_gaussian(x,1, sigma)/((sigma)*np.sqrt(2*np.pi))
    gauss2 = signal.general_gaussian(y,1, sigma)/((sigma)*np.sqrt(2*np.pi))
    # combine these two to create 2D image
    kernel = np.outer(gauss1, gauss2)
    # convolve the image using signal.fftconvolve premade function and kernel
    convolved = signal.fftconvolve(data, kernel, mode='same')

    # change the resolution if it's there
    header['FWHM'] = (lowres, '[arcsec] Resolution of the map')
              
    # saves as the fits file	

    if overwrite:
        fits.writeto(saveloc, convolved, header,overwrite=True)
    elif os.path.isfile(saveloc):
        print('******* ERROR *******')
        print(saveloc + ' already exists!')
        print('---------------------')
        sys.exit("sys.exit Error: Remove or rename "+saveloc+" before using this script again")
    else:
        fits.writeto(saveloc, convolved, header)

    print("Convolving done. File is saved in "+saveloc+". \n")

####################### BODY OF PROGRAM STARTS HERE ########################

if __name__ == "__main__":

    # Read in input arguments
    arguments       = docopt.docopt(__doc__)
    verbose         = arguments['--verbose']
    if verbose:
        print(arguments)
    saveloc         = arguments['--out']
    fitsfile        = arguments['<fitsfile>']
    initres         = float(arguments['<init_res>']) 
    finalres        = float(arguments['<final_rez>']) 
    inputres_pixel     = arguments['pixel']
    inputres_arcsec    = arguments['arcsec']
    overwrite          = arguments['--overwrite']

    convolve(saveloc, fitsfile, initres, finalres, 
            inputres_pixel=inputres_pixel, inputres_arcsec=inputres_arcsec,
            overwrite=overwrite, 
            verbose=verbose)


