#!/usr/bin/env python

"""copy_wcs.py -- read in wcs of 1 image and put it in header of another. Useful to copy wcs from original to simulated image, especially when simulated image contains extra not in the real sky sources. 

Usage: copy_wcs [-h] [-v] <copyfrom> <copyto>

Options:
    -h, --help                                  Show this screen
    -v, --verbose                               Show extra information [default: False]      

Examples:
    python copy_wcs.py withwcs.fits withoutwcs.fits 
"""

import docopt
import astropy.io.fits as fits

####################### BODY OF PROGRAM STARTS HERE ########################

if __name__ == "__main__":

    arguments = docopt.docopt(__doc__)

    # Mandatory argument
    copyfrom = arguments['<copyfrom>']
    copyto   = arguments['<copyto>']

    # Non-mandatory options without arguments
    verbose     = arguments['--verbose']

    # Read in wcs header
    h = fits.getheader(copyfrom)

    # Read in image and add wcs, save
    f           = fits.open(copyto)
    newf        = fits.PrimaryHDU()
    newf.header = f[0].header
    newf.data   = f[0].data
    newf.header['CTYPE1']   = h['CTYPE1']
    newf.header['CRPIX1']   = h['CRPIX1']
    newf.header['CRVAL1']   = h['CRVAL1']
    newf.header['CTYPE2']   = h['CTYPE2']
    newf.header['CRPIX2']   = h['CRPIX2']
    newf.header['CRVAL2']   = h['CRVAL2']
    newf.header['CD1_1']    = h['CD1_1']
    newf.header['CD1_2']    = h['CD1_2']
    newf.header['CD2_1']    = h['CD2_1']
    newf.header['CD2_2']    = h['CD2_2']
    #newf.header['RADECSYS'] = h['RADECSYS']
    newf.header['EQUINOX']  = h['EQUINOX']
    saveloc     = copyto.split('.')[0]+'_wcs.fits'
    newf.writeto(saveloc, clobber=True) 

