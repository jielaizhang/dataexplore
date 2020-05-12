#!/usr/bin/env python

""" copy_fitsfile.py -- Copy the Xth extension of input fits file to new output file, including the Xth extension header.

Usage: convolve [-h] [-v] [-o SAVELOC] [--overwrite] <fitsfile> <fits_ext_number>

Arguments:
    fitsfile (string)
    	Path to image.
    fits_ext_number (integer)
        Extension number in fits file you want to copy.

Options:
    -h, --help                          Show this screen
    -v, --verbose                       Show extra information [default: False]     
    -o SAVELOC, --out SAVELOC           The output mask is saved here [default: ./convolved.fits]
    --overwrite                         Overwrite output file [default: False]

Examples:
"""

import docopt
import astropy.io.fits as fits
import sys, os

def copy_fitsext(saveloc, fitsfile, fits_ext_number,
                overwrite=False, 
                verbose=False):

    # Grab the right extension out of input fits file
    d,h = fits.getdata(fitsfile,ext=fits_ext_number,header=True)
              
    # saves as the fits file	
    if overwrite:
        fits.writeto(saveloc,d,h,overwrite=True)
    elif os.path.isfile(saveloc):
        print('******* ERROR *******')
        print(saveloc + ' already exists!')
        print('---------------------')
        sys.exit("sys.exit Error: Remove or rename "+saveloc+" before using this script again")
    else:
        fits.writeto(saveloc,d,h)

    print("File is saved in "+saveloc+". \n")

####################### BODY OF PROGRAM STARTS HERE ########################

if __name__ == "__main__":

    # Read in input arguments
    arguments       = docopt.docopt(__doc__)
    verbose         = arguments['--verbose']
    if verbose:
        print(arguments)
    saveloc         = arguments['--out']
    fitsfile        = arguments['<fitsfile>']
    fits_ext_number = int(arguments['<fits_ext_number>'])
    overwrite       = arguments['--overwrite']

    copy_fitsext(saveloc, fitsfile, fits_ext_number,
                overwrite=overwrite, 
                verbose=verbose)


