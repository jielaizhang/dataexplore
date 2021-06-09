#!/usr/bin/env python

"""nanify_with_mask_DECam.py -- Input an image, a matching _sd file, and nanify all areas of image that are not zeros. Outputs input.replace('.fits','_nanified.fits')
This script was designed for CP trigger DECam outputs. 

Usage: nanify_with_mask.py [-h] [-v] [-o DIRECTORY] <fitsfile> <sdfitsfile>

Options:
    -h, --help                              Show this screen
    -o DIRECTORY, --out DIRECTORY           Output file into this directory [default: ./]
    -v, --verbose                           Show extra information [default: False]   

Examples:
"""

import docopt
import astropy.io.fits as fits
import numpy as np

__author__      = "Jielai Zhang"
__license__     = "MIT"
__version__     = "1.0.1"
__date__        = "2021-06-09"
__maintainer__  = "Jielai Zhang"
__email__       = "zhang.jielai@gmail.com"

##############################################################
####################### Main Function ########################
##############################################################

def nanify_with_mask_DECam(fitsfile,sdfitsfile,savedir='./',verbose=False):

    # Get data
    d_sd = fits.getdata(sdfitsfile)
    d,h  = fits.getdata(fitsfile,header=True)

    # Make pixel gaps and saturated pixels nans
    change_to_nans = np.where(d_sd>0)
    d[change_to_nans]=np.nan

    # Save new file
    f_nanified_fits = fitsfile.replace('.fits','_withnans.fits')
    fits.writeto(f_nanified_fits,d,h,overwrite=True)
    if verbose:
        print(f'Saved: {f_nanified_fits}')
    
    return f_nanified_fits

############################################################################
####################### BODY OF PROGRAM STARTS HERE ########################
############################################################################

if __name__ == "__main__":

    # Read in input arguments
    arguments       = docopt.docopt(__doc__)
    verbose         = arguments['--verbose'] 
    savedir         = arguments['--out']
    fitsfile        = arguments['<fitsfile>']
    sdfitsfile      = arguments['<sdfitsfile>']

    _ = nanify_with_mask_DECam(fitsfile,sdfitsfile,savedir=savedir,verbose=verbose)
