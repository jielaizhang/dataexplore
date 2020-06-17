#!/usr/bin/env python

""" determine_imageCentre.py -- Calculate the x,y pixel coordinates of the image center for a 2D fits image. 
Usage: determine_imageCentre.py [-h] [-v] [--debug] <fitsfile>

Arguments:
    fitsfile (string)

Options:
    -h, --help                              Show this screen
    -v, --verbose                           Show extra information, i.e. print center to screen [default: False]     
    --debug                                 Output more for debugging [default: False]

Examples:
"""

import docopt
import astropy.io.fits as fits
import sys

__author__      = "Jielai Zhang"
__license__     = "MIT"
__version__     = "1.0.1"
__date__        = "2020-06-017"

def print_verbose_string(printme,verbose=False,underscores=False):
    if verbose:
        if underscores:
            print("-" * len(f"VERBOSE: {printme}"),file=sys.stdout)
        else:
            print(f"VERBOSE: {printme}",file=sys.stdout)

def determine_imageCentre(fitsfile,verbose=False):
    '''Determine the pixel coordinates of the image center.'''

    # read in fits header
    h = fits.getheader(fitsfile)
    X = h['NAXIS1']
    Y = h['NAXIS2']
    xcenter = X/2.
    ycenter = Y/2.

    if verbose:
        print(f'The x,y centers in pixel coords are: {xcenter},{ycenter}')

    return xcenter, ycenter

####################### BODY OF PROGRAM STARTS HERE ########################

if __name__ == "__main__":

    # Read in input arguments
    arguments       = docopt.docopt(__doc__)
    verbose         = arguments['--verbose']
    debugmode       = arguments['--debug']
    if debugmode:
        print(arguments)   
    fitsfile        = arguments['<fitsfile>']

    _,_=determine_imageCentre(fitsfile,verbose=verbose)
