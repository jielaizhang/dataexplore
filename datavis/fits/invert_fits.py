#!/usr/bin/python

#!/usr/bin/env python

""" invert_fits.py -- Input list of fits files, specify output directory where inverted files are to be saved. Output files will be saved there for each file.fits as file_neg.fits.

Usage: invert_fits [-h] [-q] [-o OUTDIR] [--overwrite] <fitsfiles>...

Arguments:
    fitsfiles (fitsfiles)

Options:
    -h, --help                          Show this screen
    -q, --quietmode                     Minimize print to screen. This is useful when this function is called in another function. [default: False]  
    -o OUTDIR, --out OUTDIR             Saved output as [default: ./]
    --overwrite                         Overwrite any existing files at destination [default: False]

Examples:
"""
import docopt
import astropy.io.fits as fits
import sys, os, ntpath
from pathlib import Path

__author__      = "Jielai Zhang"
__license__     = "MIT"
__version__     = "1.0.1"
__date__        = "2020-11-07"
__maintainer__  = "Jielai Zhang"
__email__       = "zhang.jielai@gmail.com"

##############################################################
####################### Main Function ########################
##############################################################

def invert_fits(fitsfiles,out_dir='./',quietmode=False,overwrite=False):

    Inverted_Files = []

    for f in fitsfiles:

        # Determine saveloc for inverted image based on input savedir and input file name
        fname       = ntpath.basename(f)
        fname_stub  = Path(fname).stem
        saveloc     = out_dir + os.path.sep + fname_stub + '_neg.fits'

        # Read in fits file
        d,h = fits.getdata(f,header=True)

        # Invert 
        d_neg = -1*d

        # Saving fits file
        fits.writeto(saveloc,d_neg,header=h,overwrite=overwrite) 
        if not quietmode:
            printme = f'SAVED  : {saveloc}'
            print(printme)

        # Add inverted file to output list
        Inverted_Files.append(saveloc)

    return Inverted_Files

############################################################################
####################### BODY OF PROGRAM STARTS HERE ########################
############################################################################

if __name__ == "__main__":

    # Read in input arguments
    arguments       = docopt.docopt(__doc__)
    quietmode       = arguments['--quietmode']
    out_dir         = arguments['--out']
    fitsfiles       = arguments['<fitsfiles>']
    overwrite       = arguments['--overwrite']

    _ = invert_fits(fitsfiles,out_dir=out_dir,quietmode=quietmode,overwrite=overwrite)
