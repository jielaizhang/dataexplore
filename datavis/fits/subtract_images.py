#!/usr/bin/env python

""" subtract_images.py -- Input a template and many images to subtract from template. Images should be visible wavelength astronomical fits files, image subtraction is carried out using hotpants. Assume input images are already aligned.

Not implemented yet: Inserted into the header will be TEMPLATE and SUBTRACT listing the paths to the two template and file that is subtracted. 

To be improved: if badpixmapsave = True, it will be saved in same location as sub image with filename subimagename_badpixmap.fits, however, it is currently not naturally passed on from subtract_images.py right now, it is just assumed that it is done this way. Should improve such that the file names are passed on from called functions. 

Usage: subtract_images.py [-h] [-v] [--debug] [-s SAVEDIR] [--badpixmapsave] [-o] [--sextractor LOC] <template> <fitsfiles>...

Arguments:
    template (string)
    fitsfiles (string) output is template-fitsfile

Options:
    -h, --help                              Show this screen
    -v, --verbose                           Show extra information [default: False]     
    --debug                                 Output more for debugging [default: False]
    -s SAVEDIR, --save SAVEDIR              Save subtraction image as fits file at this directory, with name input_sub.fits. [default: ./] 
    --badpixmapsave                         If true, save badpixmap at input_sub_badpixmap.fits [default: False]
    -o, --overwrite                         Overwrite subtraction and badpixmap if it exists, and saveloc specified. 
    --sextractor LOC                        Indicate location of source extractor for FWHM and sky calculations. [default: /opt/local/bin/source-extractor]

Examples:
    bash: python subtract_images.py --overwrite --save '/dir/save/here' template.fits file1.fits file2.fits
    bash: python subtract_images.py template.fits file1.fits file2.fits --overwrite --save '/dir/save/here' 
    python: from datavis.fits.subtract_images import subtract_images
            subtract_images(template.fits,[fits1,fits2,fits3],saveloc='./'|False,badpixmapsave=False|True, overwrite=True|False,sextractorloc='path/to/source-extractor,verbose=False|True,debugmode=False|True)
"""


import docopt
import astropy.io.fits as fits
import numpy as np
import sys, os, ntpath
from pathlib import Path

# Jielai Zhang's modules
from datavis.fits.subtract_image import subtract_image

__author__      = "Jielai Zhang"
__license__     = "MIT"
__version__     = "1.0.1"
__date__        = "2020-06-22"
__maintainer__  = "Jielai Zhang"
__email__       = "zhang.jielai@gmail.com"

#########################################
# ======= House Keeping Functions =======
#########################################
'''These functions standardise verbose, debug printing'''
def print_verbose_string(printme,verbose=False,underscores=False):
    if verbose:
        if underscores:
            print("@" * len(f"VERBOSE: {printme}"),file=sys.stdout)
            print(f"VERBOSE: {printme}",file=sys.stdout)
            print("@" * len(f"VERBOSE: {printme}"),file=sys.stdout)
        else:
            print(f"VERBOSE: {printme}",file=sys.stdout)

def print_debug_string(printme,debugmode=False,underscores=False):
    if debugmode:
        if underscores:
            print("@" * len(f"DEBUG  : {printme}"),file=sys.stdout)
            print(f"DEBUG  : {printme}",file=sys.stdout)
            print("@" * len(f"DEBUG  : {printme}"),file=sys.stdout)
        else:
            print(f"DEBUG  : {printme}",file=sys.stdout)

'''These functions help organise or delete files'''
def clearit(fname):
    if os.path.isfile(fname):
        os.remove(fname)
    return None

##############################################################
####################### Main Function ########################
##############################################################

def subtract_images(template,fitsfiles,savedir='./',badpixmapsave=False,overwrite=False,sextractorloc='/opt/local/bin/source-extractor',verbose=False,debugmode=False):
    
    # Start array to save subtracted images
    subtracted_images           = []
    badpixmap_images            = []
    paths_savedimages           = []
    paths_savedbadpixmapimages  = []

    # For each fitspath in fitsfiles, do template - fitspath 
    for fitspath in fitsfiles:

        print_verbose_string(f'Operating on {fitspath}.',underscores=True,verbose=verbose)

        # Determine saveloc for subtracted image based on input savedir and input file name
        fname       = ntpath.basename(fitspath)
        fname_stub  = Path(fname).stem
        saveloc     = savedir + os.path.sep + fname_stub + '_sub.fits'

        # Assume path to saved badpixmap if badpixmapsave is True.
        if badpixmapsave:
            badpixsaveloc = saveloc.replace('.fits','_badpixmap.fits')

        # Do image subtraction
        subtracted_image, badpixmap = subtract_image(template,fitspath,saveloc=saveloc,overwrite=overwrite,sextractorloc=sextractorloc,verbose=verbose,debugmode=debugmode)

        # Append to array of subtracted images
        subtracted_images.append(subtracted_image)
        badpixmap_images.append(badpixmap)

        # Append to array of subtracted image and badpixmap saved file locations.
        paths_savedimages.append(saveloc)
        if badpixmapsave:
            paths_savedbadpixmapimages.append(badpixsaveloc)

    # After all is done, print out all saved files again:
    print('\nSummarise what is saved:')
    if badpixmapsave:
        for f,b in zip(paths_savedimages,badpixmapsave):
            print(f'SAVED  :{f},{b}')
    else:
        for f in paths_savedimages:
            print(f'SAVED  :{f}')

    return subtracted_images

############################################################################
####################### BODY OF PROGRAM STARTS HERE ########################
############################################################################

if __name__ == "__main__":

    # Read in input arguments
    arguments       = docopt.docopt(__doc__)
    verbose         = arguments['--verbose']
    debugmode       = arguments['--debug']
    if debugmode:
        print(arguments)  
 
    template        = arguments['<template>']
    fitsfiles       = arguments['<fitsfiles>']

    savedir         = arguments['--save']
    badpixmapsave   = arguments['--badpixmapsave']
    overwrite       = arguments['--overwrite']
    sextractorloc   = arguments['--sextractor']

    _ = subtract_images(template,fitsfiles,savedir=savedir,badpixmapsave=badpixmapsave,overwrite=overwrite,sextractorloc=sextractorloc,verbose=verbose,debugmode=debugmode)

