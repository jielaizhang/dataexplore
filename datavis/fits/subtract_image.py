#!/usr/bin/python

#!/usr/bin/env python

""" subtract_image.py -- Input two visible wavelength astronomical fits files, image subtraction is carried out using hotpants. Assume images are aligned.
Usage: subtract_image.py [-h] [-v] [--debug] [-s SAVELOC] [--badpixmapsave] [-o] [--sextractor LOC] <fitsfile1> <fitsfile2>

Arguments:
    fitsfile1 (string)
    fitsfile2 (string) output is fitsfile2-fitsfile1

Options:
    -h, --help                              Show this screen
    -v, --verbose                           Show extra information [default: False]     
    --debug                                 Output more for debugging [default: False]
    -s SAVELOC, --save SAVELOC              Save subtraction image as fits file at this path. 
    --badpixmapsave                         If input, badpix map will be saved at saveloc_badpixmap.fits. This only happens if there is a -s saveloc [default: False]
    -o, --overwrite                         Overwrite subtraction and badpixmap if it exists, and saveloc specified. 
    --sextractor LOC                        Indicate location of source extractor for FWHM and sky calculations. [default: /opt/local/bin/source-extractor]

Examples:
    bash: python subtract_image.py file1.fits file2.fits --overwrite --save file1_sub_file2.fits 
    python: from datavis.fits.subtract_image import subtract_image
            subtract_image(fitsfile1,fitsfile2,saveloc=False|'./sub.fits',badpixmapsave=False|True,verbose=False|True,debugmode=False|True)
"""
import docopt
import astropy.io.fits as fits
import sys, os
import numpy as np
import subprocess

# Jielai written modules
from datastats.calculate_skyStats2 import calculate_skyStats
from datastats.calculate_FWHM import calculate_FWHM

__author__      = "Jielai Zhang"
__license__     = "MIT"
__version__     = "1.0.1"
__date__        = "2020-06-18"
__maintainer__  = "Jielai Zhang"
__email__       = "zhang.jielai@gmail.com"

#########################################
# ======= House Keeping Functions =======
#########################################
'''These functions standardise verbose, debug printing'''
def print_verbose_string(printme,verbose=False,underscores=False):
    if verbose:
        if underscores:
            print("#" * len(f"VERBOSE: {printme}"),file=sys.stdout)
            print(f"VERBOSE: {printme}",file=sys.stdout)
            print("#" * len(f"VERBOSE: {printme}"),file=sys.stdout)
        else:
            print(f"VERBOSE: {printme}",file=sys.stdout)

def print_debug_string(printme,debugmode=False,underscores=False):
    if debugmode:
        if underscores:
            print("#" * len(f"DEBUG  : {printme}"),file=sys.stdout)
            print(f"DEBUG  : {printme}",file=sys.stdout)
            print("#" * len(f"DEBUG  : {printme}"),file=sys.stdout)
        else:
            print(f"DEBUG  : {printme}",file=sys.stdout)

'''These functions help organise or delete files'''
def clearit(fname):
    if os.path.isfile(fname):
        os.remove(fname)
    return None

def ifexistexit(fname):
    # use exists instead of isfile because if it's a directory, also exit.
    if os.path.exists(fname):
        exitmsg = f'ERROR: file exists, delete/rename/choose different output before trying subtract_image again: {fname}'
        sys.exit(exitmsg)
    return None

##############################################################
####################### Main Function ########################
##############################################################

def subtract_image(fitsfile1,fitsfile2,saveloc=False,badpixmapsave=False,overwrite=False,sextractorloc='/opt/local/bin/source-extractor',verbose=False,debugmode=False):

    print_debug_string(f'sextractor loc specified in subtract_image is: {sextractorloc}',debugmode=debugmode)

    # Specify output image, system exit if image already exists (do not overwrite)
    # If saveloc provided, also set location of badpixmap output 
    if saveloc:
        saveit            = True
        saveloc_badpixmap = saveloc.replace('.fits','_badpixmap.fits')
        if not overwrite:
            ifexistexit(saveloc)
            ifexistexit(saveloc_badpixmap)
    else:
        saveit             = False
        saveloc            = './sub_temp.fits'
        saveloc_badpixmap  = './sub_badpixmap_temp.fits'
        print_verbose_string(f'Temporary hotpants output image will be saved at: {saveloc}.',verbose=verbose)
        print_verbose_string(f'Temporary hotpants output image will be saved at: {saveloc_badpixmap}.',verbose=verbose)

    # Set setting for printing out from subfunctions:
    if verbose:
        quietmode=False
    else:
        quietmode=True

    # Get parameters required by hotpants: FWHM_larger
    print_verbose_string('### Calculating FWHM of two input images ###',verbose=verbose,underscores=True)
    [FWHM1,FWHM2]   = calculate_FWHM([fitsfile1,fitsfile2],sextractorloc=sextractorloc,quietmode=quietmode,verbose=verbose)
    FWHM_larger     = max([FWHM1,FWHM2])
    if debugmode:
        print(f'DEBUG: FWHM1 & 2: {FWHM1},{FWHM2}')

    # Get parameters required by hotpants: sky values of two images
    print_verbose_string('### Calculating sky params of two input images ###',verbose=verbose,underscores=True)
    d                                     = fits.getdata(fitsfile1)
    naxis1,naxis2                         = np.shape(d)
    box_size                              = int(min(naxis1,naxis2)/6.)
    print_debug_string(f'Box size for calculating sky: {box_size}',debugmode=debugmode)
    sky1, sky_error1, sky_pixel_std1 = calculate_skyStats(fitsfile1, place_boxes = box_size, n_iterations = 10, sextractorloc=sextractorloc, verbose=verbose, debugmode=debugmode, quietmode=quietmode)
    sky2, sky_error2, sky_pixel_std2 = calculate_skyStats(fitsfile2, place_boxes = box_size, n_iterations = 10, sextractorloc=sextractorloc, verbose=verbose, debugmode=debugmode, quietmode=quietmode)    

    # Get parameters required by hotpants: actual input parameters calculated via FWHM_larger, and sky values of input images
    tl  = sky1-5*sky_pixel_std1
    il  = sky2-5*sky_pixel_std2
    r   = 2.5*FWHM_larger
    rss = 2*2.5*FWHM_larger

    if debugmode:
        print(f'DEBUG: sky1: {sky1}')
        print(f'DEBUG: sky2: {sky2}')
        print(f'DEBUG: tl: {tl}')
        print(f'DEBUG: il: {il}')
        print(f'DEBUG: r: {r}')
        print(f'DEBUG: rss: {rss}')

    
    # Run hotpants to do image subtraction
    print_verbose_string('### Running Hotpants to do Image Subtraction ###',verbose=verbose,underscores=True)
    if verbose:
        hotpants_command = f"hotpants -tmplim {fitsfile1} -inim {fitsfile2} -outim {saveloc} -il {il} -tl {tl} -r {r} -rss {rss} -n t -omi {saveloc_badpixmap}"
    else:
        hotpants_command = f"hotpants -tmplim {fitsfile1} -inim {fitsfile2} -outim {saveloc} -il {il} -tl {tl} -r {r} -rss {rss} -n t -omi {saveloc_badpixmap} -v 0"
    
    if debugmode:
        print(f'DEBUG: {hotpants_command}')
    subprocess.call(hotpants_command,shell=True)

    # Read in hotpants output files for returning
    subtraction_image = fits.getdata(saveloc)
    bad_pix_map       = fits.getdata(saveloc_badpixmap)

    # If don't need to save it, clear hotpnats output files
    print_verbose_string('### Image Subtraction Completed ###',verbose=verbose,underscores=True)
    if not saveit:
        clearit(saveloc)
        clearit(saveloc_badpixmap)
        print_verbose_string(f'Temporary hotpants output image Removed: {saveloc}.',verbose=verbose)
        print_verbose_string(f'Temporary hotpants output image Removed: {saveloc_badpixmap}.',verbose=verbose)
        print('No output path specified, so no subtraction image will be saved.')
    elif saveit and not badpixmapsave:
        clearit(saveloc_badpixmap)
        print_verbose_string(f'Temporary hotpants output image Removed: {saveloc_badpixmap}.',verbose=verbose)
        print(f'SAVED  : {saveloc}')
    else:
        print(f'SAVED  : {saveloc}')
        print(f'SAVED  : {saveloc_badpixmap}')

    return subtraction_image, bad_pix_map

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
 
    fitsfile1       = arguments['<fitsfile1>']
    fitsfile2       = arguments['<fitsfile2>']

    saveloc         = arguments['--save']
    badpixmapsave   = arguments['--badpixmapsave']
    overwrite       = arguments['--overwrite']
    sextractorloc   = arguments['--sextractor']

    _ = subtract_image(fitsfile1,fitsfile2,
                             saveloc=saveloc,badpixmapsave=badpixmapsave,overwrite=overwrite,
                             sextractorloc=sextractorloc,verbose=verbose,debugmode=debugmode)
