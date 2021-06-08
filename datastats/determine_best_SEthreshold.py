#!/usr/bin/env python

""" determine_best_SEthreshold.py -- input a fits file, make a plot to help decide what SE detection threshold to use. Outputs a png.
Temporary outputs: ./basename(fitsfile)_temporarydirectoryXXXXX/basename(fitsfile).replace('.fits','_inv_temporary.fits').
Temporary outputs: ./basename(fitsfile)_temporarydirectoryXXXXX/ a lot of SE cat files.
                    Where XXX will be 5 random letters.
Temporary outputs: ./run_sourceextractor* SE required files
Temporary outputs will be deleted after the script finishes running.

!!! HARD CODED in, the FWHM in arc sec of the input image (1.0 for DECam)
!!! HARD CODED in, the min det area (10)
!!! HARD CODED in, the SE thresholds to try

Usage: determine_best_SEthreshold.py [-h] [-v] [--debug] [-q] [-o SAVELOC] [--SE_loc PATH] <fitsfile>

Options:
    -h, --help                              Show this screen
    -v, --verbose                           Show extra information [default: False]   
    -q, --quietmode                         Minimize print to screen. This is useful when this function is called in another function. [default: False]  
    --debug                                 Output more for debugging; Won't run through all thresholds, just 2 [default: False]
    -o SAVELOC, --out SAVELOC               Saved output as [default: ./SEthreshold.png]
    --SE_loc PATH                           Path to SE [default: /opt/local/bin/source-extractor]

Examples:
"""

__author__      = "Jielai Zhang"
__license__     = "MIT"
__version__     = "1.0.1"
__date__        = "2021-06-08"
__maintainer__  = "Jielai Zhang"
__email__       = "zhang.jielai@gmail.com"

import docopt
import shutil
import astropy.io.ascii as ascii
import astropy.io.fits as fits
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import ntpath
import os
from pathlib import Path
import string    
import random

font = {'family' : 'normal',
        'weight' : 'normal',
        'size'   : 20}
matplotlib.rc('font', **font)


# Jielai's dataexplore modules
from datastats.run_sourceextractor import run_sourceExtractor

##############################################################
####################### Main Function ########################
##############################################################
def determine_best_SEthreshold(fitsfile,saveloc='./SEthreshold.png',
                                SE_loc='/opt/local/bin/source-extractor',
                                verbose=False,debugmode=False,quietmode=False):

    # Set thresholds to try:
    # !!! HARD CODED IN
    thresholds = np.array([0.7,0.9,0.95,1.0,1.05,1.1,1.3,1.5,2.1,2.5,3.5,4.5,5.5,6.5])
    if debugmode:
        thresholds=np.array([1.0,6.5])

    # Determine save directory
    savedir = os.sep.join(saveloc.split(os.sep)[0:-1]) + os.sep
    # If savedir does not exist, create
    if not os.path.isdir(savedir):
        os.makedirs(savedir)

    # Create temporary directory
    fname       = ntpath.basename(fitsfile)
    fname_stub  = Path(fname).stem
    S=5
    ran = ''.join(random.choices(string.ascii_uppercase + string.digits, k = S)) 
    temp_dir    = '.'+os.sep+fname_stub+'_temporarydirectory'+ran
    if not os.path.isdir(temp_dir):
        os.makedirs(temp_dir)
        if debugmode:
            print(f'DEBUG: created: {temp_dir}')
    
    # Create temporary negative of input fits image
    d = fits.getdata(fitsfile)
    new_d = d*-1
    f_inv = temp_dir+os.sep+fname.replace('.fits','_inv_temporary.fits')
    fits.writeto(f_inv,new_d,overwrite=True)
    if debugmode:
        print(f'DEBUG: saved: {f_inv}')

    # Create catalogues for input image
    if verbose:
        print(f'DEBUG: running SE {len(thresholds)} times for diff det thresh on input image...')
    for DT in thresholds:
        cat, catted_f = run_sourceExtractor([fitsfile],spreadmodel=False,
                                savecats_dir = None, catending = None,
                                sextractorloc=SE_loc,
                                fwhm = 1.0, detect_minarea = 10, detect_thresh = DT, 
                                verbose=False,quietmode=True,debugmode=False)
        cat_basename       = ntpath.basename(cat[0])
        cat_basename_new   = cat_basename.replace('.cat',f'_{DT}.cat')
        cat_new            = temp_dir+os.sep+cat_basename_new
        shutil.move(cat[0],cat_new)

    # Create catalogues for negative of input image
    if verbose:
        print(f'DEBUG: running SE {len(thresholds)} times for diff det thresh on input image*-1...')
    for DT in thresholds:
        cat, catted_f = run_sourceExtractor([f_inv],spreadmodel=False,
                                savecats_dir = None, catending = None,
                                sextractorloc=SE_loc,
                                fwhm = 1.0, detect_minarea = 10, detect_thresh = DT, 
                                verbose=False,quietmode=True,debugmode=False)
        cat_basename       = ntpath.basename(cat[0])
        cat_basename_new   = cat_basename.replace('.cat',f'_{DT}.cat')
        cat_new            = temp_dir+os.sep+cat_basename_new
        shutil.move(cat[0],cat_new)


    # Get how many sources are extracted for each SE run
    ns         = []
    ns_inv     = []
    for DT in thresholds:
        cat_basename       = ntpath.basename(fitsfile)
        cat_basename_new   = cat_basename.replace('.fits',f'_{DT}.cat')
        cat_new            = temp_dir+os.sep+cat_basename_new
        f_open=open(cat_new,'rU')
        lines = f_open.readlines()
        n = len(lines)
        f_open.close()
        ns.append(n)
        cat_basename       = ntpath.basename(f_inv)
        cat_basename_new   = cat_basename.replace('.fits',f'_{DT}.cat')
        cat_new            = temp_dir+os.sep+cat_basename_new
        f_open=open(cat_new,'rU')
        lines = f_open.readlines()
        n = len(lines)
        f_open.close()
        ns_inv.append(n)

    # plot

    index = np.where(thresholds==1.0)[0][0]

    fig, axs = plt.subplots(3,figsize=(18,15))
    axs[0].plot(thresholds,ns,'x')
    axs[0].set_xlabel('SE detection threshold')
    axs[0].set_ylabel('#SEsrcs, image')
    axs[0].axvline(thresholds[index],label=f'det thresh ={thresholds[index]}')
    axs[0].legend()

    axs[1].plot(thresholds,ns_inv,'x')
    axs[1].set_xlabel('SE detection threshold')
    axs[1].set_ylabel('#SEsrcs,NEG image')
    axs[1].axvline(thresholds[index],label=f'det thresh ={thresholds[index]}')
    axs[1].legend()

    axs[2].plot(ns,ns_inv,'x')
    axs[2].set_xlabel('#SEsrcs in image')
    axs[2].set_ylabel('##SEsrcs,NEG image')
    axs[2].axvline(ns[index],label=f'det thresh ={thresholds[index]}')
    axs[2].legend()

    plt.savefig(saveloc)
    if not quietmode:
        print(f'Saved: {saveloc}')

    # Clean up and remove everything in temporary folder and the temporary folder itself
    shutil.rmtree(temp_dir) 

    return None


############################################################################
####################### BODY OF PROGRAM STARTS HERE ########################
############################################################################

if __name__ == "__main__":

    # Read in input arguments
    arguments       = docopt.docopt(__doc__)
    verbose         = arguments['--verbose']
    debugmode       = arguments['--debug']
    quietmode       = arguments['--quietmode']
    if debugmode:
        print(arguments)   
    saveloc         = arguments['--out']
    SE_loc          = arguments['--SE_loc']
    fitsfile        = arguments['<fitsfile>']

    _ = determine_best_SEthreshold(fitsfile,saveloc=saveloc,
                                    SE_loc=SE_loc,
                                    verbose=verbose,debugmode=debugmode,quietmode=quietmode)
