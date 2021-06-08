#!/usr/bin/env python

"""run_sourceextractor.py -- read in fits files, get catalog with spread model.
Currently uses -PIXEL_SCALE 0, which assumes there is WCS in input fits, so can get pixel scale from header wcs.

Usage: run_sourceextractor [-h] [-q] [--debug] [-v] [-s LOC] [-p LOC] [--spreadmodel] [--savecats LOC] [--catending STRING] [--fwhm FLOAT] [--detect_thresh FLOAT] [--detect_minarea INT] <fitsfiles>...

Options:
    -h, --help                  Show this screen
    -q, --quiet                 Quiet mode, suppress printout of parameters measured.
    --debug                     Show extra information and save extra files, including the .cat files [default: False]
    -v, --verbose               Show extra information [default: False]     
    -s LOC, --sextractor LOC    Location of source extractor [default: /opt/local/bin/source-extractor] 
    -p LOC, --psfex LOC         Location of PSFEx [default: /opt/local/bin/psfex]
    --spreadmodel               If on, spread_model, mag_model, magerr_model will be in output catalog [default: False]
    --savecats LOC              Output cat files to this directory. If not input, outputs to same directory as input file.
    --catending STRING          Output cats for each file.fits will be file_STRING.cat, if empty, output will be file.cat.
    --fwhm FLOAT                The SEEING_FWHM config param for Source Extraction in arcsec [default: 1.2]
    --detect_thresh FLOAT       The Source Extraction detection threshold [default: 1.5]
    --detect_minarea INT        The Source Extraction detection minimum area [default: 5]

Examples:
    python run_sourceextractor.py *fits -o ./outfile.txt --savecats ../cats
    from datastats.run_sourceextractor import run_sourceextractor
        cats = run_sourceextractor(fitsfiles,
                                    extractorloc='/opt/local/bin/source-extractor',
                                    psfexloc='/opt/local/bin/psfex',
                                    verbose=False,quietmode=False,
                                    debugmode=False)
"""

import docopt, os
import astropy.io.fits as fits
import subprocess
from astropy.io import ascii
import pandas as pd
import numpy as np
import ntpath, shutil
import time
import string    
import random
import shutil

# Jielai Modules
from datastats.measure_psf import measure_psf
from misc.write_sourceextractor_files import write_sourceextractor_files

__author__      = "Jielai Zhang"
__license__     = "MIT"
__version__     = "1.0.1"
__date__        = "2020-11-12"
__maintainer__  = "Jielai Zhang"
__email__       = "zhang.jielai@gmail.com"

####################### Auxiliary Functions #######################
def remove_temp_files(fs):
    for f in fs:
        os.remove(f)
    return None

def mkdir_ifnotexist(directory):
    if not os.path.isdir(directory):
        os.makedirs(directory)
    return None

####################### MAIN function #######################

def run_sourceExtractor(fitsfiles,spreadmodel=False,
                        savecats_dir = None, catending = None,
                        sextractorloc='/opt/local/bin/source-extractor',
                        psfexloc='/opt/local/bin/psfex',
                        fwhm = 1.2, detect_minarea = 5, detect_thresh = 1.5, # Default settings in Source Extractor
                        verbose=False,quietmode=False,debugmode=False):

    if savecats_dir == None:
        savecats_dir = './'

    # Create temporary folder for source extractor input files
    S=5
    ran = ''.join(random.choices(string.ascii_uppercase + string.digits, k = S))
    temp_dir = './SE_temporary_'+ran
    os.makedirs(temp_dir)

    # Create temporary files required for Source Extractor
    (nnw_path, 
    conv_path, 
    params_path, 
    config_path) = write_sourceextractor_files(temp_dir,'run_sourceextractor',spreadmodel=spreadmodel,
                                                sextractorloc=sextractorloc,quietmode=quietmode)


    # Print useful information if verbose
    if verbose:
        if spreadmodel:
            print('The Steps for 1 file at a time are:')
            print('1. Measure PSF.')
            print('1a. Run Source Extractor')
            print('1b. Run PSFEx')
            print('2. Run Source Extractor with spread_model as a parameter output')
        print('#####################################')
        print('###### RUNNING SOURCE EXTRACTOR #####')
        print('#####################################\n')

    # Set Source Extractor verbose type
    if verbose:
        VERBOSE_TYPE = 'NORMAL'
    else:
        VERBOSE_TYPE = 'QUIET'

    # Run Source Extractor one image at a time
    catfiles = []
    psffiles = []
    catted_fitsfiles = []
    PSF_success = [False]*len(fitsfiles)
    for ii, f in enumerate(fitsfiles):

        # Determine where cat file is to be saved (full path including file name)
        fname       = ntpath.basename(f)
        if catending:
            catalog_name = savecats_dir+os.path.sep+fname.replace('.fits','_'+catending+'.cat')
        else:
            catalog_name = savecats_dir+os.path.sep+fname.replace('.fits','.cat')

        # Print separator
        if verbose:
            print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
            print(f'Input {f}')
            print(f'Intended output {catalog_name}')
            print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

        # Run Source Extractor then PSFEx
        if spreadmodel:
            try:
                if verbose:
                    print('Currently measuring PSF using source extractor and psfex...')
                [f_psf] = measure_psf([f], outdir=savecats_dir, savepsffits=False,
                                        sextractorloc=sextractorloc,
                                        psfexloc=psfexloc,
                                        verbose=verbose,quietmode=True)
                PSF_success[ii] = True
            except:
                print(f'\nSKIPPED: PSF measurement unsuccessful for {f}')
                continue
        else:
            # If spread_model is false, then don't need to run it, so set as true to run below source extraction. 
            PSF_success[ii] = True

        if PSF_success[ii] == True:
            # Run Source Extractor on image to get final catalog output
            if spreadmodel:
                command =   f'{sextractorloc} -c {config_path} '\
                            f'-CATALOG_NAME {catalog_name} '\
                            f'-CATALOG_TYPE ASCII_HEAD '\
                            f'-PARAMETERS_NAME {params_path} -FILTER_NAME {conv_path} '\
                            f'-STARNNW_NAME {nnw_path} -PIXEL_SCALE 0 -MAG_ZEROPOINT 25.0 '\
                            f'-PSF_NAME {f_psf} -PSF_NMAX 1 -PATTERN_TYPE GAUSS-LAGUERRE '\
                            f'-VERBOSE_TYPE {VERBOSE_TYPE} '\
                            f'-SEEING_FWHM {fwhm} -DETECT_MINAREA {detect_minarea} -DETECT_THRESH {detect_thresh} '\
                            f'{f}'
            else:
                command =   f'{sextractorloc} -c {config_path} '\
                            f'-CATALOG_NAME {catalog_name} '\
                            f'-CATALOG_TYPE ASCII_HEAD '\
                            f'-PARAMETERS_NAME {params_path} -FILTER_NAME {conv_path} '\
                            f'-STARNNW_NAME {nnw_path} -PIXEL_SCALE 0 -MAG_ZEROPOINT 25.0 '\
                            f'-VERBOSE_TYPE {VERBOSE_TYPE} '\
                            f'-SEEING_FWHM {fwhm} -DETECT_MINAREA {detect_minarea} -DETECT_THRESH {detect_thresh} '\
                            f'{f}'
            if verbose:
                print('Currently running source extractor to output required catalog...')
                print('=============================================')
                print('Executing command: %s\n' % command)
            try:
                rval = subprocess.run(command.split(), check=True)
                catfiles.append(catalog_name)
                catted_fitsfiles.append(f)
                if spreadmodel:
                    os.remove(f_psf)
                if verbose:
                    print(f'Success! Cat saved: {catalog_name}')
            except subprocess.CalledProcessError as err:
                print('\nCould not run SExtractor with exit error %s\n'%err)
                print('Command used:\n%s\n'%command)

    # Remove temporary files required for Source Extractor
    shutil.rmtree(temp_dir)
    
    return catfiles, catted_fitsfiles 

############################################################################
####################### BODY OF PROGRAM STARTS HERE ########################
#######################              MAIN           ########################
############################################################################

if __name__ == "__main__":

    arguments = docopt.docopt(__doc__)

    # Mandatory argument
    fitsfiles = arguments['<fitsfiles>']

    # Non-mandatory options without arguments
    quietmode       = arguments['--quiet']
    verbose         = arguments['--verbose']
    debugmode       = arguments['--debug']
    savecats_dir    = arguments['--savecats']
    catending       = arguments['--catending']
    sextractorloc   = arguments['--sextractor']
    psfexloc        = arguments['--psfex']
    spreadmodel     = arguments['--spreadmodel']
    fwhm            = arguments['--fwhm']
    detect_thresh   = arguments['--detect_thresh'] 
    detect_minarea  = arguments['--detect_minarea']

    # Calculate
    _,_ = run_sourceExtractor(fitsfiles,spreadmodel=spreadmodel,
                                sextractorloc=sextractorloc,
                                psfexloc=psfexloc,
                                savecats_dir = savecats_dir, catending = catending,
                                fwhm=fwhm, detect_thresh=detect_thresh,detect_minarea=detect_minarea,
                                verbose=verbose,quietmode=quietmode,debugmode=debugmode)
