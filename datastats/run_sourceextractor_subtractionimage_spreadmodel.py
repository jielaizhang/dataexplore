#!/usr/bin/env python

"""run_sourceextractor_subtractionimage_spreadmodel.py -- Create subtraction image catalog with spread_model. Read in subtraction image fits files, and corresponding science image fits file, get catalog with spread model. PSF used for spread model is obtained from science image and not subtraction image. 
Currently uses -PIXEL_SCALE 0, which assumes there is WCS in input fits, so can get pixel scale from header wcs.

If input list (with --list option) the two lists must be matched. E.g. 5th in sci corresponds to 5th in sub list. 

Usage: run_sourceextractor_subtractionimage_spreadmodel [-h] [-q] [--debug] [-v] [-s LOC] [-p LOC] [--savecats LOC] [--catending STRING] [--fwhm FLOAT] [--detect_thresh FLOAT] [--detect_minarea INT] [--list] <sub> <sci> 

Options:
    -h, --help                  Show this screen
    -q, --quiet                 Quiet mode, suppress printout of parameters measured.
    --debug                     Show extra information and save extra files, including the .cat files [default: False]
    -v, --verbose               Show extra information [default: False]     
    -s LOC, --sextractor LOC    Location of source extractor [default: /opt/local/bin/source-extractor] 
    -p LOC, --psfex LOC         Location of PSFEx [default: /opt/local/bin/psfex]
    --savecats LOC              Output cat files to this directory. If not input, outputs to same directory as input file.
    --catending STRING          Output cats for each file.fits will be file_STRING.cat, if empty, output will be file.cat.
    --fwhm FLOAT                The SEEING_FWHM config param for Source Extraction in arcsec [default: 1.2]
    --detect_thresh FLOAT       The Source Extraction detection threshold [default: 1.5]
    --detect_minarea INT        The Source Extraction detection minimum area [default: 5]
    --list                      If specified, the input sub, sci should be text time lists [default: False]

Examples:
    python run_sourceextractor_subtractionimage_spreadmodel.py sub.fits sci.fits -o ./outfile.txt --savecats ../cats
    from datastats.run_sourceextractor_subtractionimage_spreadmodel import run_sourceextractor_subtractionimage_spreadmodel
        cats = run_sourceextractor_subtractionimage_spreadmodel(sub_list, sci_list, list=True
                                                                extractorloc='/opt/local/bin/source-extractor',
                                                                psfexloc='/opt/local/bin/psfex',
                                                                verbose=False,quietmode=False,
                                                                debugmode=False)
        cat = run_sourceextractor_subtractionimage_spreadmodel(subfits, scifits, list=False
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
import astropy.io.ascii as ascii

# Jielai Modules
from datastats.measure_psf import measure_psf
from misc.write_sourceextractor_files import write_sourceextractor_files

__author__      = "Jielai Zhang"
__license__     = "MIT"
__version__     = "1.0.1"
__date__        = "2020-11-16"
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

def run_sourceextractor_subtractionimage_spreadmodel_(sub, sci,
                        savecats_dir = None, catending = None,
                        sextractorloc='/opt/local/bin/source-extractor',
                        psfexloc='/opt/local/bin/psfex',
                        fwhm = 1.2, detect_minarea = 5, detect_thresh = 1.5, # Default settings in Source Extractor
                        verbose=False,quietmode=False,debugmode=False):

    # Track if PSF was obtained.
    PSF_success = False

    # Determine where cat file is to be saved (full path including file name)
    fname       = ntpath.basename(sub)
    if catending:
        catalog_name = savecats_dir+os.path.sep+fname.replace('.fits','_'+catending+'.cat')
    else:
        catalog_name = savecats_dir+os.path.sep+fname.replace('.fits','.cat')

    # Print separator
    if verbose:
        print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
        print(f'%% Input {sub}')
        print(f'%% Input {sci}')
        print(f'%% Intended output {catalog_name}')
        print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

    # Run Source Extractor then PSFEx
    try:
        if verbose:
            print('Currently measuring PSF using source extractor and psfex...')
        [f_psf] = measure_psf([sci], outdir=savecats_dir, savepsffits=False,
                                sextractorloc=sextractorloc,
                                psfexloc=psfexloc,
                                verbose=verbose,quietmode=True)
        PSF_success = True
    except:
        print(f'\nSKIPPED: PSF measurement unsuccessful for {f}')

    if PSF_success == True:

        # Create temporary files required for Source Extractor
        (nnw_path, 
        conv_path, 
        params_path, 
        config_path) = write_sourceextractor_files('./','run_sourceextractor',spreadmodel=True,
                                                    sextractorloc=sextractorloc,quietmode=quietmode)

        # Set Source Extractor verbose type
        if verbose:
            VERBOSE_TYPE = 'NORMAL'
        else:
            VERBOSE_TYPE = 'QUIET'

        # Run Source Extractor on image to get final catalog output
        command =   f'{sextractorloc} -c {config_path} '\
                    f'-CATALOG_NAME {catalog_name} '\
                    f'-CATALOG_TYPE ASCII_HEAD '\
                    f'-PARAMETERS_NAME {params_path} -FILTER_NAME {conv_path} '\
                    f'-STARNNW_NAME {nnw_path} -PIXEL_SCALE 0 -MAG_ZEROPOINT 25.0 '\
                    f'-PSF_NAME {f_psf} -PSF_NMAX 1 -PATTERN_TYPE GAUSS-LAGUERRE '\
                    f'-VERBOSE_TYPE {VERBOSE_TYPE} '\
                    f'-SEEING_FWHM {fwhm} -DETECT_MINAREA {detect_minarea} -DETECT_THRESH {detect_thresh} '\
                    f'{sub}'
        if verbose:
            print('Currently running source extractor to output required catalog...')
            print('=============================================')
            print('Executing command: %s\n' % command)
        try:
            rval = subprocess.run(command.split(), check=True)
            os.remove(f_psf)
            if verbose:
                print(f'Success! Cat saved: {catalog_name}')
        except subprocess.CalledProcessError as err:
            print('\nCould not run SExtractor with exit error %s\n'%err)
            print('Command used:\n%s\n'%command)

        # Remove temporary files required for Source Extractor
        remove_temp_files([nnw_path,conv_path,params_path,config_path])


    return catalog_name

####################### MAIN function #######################

def run_sourceextractor_subtractionimage_spreadmodel(sub, sci,input_list=False,
                        savecats_dir = None, catending = None,
                        sextractorloc='/opt/local/bin/source-extractor',
                        psfexloc='/opt/local/bin/psfex',
                        fwhm = 1.2, detect_minarea = 5, detect_thresh = 1.5, # Default settings in Source Extractor
                        verbose=False,quietmode=False,debugmode=False):

    if savecats_dir == None:
        savecats_dir = './'
    mkdir_ifnotexist(savecats_dir)

    # Print useful information if verbose
    if verbose:
        print('The Steps for 1 file at a time are:')
        print('1. Measure PSF of science input.')
        print('1a. Run Source Extractor')
        print('1b. Run PSFEx')
        print('2. Run Source Extractor with spread_model as a parameter output on subtraction input, using science PSF')
        print('#####################################')
        print('###### RUNNING SOURCE EXTRACTOR #####')
        print('#####################################\n')

    # Do Source Extraction

    if input_list == True:
        # Start catalog path list
        catfiles = []
        
        # Read in lists
        list_sub = [list(x) for x in ascii.read(sub,format='no_header')]
        list_sub = np.transpose(list_sub)[0]
        list_sci = [list(x) for x in ascii.read(sci,format='no_header')]
        list_sci = np.transpose(list_sci)[0]

        # For each input sub/sci pair, get PSF of sci, then get cat of sub
        for sub_,sci_ in zip(list_sub,list_sci):
            cat = run_sourceextractor_subtractionimage_spreadmodel_(sub_,sci_,
                                                                    sextractorloc=sextractorloc,
                                                                    psfexloc=psfexloc,
                                                                    savecats_dir = savecats_dir, catending = catending,
                                                                    fwhm=fwhm, detect_thresh=detect_thresh,detect_minarea=detect_minarea,
                                                                    verbose=verbose,quietmode=quietmode,debugmode=debugmode)     
            catfiles.append(cat)       

    elif input_list == False:
        # get PSF of sci, then get cat of sub
        cat = run_sourceextractor_subtractionimage_spreadmodel_(sub,sci,
                                                                sextractorloc=sextractorloc,
                                                                psfexloc=psfexloc,
                                                                savecats_dir = savecats_dir, catending = catending,
                                                                fwhm=fwhm, detect_thresh=detect_thresh,detect_minarea=detect_minarea,
                                                                verbose=verbose,quietmode=quietmode,debugmode=debugmode)

    # Return output cat file path(s)
    if input_list == True:
        return catfiles
    else:
        return cat

############################################################################
####################### BODY OF PROGRAM STARTS HERE ########################
#######################              MAIN           ########################
############################################################################

if __name__ == "__main__":

    arguments = docopt.docopt(__doc__)

    # Mandatory argument
    sub = arguments['<sub>']
    sci = arguments['<sci>']

    # Non-mandatory options without arguments
    quietmode       = arguments['--quiet']
    verbose         = arguments['--verbose']
    debugmode       = arguments['--debug']
    savecats_dir    = arguments['--savecats']
    catending       = arguments['--catending']
    sextractorloc   = arguments['--sextractor']
    psfexloc        = arguments['--psfex']
    fwhm            = arguments['--fwhm']
    detect_thresh   = arguments['--detect_thresh'] 
    detect_minarea  = arguments['--detect_minarea']
    input_list      = arguments['--list'] 

    # Calculate
    _ = run_sourceextractor_subtractionimage_spreadmodel(sub,sci,input_list=input_list,
                                                            sextractorloc=sextractorloc,
                                                            psfexloc=psfexloc,
                                                            savecats_dir = savecats_dir, catending = catending,
                                                            fwhm=fwhm, detect_thresh=detect_thresh,detect_minarea=detect_minarea,
                                                            verbose=verbose,quietmode=quietmode,debugmode=debugmode)
