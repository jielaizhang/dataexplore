#!/usr/bin/env python

"""measure_psf.py -- read in fits files and use PSFEx to measure the average PSF of the images.  

Usage: measure_psf [-q] [-h] [-v] [-s LOC] [-p LOC] [-o DIR] [--fitssave] <fitsfiles>...

Output:
    For each inputfile.fits output inputfile.psf, a single extension file of the PSF and nothing else.

Options:
    -h, --help                  Show this screen
    -q, --quiet                 Quiet mode, suppress printout FWHM measured.
    -v, --verbose               Show extra information [default: False]     
    -s LOC, --sextractor LOC    Location of source extractor [default: /opt/local/bin/source-extractor] 
    -p LOC, --psfex LOC         Location of psfex [default: /opt/local/bin/psfex]
    -o DIR, --out DIR           Location of PSF output. [default: ./]
    --fitssave                  Output not just .psf file, but also [default: False]

Examples:
    bash: python measure_psf.py -v -o ../psfs one.fits two.fits 
    python: from datastats.measure_psf import measure_psf
            PSFs = measure_psf(fitsfiles,sextractorloc='/opt/local/bin/source-extractor,psfloc='/opt/local/bin/psfex',verbose=False|True)
"""

import docopt, os
import astropy.io.fits as fits
import subprocess
from astropy.io import ascii
import pandas as pd
import numpy as np
import ntpath
from pathlib import Path

####################### Source Extractor Files Templates #######################
conv_name = "./temp_default.conv"
params_name = "./temp_params.txt"
config_name = "./temp_default.sex"

f_conv = '''CONV NORM
# 3x3 ``all-ground'' convolution mask with FWHM = 2 pixels.
1 2 1
2 4 2
1 2 1'''

f_params='''X_IMAGE
Y_IMAGE
FLUX_RADIUS
FLUX_APER(1)
FLUXERR_APER(1)
ELONGATION
FLAGS
SNR_WIN
VIGNET(35,35)
'''

####################### Auxiliary Functions #######################
def create_temp_files_for_SourceExtractor(f_conv,f_params,conv_name,params_name,
                                          config_name,
                                          sextractorloc='/opt/local/bin/source-extractor'):

    fp = open(params_name, "w")
    fp.write(f_params)
    fp.close()

    fp = open(conv_name, "w")
    fp.write(f_conv)
    fp.close()

    command = sextractorloc+' -d > '+ config_name 
    subprocess.call(command,shell=True)

    return None

def remove_temp_files(fs):
    for f in fs:
        os.remove(f)
    return None

####################### calculate_FWHM function #######################

def measure_psf(fitsfiles, outdir='./', savepsffits=False,
                sextractorloc='/opt/local/bin/source-extractor',
                psfexloc='/opt/local/bin/psfex',
                verbose=False,quietmode=False):

    # Create temporary files required for Source Extractor
    create_temp_files_for_SourceExtractor(f_conv,f_params,conv_name,params_name,
                                          config_name,sextractorloc=sextractorloc)

    # Set Source Extractor verbose type
    if verbose:
        VERBOSE_TYPE = 'NORMAL'
    else:
        VERBOSE_TYPE = 'QUIET'

    PSFs = []

    # Run Source Extractor and PSFEx one image at a time
    for f in fitsfiles:

        # Name temporary Source Extractor Output (this enables PSF output to be named appropriately)
        cat_out_name_temp = f.replace('.fits','.psfcat')

        # Print separator
        if verbose:
            print('==============================================')

        # Run Source Extractor on image
        try:
            command = (f"{sextractorloc} -c {config_name} " 
                      f"-MAG_ZEROPOINT 25.0 "
                      f"-CATALOG_TYPE FITS_LDAC "
                      f"-FILTER_NAME {conv_name} "
                      f"-VERBOSE_TYPE {VERBOSE_TYPE} "
                      f"-PARAMETERS_NAME {params_name} " 
                      f"-CATALOG_NAME {cat_out_name_temp} "
                      f"{f}")
            if verbose:
                print('Executing command: %s\n' % command)
            rval = subprocess.run(command.split(), check=True)
            if verbose:
                print('Success!')
        except subprocess.CalledProcessError as err:
            print('\nCould not run SExtractor with exit error %s\n'%err)
            print('Command used:\n%s\n'%command)

        # Run PSFEx on image
        try:
            if savepsffits:
                command = (f"{psfexloc} "
                           f"-PSF_DIR {outdir} "
                           f"-CHECKIMAGE_TYPE PROTOTYPES "
                           f"-CHECKIMAGE_NAME  proto.fits "
                           f"-PSF_SUFFIX .psf "
                           f"{cat_out_name_temp}")
            else:
                command = (f"{psfexloc} "
                           f"-PSF_DIR {outdir} "
                           f"-CHECKIMAGE_TYPE NONE "
                           f"-CHECKIMAGE_NAME  NONE "
                           f"-PSF_SUFFIX .psf "
                           f"{cat_out_name_temp}")
            if verbose:
                print('Executing command: %s\n' % command)
            rval = subprocess.run(command.split(), check=True)
            if verbose:
                print('Success!')
        except subprocess.CalledProcessError as err:
            print('Could not run psfex with exit error %s'%err) 

        remove_temp_files([cat_out_name_temp])    

        # Save psf as a fits file if requested  
        if savepsffits:
            f_filestub = Path(ntpath.basename(f)).stem
            proto_file = './proto_'+f_filestub+'.fits'
            # Read out the PSF 
            d_psf      = fits.getdata(proto_file)[0:25,0:25]
            # Save 
            f_psf      = outdir+'/'+f_filestub+'_psf.fits'
            fits.writeto(f_psf,d_psf)
            # Remove proto file
            os.remove(proto_file)
            
    # Remove temporary files required for Source Extractor
    remove_temp_files([params_name,conv_name,config_name,'psfex.xml'])
    
    return PSFs

####################### BODY OF PROGRAM STARTS HERE ########################

if __name__ == "__main__":

    arguments = docopt.docopt(__doc__)

    # Mandatory argument
    fitsfiles = arguments['<fitsfiles>']

    # Non-mandatory options without arguments
    quietmode       = arguments['--quiet']
    verbose         = arguments['--verbose']
    sextractorloc   = arguments['--sextractor']
    psfexloc        = arguments['--psfex']
    outdir          = arguments['--out']
    savepsffits     = arguments['--fitssave']
    
    # Calculate
    _ = measure_psf(fitsfiles,outdir=outdir,
                    savepsffits = savepsffits,
                    sextractorloc=sextractorloc,
                    psfexloc=psfexloc,
                    verbose=verbose,quietmode=quietmode)
