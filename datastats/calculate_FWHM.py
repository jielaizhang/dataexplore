#!/usr/bin/env python

"""calculate_FWHM.py -- read in fits files and calculate the average FWHM of the image by using source extractor. Conditions for which sources to include in the average calculation can be entered as options. 

Usage: calculate_FWHM [-q] [-h] [-v] [-s LOC] <fitsfiles>...

Options:
    -h, --help                  Show this screen
    -q, --quiet                 Quiet mode, suppress printout FWHM measured.
    -v, --verbose               Show extra information [default: False]     
    -s LOC, --sextractor LOC    Location of source extractor [default: /opt/local/bin/source-extractor] 

Examples:
    bash: python calculate_FWHM.py -v one.fits two.fits 
    python: from datastats.calculate_FWHM import calculate_FWHM
            FWHMs = calculate_FWHM(fitsfiles,sextractorloc='/opt/local/bin/source-extractor,verbose=False|True)
"""

import docopt, os
import astropy.io.fits as fits
import subprocess
from astropy.io import ascii
import pandas as pd
import numpy as np

####################### Source Extractor Files Templates #######################
conv_name = "./temp_default.conv"
nnw_name = "./temp_default.nnw"
params_name = "./temp_params.txt"
config_name = "./temp_default.sex"

f_conv = '''CONV NORM
# 3x3 ``all-ground'' convolution mask with FWHM = 2 pixels.
1 2 1
2 4 2
1 2 1'''

f_nnw = '''NNW
# Neural Network Weights for the SExtractor star/galaxy classifier (V1.3)
# inputs:   9 for profile parameters + 1 for seeing.
# outputs:  ``Stellarity index'' (0.0 to 1.0)
# Seeing FWHM range: from 0.025 to 5.5'' (images must have 1.5 < FWHM < 5 pixels)
# Optimized for Moffat profiles with 2<= beta <= 4.

 3 10 10  1

-1.56604e+00 -2.48265e+00 -1.44564e+00 -1.24675e+00 -9.44913e-01 -5.22453e-01  4.61342e-02  8.31957e-01  2.15505e+00  2.64769e-01
 3.03477e+00  2.69561e+00  3.16188e+00  3.34497e+00  3.51885e+00  3.65570e+00  3.74856e+00  3.84541e+00  4.22811e+00  3.27734e+00

-3.22480e-01 -2.12804e+00  6.50750e-01 -1.11242e+00 -1.40683e+00 -1.55944e+00 -1.84558e+00 -1.18946e-01  5.52395e-01 -4.36564e-01 -5.30052e+00
 4.62594e-01 -3.29127e+00  1.10950e+00 -6.01857e-01  1.29492e-01  1.42290e+00  2.90741e+00  2.44058e+00 -9.19118e-01  8.42851e-01 -4.69824e+00
-2.57424e+00  8.96469e-01  8.34775e-01  2.18845e+00  2.46526e+00  8.60878e-02 -6.88080e-01 -1.33623e-02  9.30403e-02  1.64942e+00 -1.01231e+00
 4.81041e+00  1.53747e+00 -1.12216e+00 -3.16008e+00 -1.67404e+00 -1.75767e+00 -1.29310e+00  5.59549e-01  8.08468e-01 -1.01592e-02 -7.54052e+00
 1.01933e+01 -2.09484e+01 -1.07426e+00  9.87912e-01  6.05210e-01 -6.04535e-02 -5.87826e-01 -7.94117e-01 -4.89190e-01 -8.12710e-02 -2.07067e+01
-5.31793e+00  7.94240e+00 -4.64165e+00 -4.37436e+00 -1.55417e+00  7.54368e-01  1.09608e+00  1.45967e+00  1.62946e+00 -1.01301e+00  1.13514e-01
 2.20336e-01  1.70056e+00 -5.20105e-01 -4.28330e-01  1.57258e-03 -3.36502e-01 -8.18568e-02 -7.16163e+00  8.23195e+00 -1.71561e-02 -1.13749e+01
 3.75075e+00  7.25399e+00 -1.75325e+00 -2.68814e+00 -3.71128e+00 -4.62933e+00 -2.13747e+00 -1.89186e-01  1.29122e+00 -7.49380e-01  6.71712e-01
-8.41923e-01  4.64997e+00  5.65808e-01 -3.08277e-01 -1.01687e+00  1.73127e-01 -8.92130e-01  1.89044e+00 -2.75543e-01 -7.72828e-01  5.36745e-01
-3.65598e+00  7.56997e+00 -3.76373e+00 -1.74542e+00 -1.37540e-01 -5.55400e-01 -1.59195e-01  1.27910e-01  1.91906e+00  1.42119e+00 -4.35502e+00

-1.70059e+00 -3.65695e+00  1.22367e+00 -5.74367e-01 -3.29571e+00  2.46316e+00  5.22353e+00  2.42038e+00  1.22919e+00 -9.22250e-01 -2.32028e+00


 0.00000e+00 
 1.00000e+00'''

f_params='''FLUX_AUTO
CLASS_STAR
FWHM_IMAGE
ELLIPTICITY'''

####################### Auxiliary Functions #######################
def read_cat2pd(cat_file):
    cat = ascii.read(cat_file)

    df = pd.DataFrame({'CLASS_STAR': np.array(cat['CLASS_STAR']),
                       'ELLIPTICITY': np.array(cat['ELLIPTICITY']),
                       'FLUX_AUTO':np.array(cat['FLUX_AUTO']),
                       'FWHM_IMAGE':np.array(cat['FWHM_IMAGE']),
                      })
            
    return df

def create_temp_files_for_SourceExtractor(f_nnw,f_conv,f_params,
                                          nnw_name,conv_name,params_name):

    fp = open(nnw_name, "w")
    fp.write(f_nnw)
    fp.close()
    fp = open(conv_name, "w")
    fp.write(f_conv)
    fp.close()
    fp = open(params_name, "w")
    fp.write(f_params)
    fp.close()
    return None

def remove_temp_files(fs):
    for f in fs:
        os.remove(f)
    return None

####################### calculate_FWHM function #######################

def calculate_FWHM(fitsfiles,sextractorloc='/opt/local/bin/source-extractor',verbose=False,quietmode=False):

    # Create temporary files required for Source Extractor
    create_temp_files_for_SourceExtractor(f_nnw,f_conv,f_params,nnw_name,conv_name,params_name)

    # Create temporary default.sex file
    command = sextractorloc+' -d > '+ config_name # config name comes from global parameter set above.
    subprocess.call(command,shell=True)

    # Print useful information if verbose
    if verbose:
        print('CLASS_STAR  1=star\nEllipticity 0=round')

    # Set Source Extractor verbose type
    if verbose:
        VERBOSE_TYPE = 'NORMAL'
    else:
        VERBOSE_TYPE = 'QUIET'

    FWHMs = []

    # Run Source Extractor and Calculate Average FWHM one image at a time
    for f in fitsfiles:

        # Print separator
        if verbose:
            print('==============================================')

        # Run Source Extractor on image
        try:
            command = sextractorloc+" %s %s %s %s %s %s %s %s %s" % (f,
                                         '-c '+config_name,
                                         '-CATALOG_NAME ./temp.cat',
                                         '-CATALOG_TYPE ASCII_HEAD',
                                         '-PARAMETERS_NAME '+params_name,
                                         '-FILTER_NAME '+conv_name,
                                         '-STARNNW_NAME '+nnw_name,
                                         '-MAG_ZEROPOINT 25.0',
                                         '-VERBOSE_TYPE '+VERBOSE_TYPE
                                         )
            if verbose:
                print('Executing command: %s\n' % command)
            rval = subprocess.run(command.split(), check=True)
            if verbose:
                print('Success!')
        except subprocess.CalledProcessError as err:
            print('\nCould not run SExtractor with exit error %s\n'%err)
            print('Command used:\n%s\n'%command)

        # Read in catalog
        df = read_cat2pd('./temp.cat') 
        if verbose:
            print('\nFor file                     : ',f)
            print('Number of sources extracted  : ',len(df['CLASS_STAR']))
            print('Average CLASS_STAR           : %.5f'%np.average(df['CLASS_STAR']))
            print('Average ELLIPTICITY          : %.5f'%np.average(df['ELLIPTICITY']))

        # Pick catalog entries most likely to be stars
        # CLASS_STAR  1=star
        # Ellipticity 0=round
        df_stars = df[(df['ELLIPTICITY']<0.13) & (df['CLASS_STAR']>0.01)]
        
        # Calculate Average FWHM and print
        FWHM_average = np.average(df_stars['FWHM_IMAGE'])
        if not quietmode:
            print('The Average FWHM (pixels) is : %.5f for file %s'%(FWHM_average,f))
        FWHMs.append(FWHM_average)

    # Remove temporary files required for Source Extractor
    remove_temp_files([nnw_name,conv_name,params_name,config_name,'./temp.cat'])

    # Print useful information if verbose
    if verbose:
        print('\nCLASS_STAR  1=star\nEllipticity 0=round')
    
    return FWHMs

####################### BODY OF PROGRAM STARTS HERE ########################

if __name__ == "__main__":

    arguments = docopt.docopt(__doc__)

    # Mandatory argument
    fitsfiles = arguments['<fitsfiles>']

    # Non-mandatory options without arguments
    quietmode       = arguments['--quiet']
    verbose         = arguments['--verbose']
    sextractorloc   = arguments['--sextractor']
    
    # Calculate
    _ = calculate_FWHM(fitsfiles,sextractorloc=sextractorloc,verbose=verbose,quietmode=quietmode)
