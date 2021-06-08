#!/usr/bin/env python

"""write_sourceextractor_files.py -- input directory and write out .nnw, .conv, .config. A .param is written, but check it so that it has the parameters you need, and edit as appropriate.

Usage: write_sourceextractor_files [-h] [-q] [--debug] [-s LOC] [--name STRING] [--spreadmodel] <outdirectory>

Options:
    -h, --help                  Show this screen
    -q, --quiet                 Quiet mode, suppress printout of parameters measured.
    --name STRING               Output files are outdirectory/STRING_default.conv etc, if not set, output outdirectory/default.conv
    -s LOC, --sextractor LOC    Location of source extractor [default: /opt/local/bin/source-extractor] 
    --spreadmodel               If true, parameters output contains spread_model and mag_model and mag model err [default: False]
    --debug                     Print source extractor location if true. [default: False]

Examples:
    from misc.write_sourceextractor_files import write_sourceextractor_files
    (nnw_path, 
    conv_path, 
    params_path, 
    config_path) = write_sourceextractor_files(outdirectory,file_start_string,spreadmodel=True,
                                    sextractorloc=sextractorloc,quietmode=quietmode)

"""

import docopt, os
import subprocess
import sys

####################### Source Extractor Files Templates #######################

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

f_params='''NUMBER
FLUX_AUTO
FLUXERR_AUTO
FLUX_RADIUS
FLUX_APER
FLUXERR_APER
X_IMAGE
Y_IMAGE
X_WORLD
Y_WORLD
FLUX_RADIUS
FLAGS
CLASS_STAR
MAG_AUTO
MAGERR_AUTO
MAG_ISO
MAGERR_ISO
BACKGROUND
A_IMAGE
B_IMAGE
THETA_IMAGE
THETA_SKY
ISOAREA_IMAGE
FWHM_IMAGE
FWHM_WORLD
ISOAREAF_IMAGE
ELLIPTICITY
'''

f_params_spreadmodel='''NUMBER
FLUX_AUTO
FLUXERR_AUTO
FLUX_RADIUS
FLUX_APER
FLUXERR_APER
X_IMAGE
Y_IMAGE
X_WORLD
Y_WORLD
FLUX_RADIUS
FLAGS
CLASS_STAR
MAG_AUTO
MAGERR_AUTO
MAG_ISO
MAGERR_ISO
BACKGROUND
A_IMAGE
B_IMAGE
THETA_IMAGE
THETA_SKY
ISOAREA_IMAGE
FWHM_IMAGE
FWHM_WORLD
ISOAREAF_IMAGE
ELLIPTICITY
SPREAD_MODEL
MAG_MODEL
MAGERR_MODEL
'''

def mkdir_ifnotexist(directory):
    if not os.path.isdir(directory):
        os.makedirs(directory)
    return None

def write_sourceextractor_files(outdirectory,file_start_string,spreadmodel=False,
                                sextractorloc='/opt/local/bin/source-extractor',quietmode=False,debug=False):

    if debug:
        print(f'This is the Source Extractor Loc used: {sextractorloc}')
    mkdir_ifnotexist(outdirectory)    

    if file_start_string:
        nnw_path    = outdirectory+os.sep+file_start_string+'_default.nnw'
        conv_path   = outdirectory+os.sep+file_start_string+'_default.conv'
        params_path = outdirectory+os.sep+file_start_string+'_default.param'
        config_path = outdirectory+os.sep+file_start_string+'_default.config'
    else:
        nnw_path    = outdirectory+os.sep+'default.nnw'
        conv_path   = outdirectory+os.sep+'default.conv'
        params_path = outdirectory+os.sep+'default.param'
        config_path = outdirectory+os.sep+'default.config'
        
    
    fp = open(nnw_path, "w")
    fp.write(f_nnw)
    fp.close()

    fp = open(conv_path, "w")
    fp.write(f_conv)
    fp.close()
    
    if spreadmodel:
        fp = open(params_path, "w")
        fp.write(f_params_spreadmodel)
        fp.close()
    else:
        fp = open(params_path, "w")
        fp.write(f_params)
        fp.close()

    # Create temporary default.sex file
    try:
        command = sextractorloc+' -d > '+  config_path
        if debug:
            print(f'This is the source extractor command: {command}')
        subprocess.check_call(command,shell=True)
        if not quietmode:
            print('Saved: ',config_path)
    except: 
        print('ERROR: NOT SAVED: ',config_path)
        print(f'ERROR: Is this the right source extractor loc: {sextractorloc}')
        print(f'ERROR: If not, add this option: sextractorloc="path", or -s path in command line')
        sys.exit('Cannot continue')

    if not quietmode:
        print('Saved: ',nnw_path)
        print('Saved: ',conv_path)
        print('Saved: ',params_path)
        print('Check ',params_path,' to see if it has the parameters you want. Edit for your purposes.')

    return nnw_path, conv_path, params_path, config_path

############################################################################
####################### BODY OF PROGRAM STARTS HERE ########################
#######################              MAIN           ########################
############################################################################

if __name__ == "__main__":

    arguments = docopt.docopt(__doc__)

    # Mandatory argument
    outdirectory        = arguments['<outdirectory>']

    # Non-mandatory options without arguments
    quietmode           = arguments['--quiet']
    sextractorloc       = arguments['--sextractor']
    file_start_string   = arguments['--name']
    spreadmodel         = arguments['--spreadmodel']
    debug               = arguments['--debug']

    # Calculate
    _ = write_sourceextractor_files(outdirectory,file_start_string,spreadmodel=spreadmodel,
                                    sextractorloc=sextractorloc,quietmode=quietmode,debug=debug)
