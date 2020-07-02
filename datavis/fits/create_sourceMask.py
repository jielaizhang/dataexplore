#!/usr/bin/python

"""create_sourceMask.py -- Create a fits file that can act as a source mask for input fits image. Mask will have zeros and ones, ones are the locations where there are stars. The mask is created by growing the source extractor segmentation map. 

Required input 
    fitsimage - image for which a source mask is to be created.

Usage:
    create_sourceMask.py [-h] [-v] [--debug] [-o LOC] [-s LOC] <fitsimage>

Options:
    -h, --help                      Print this screen.
    -v, --verbose                   Print extra information [default: False]
    --debug                         Print extra extra information and save extra files [default: False]
    -o LOC, --out LOC               Output mask fits file here [default: ./source_mask.fits]
    -s LOC, --sextractor LOC        SExtractor location [default: /opt/local/bin/source-extractor]   

Example:
    python create_sourceMask.py -v fitsimage.fits
"""

import docopt
import sys, os
import astropy.io.fits as fits
import subprocess
from scipy import ndimage

##################################################
# ======= Source Extractor Required Inputs =======
##################################################
'''Source extractor is used to produce the star mask. 
The segmentation map output by source extractor is grown
a little larger and used as the aggressive source mask.
'''

sextractor_params = """NUMBER
FLUX_AUTO
FLUXERR_AUTO
FLUX_APER
FLUXERR_APER
X_IMAGE
Y_IMAGE
X_WORLD
Y_WORLD
FLUX_RADIUS
FLAGS
CLASS_STAR
BACKGROUND
ELLIPTICITY
FWHM_IMAGE
"""

sextractor_config = """
    ANALYSIS_THRESH 3
        BACK_FILTERSIZE 3
        BACKPHOTO_TYPE LOCAL
        BACK_SIZE 32
        CATALOG_NAME test.cat
        CATALOG_TYPE ASCII_HEAD
        CHECKIMAGE_TYPE SEGMENTATION
        CHECKIMAGE_NAME {check_name}
        CLEAN Y
        CLEAN_PARAM 1.
        DEBLEND_MINCONT 0.001
        DEBLEND_NTHRESH 32
        DETECT_MINAREA 5
        DETECT_THRESH 3
        DETECT_TYPE CCD
        FILTER Y
        FILTER_NAME {filter_name}
        FLAG_IMAGE flag.fits
        GAIN 1.0
        MAG_GAMMA 4.
        MAG_ZEROPOINT 0.0
        MASK_TYPE CORRECT
        MEMORY_BUFSIZE 1024
        MEMORY_OBJSTACK 3000
        MEMORY_PIXSTACK 300000
        PARAMETERS_NAME {parameters_name}
        PHOT_APERTURES 5
        PHOT_AUTOPARAMS 2.5, 3.5
        PIXEL_SCALE 2.85
        SATUR_LEVEL 50000.
        SEEING_FWHM 2.5
        STARNNW_NAME {starnnw_name}
        VERBOSE_TYPE {verbose_type}
"""

default_conv = """CONV NORM
# 3x3 ``all-ground'' convolution mask with FWHM = 2 pixels.
1 2 1
2 4 2
1 2 1
"""

default_nnw = """NNW
# Neural Network Weights for the SExtractor star/galaxy classifier (V1.3)
# inputs:       9 for profile parameters + 1 for seeing.
# outputs:      ``Stellarity index'' (0.0 to 1.0)
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
 1.00000e+00
"""

#########################################
# ======= House Keeping Functions =======
#########################################
'''These functions standardise verbose, debug printing'''
def print_verbose_string(printme,verbose=False,underscores=False):
    if verbose:
        if underscores:
            print("-" * len(f"VERBOSE: {printme}"),file=sys.stdout)
        else:
            print(f"VERBOSE: {printme}",file=sys.stdout)

def print_debug_string(printme,debugmode=False,underscores=False):
    if debugmode:
        if underscores:
            print("-" * len(f"DEBUG  : {printme}"),file=sys.stdout)
        else:
            print(f"DEBUG  : {printme}",file=sys.stdout)

'''These functions help organise or delete files'''
def clearit(fname):
    if os.path.isfile(fname):
        os.remove(fname)
    return None

#################################
# ======= Other Functions =======
#################################

def run_SExtractor(image_name,verbose=False,sextractorloc='/opt/local/bin/source-extractor'):

    'Create temporary directory'
    if not os.path.exists('tmp_source_extractor_files'):
        os.makedirs('tmp_source_extractor_files')    
    else:
        print('./tmp_source_extractor_files directory existed already')

    'Names of required config files'
    sextractor_config_name = './tmp_source_extractor_files/scamp.sex'
    params_name = './tmp_source_extractor_files/scamp.param'
    conv_name = './tmp_source_extractor_files/default.conv'
    nnw_name = './tmp_source_extractor_files/default.nnw'

    catalog_name = image_name.split('.fits')[0]+'_bkg.cat'
    check_name = image_name.split('.fits')[0]+'_bkg_segmap.fits'

    if verbose:
        verbose_type = 'NORMAL'
    else:
        verbose_type = 'QUIET'

    'Stick content in config files'
    configs = zip([sextractor_config_name,params_name,conv_name,nnw_name],[sextractor_config,sextractor_params,default_conv,default_nnw])

    for fname,fcontent in configs:
        fout = open(fname,'w')

        if 'scamp.sex' in fname:
            fout.write(fcontent.format(filter_name=conv_name,
            parameters_name=params_name,starnnw_name=nnw_name,
            verbose_type=verbose_type,check_name=check_name))

        else:
            fout.write(fcontent)

        fout.close()

    if verbose:
        print('SExtracting...')


    'SExtractor command'
    command = sextractorloc + ' -c {config} -CATALOG_NAME {catalog} {image}'.format(config=sextractor_config_name,catalog=catalog_name,image=image_name)

    if verbose:
        print('Running this command:')
        print(command+'\n')

    'Run SExtractor'
    subprocess.call(command,shell=True)

    'Clear unnecessary files'
    for fname in [sextractor_config_name,params_name,conv_name,nnw_name,catalog_name]:
        clearit(fname)
    
    'Remove temp directory if it is empty'
    try:
        os.rmdir('tmp_source_extractor_files')
    except OSError as ex:
        if ex.errno == errno.ENOTEMPTY:
            print(f"directory not empty, cannot delete: ./tmp_source_extractor_files")

    return check_name


###############################
# ======= Main Function =======
###############################

def create_sourceMask(fitsimage,output=False,sextractorloc='/opt/local/bin/source-extractor',verbose=False,debugmode=False):
    '''Create a fits file that can act as a source mask for input fits image. 
    Mask will have zeros and ones, ones are the locations where there are stars. 
    The mask is created by growing the source extractor segmentation map.'''

    print_debug_string(f'sextractor loc is specified in create_sourceMask to be: {sextractorloc}',debugmode=debugmode)

    # If write fits file specified, check that file doesn't already exist. If yes, exit.
    if output:
        if os.path.exists(output):
            printme = f'ERROR  : cannot overwrite existing file: {output}. Delete/rename/specify a different output name before trying again.'
            sys.exit(printme)

    # Read in fits image
    d = fits.getdata(fitsimage)
    
    # Run Source Extractor, save segmentation map
    segname = run_SExtractor(fitsimage, sextractorloc=sextractorloc,verbose=verbose)
    printme = f'SAVED: {segname}'
    print_debug_string(printme,debugmode=debugmode)

    # Read in segmentation map and grow it so sources are definitely covered
    # Seg map has each segmentation (or source) marked with a different number.
    segmap  = fits.getdata(segname)
    grownsegmap = ndimage.grey_dilation(segmap,size=(15,15))

    # Change all non-zero pixels in grown seg map to 1s
    grownsegmap[grownsegmap!=0] = 1
    source_mask = grownsegmap

    # Save mask if option specifies
    if output:
        fits.writeto(output,source_mask) 
        printme = f'SAVED  : {output}'
        print(printme)

    # Remove segmentation map save if not verbose
    # No need to print removal, as if not debug, print save doesn't happen either.
    if not debugmode:
        clearit(segname)
        
    return source_mask

###########################################################
###################### Start of main ######################
###########################################################
if __name__=='__main__':

    # Import arguments
    arguments = docopt.docopt(__doc__)
    
    fitsimage           = arguments['<fitsimage>']
    verbose             = arguments['--verbose']
    debugmode           = arguments['--debug']
    output              = arguments['--out']
    sextractorloc       = arguments['--sextractor']

    if debugmode:
        print(arguments)

    _ = create_sourceMask(  fitsimage, output = output, sextractorloc = sextractorloc, verbose = verbose, debugmode = debugmode)


