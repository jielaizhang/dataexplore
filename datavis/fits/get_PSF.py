#!/usr/bin/env python

""" get_PSF.py -- Input fits file, run source extractor and psfex to get the PSF. Outputs psf fits file, and .psf file for input into Source Extractor to do model-based star-galaxy separation using SPREAD_MODEL. 
 
Usage: get_PSF.py [-h] [-v] [--debug] [-q] [-s LOC] [--psfex LOC] <fitsfile>

Arguments:
    fitsfile (string)
        fitsfile

Options:
    -h, --help                              Show this screen
    -v, --verbose                           Show extra information [default: False]   
    -q, --quietmode                         Minimize print to screen. This is useful when this function is called in another function. [default: False]  
    --debug                                 Output more for debugging [default: False]
    -s LOC, --sextractor LOC                SExtractor location [default: /opt/local/bin/source-extractor]   
    --psfex LOC                             PSFex location [default: /opt/local/bin/psfex]

Examples:
"""
import docopt
import astropy.io.fits as fits
import sys, os, ntpath
from pathlib import Path
import shutil
import subprocess


__author__      = "Jielai Zhang"
__license__     = "MIT"
__version__     = "1.0.1"
__date__        = "2020-07-02"

#########################################
# ======= House Keeping Functions =======
#########################################
'''These functions standardise verbose, debug printing'''
def print_verbose_string(printme,verbose=False,underscores=False):
    if verbose:
        if underscores:
            print("-" * len(f"VERBOSE: {printme}"),file=sys.stdout)
            print(f"VERBOSE: {printme}",file=sys.stdout)
            print("-" * len(f"VERBOSE: {printme}"),file=sys.stdout)
        else:
            print(f"VERBOSE: {printme}",file=sys.stdout)

def print_debug_string(printme,debugmode=False,underscores=False):
    if debugmode:
        if underscores:
            print("-" * len(f"DEBUG  : {printme}"),file=sys.stdout)
            print(f"DEBUG  : {printme}",file=sys.stdout)
            print("-" * len(f"DEBUG  : {printme}"),file=sys.stdout)
        else:
            print(f"DEBUG  : {printme}",file=sys.stdout)

'''These functions help organise or delete files'''
def clearit(fnames,debugmode=False):
    for fname in fnames:
        if os.path.isfile(fname):
            os.remove(fname)
            printme = f'Temporary file deleted: {fname}.'
            print_debug_string(printme,debugmode=debugmode)
    return None

##################################################
# ======= Source Extractor Required Inputs =======
##################################################
'''Source extractor is used to produce catalog taken in by
PSFex to create the PSF.fits file.
'''
params="""X_IMAGE
Y_IMAGE
FLUX_RADIUS
FLUX_APER(1)
FLUXERR_APER(1)
ELONGATION
FLAGS
SNR_WIN
VIGNET(35,35)"""

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

#################################
# ======= Other Functions =======
#################################

def run_SExtractor(image_name,verbose=False,debugmode=False,sextractorloc='/opt/local/bin/source-extractor'):

    # Create temporary directory
    if not os.path.exists('tmp_source_extractor_files'):
        os.makedirs('tmp_source_extractor_files')    
    else:
        print('./tmp_source_extractor_files directory existed already')

    # Set name of catalog output by Source Extractor and required by PSFex
    catalog_name = os.path.abspath(image_name.split('.fits')[0]+'_psfex.cat')

    # Create default Source Extractor config file, so that no error is printed
    config_name = './tmp_source_extractor_files/SEx.sex'
    command      = f'{sextractorloc} -d > {config_name}'
    subprocess.call(command,shell=True)
    printme = f'Temporarily saved: {config_name}.'
    print_debug_string(printme,debugmode=debugmode)

    # Names of required config files for Source Extractor
    params_name = './tmp_source_extractor_files/SEx.param'
    conv_name   = './tmp_source_extractor_files/default.conv'
    nnw_name    = './tmp_source_extractor_files/default.nnw'

    # Write temporary config files for Source Extractor
    configs = zip([params_name,conv_name,nnw_name],[params,default_conv,default_nnw])
    for fname,fcontent in configs:
        fout = open(fname,'w')
        fout.write(fcontent)
        fout.close()

    # Run Source Extractor
    if verbose:
        verbose_type = 'NORMAL'
    else:
        verbose_type = 'QUIET'
    print_verbose_string('SExtracting...',verbose=verbose)
    command = f'{sextractorloc} -c {config_name} '\
              f'-CATALOG_NAME {catalog_name} -CATALOG_TYPE FITS_LDAC '\
              f'-PARAMETERS_NAME {params_name} '\
              f'-FILTER_NAME {conv_name} '\
              f'-STARNNW_NAME {nnw_name} '\
              f'-VERBOSE_TYPE {verbose_type} '\
              f'{image_name}'

    print_verbose_string('Running this command:',verbose=verbose)
    print_verbose_string(command+'\n',verbose=verbose)

    subprocess.call(command,shell=True)

    # Clear unnecessary files
    clearit([config_name,params_name,conv_name,nnw_name])

    # Remove temp directory if it is empty
    try:
        os.rmdir('tmp_source_extractor_files')
    except OSError as ex:
        if ex.errno == errno.ENOTEMPTY:
            print(f"directory not empty, cannot delete: ./tmp_source_extractor_files")

    # Print what is saved:
    print_debug_string(f'SAVED  : {os.path.abspath(catalog_name)}',debugmode=debugmode)
    
    return catalog_name

def run_PSFex(psfex_cat,psfexloc='opt/local/bin/psfex',verbose=False,debugmode=False):

    # Create temporary directory
    if not os.path.exists('tmp_PSFEx_files'):
        os.makedirs('tmp_PSFEx_files')    
    else:
        print('./tmp_PSFEx_files directory existed already')

    # Create default PSFEx config file, so that no error is printed
    config_name = './tmp_PSFEx_files/PSFEx.sex'
    command      = f'{psfexloc} -d > {config_name}'
    subprocess.call(command,shell=True)
    printme = f'Temporarily saved: {config_name}.'
    print_debug_string(printme,debugmode=debugmode)

    #Run PSFex to compute PSF
    if verbose:
        verbose_type = 'NORMAL'
    else:
        verbose_type = 'QUIET'
    try:
        command = f"psfex {psfex_cat} -c {config_name} "\
                  f"-CHECKIMAGE_TYPE PROTOTYPES -CHECKIMAGE_NAME proto.fits "\
                  f"-VERBOSE_TYPE {verbose_type}"
        print_verbose_string(f'Executing command: {command}',verbose=verbose)
        rval = subprocess.run(command.split(), check=True)
        if not quietmode:
            print('Success!')
    except subprocess.CalledProcessError as err:
        print(f'command: {command}')
        sys.exit('ERROR  : Could not run psfex with exit error %s'%err)

    # Determine psf output name
    psfex_namestub     = Path(ntpath.basename(psfex_cat)).stem
    psfex_dir          = ntpath.dirname(os.path.abspath(psfex_cat))
    psfex_psf_file     = psfex_dir + os.sep + psfex_namestub + '.psf'

    # Determine psf fits file checkimage output, and move it to where psf output is (same dir as input file)
    psfex_proto_file      = psfex_dir + os.sep + 'proto_'+ psfex_namestub + '.fits'
    psfex_proto_outputloc = '.'+os.sep+'proto_'+psfex_namestub+'.fits'
    shutil.move(psfex_proto_outputloc,psfex_proto_file)

    # Clean temp files
    xml_file = 'psfex.xml'
    clearit([config_name,xml_file])

    # Remove temp directory if it is empty
    try:
        os.rmdir('tmp_PSFEx_files')
    except OSError as ex:
        if ex.errno == errno.ENOTEMPTY:
            print(f"directory not empty, cannot delete: ./tmp_PSFEx_files")

    return psfex_psf_file, psfex_proto_file

##############################################################
####################### Main Function ########################
##############################################################

def get_PSF(fitsfile,
            sextractorloc='/opt/local/bin/source-extractor',psfexloc='/opt/local/bin/psfex',
            verbose=False,debugmode=False,quietmode=False):
    '''Input fits file, run source extractor and psfex to get the PSF. 
    Outputs psf fits file, and .psf file for input into Source Extractor to do model-based star-galaxy separation using SPREAD_MODEL.
    '''

    print_debug_string(f'sextractor loc is specified in get_PSF.py to be: {sextractorloc}',debugmode=debugmode)
    print_debug_string(f'PSFEx loc is specified in get_PSF.py to be: {psfexloc}',debugmode=debugmode)

    # Run Source Extractor on fits file
    print_verbose_string('Running Source Extractor',underscores=True,verbose=verbose)
    psfex_cat = run_SExtractor(fitsfile, sextractorloc=sextractorloc,verbose=verbose,debugmode=debugmode)

    # Run PSFex on source extractor output
    print_verbose_string('Running PSFEx',underscores=True,verbose=verbose)
    psfex_psf_file, psfex_proto_file = run_PSFex(psfex_cat,psfexloc=psfexloc,verbose=verbose,debugmode=debugmode)

    # Print what was saved if not quietmode
    if not quietmode:
        print(f'SAVED  : {psfex_psf_file}')
        print(f'SAVED  : {psfex_proto_file}')

    # If not debug mode, remove following files
    if not debugmode:
        clearit([psfex_cat],debugmode=debugmode)
    else:
        print_debug_string(f'SAVED  : {psfex_cat}',debugmode=debugmode)

    return psfex_psf_file, psfex_proto_file

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
    sextractorloc       = arguments['--sextractor']
    psfexloc            = arguments['--psfex']
    fitsfile            = arguments['<fitsfile>']

    _,_ = get_PSF(fitsfile,
                    sextractorloc=sextractorloc,psfexloc=psfexloc,
                    verbose=verbose,debugmode=debugmode,quietmode=quietmode)
