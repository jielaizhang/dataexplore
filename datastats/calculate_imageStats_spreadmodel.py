#!/usr/bin/env python

"""calculate_imageStats_spreadmodel.py -- read in fits files and calculate the average FWHM, average ELLIPTICITY, # sources, background of the image by using source extractor. 
To do: write another one that is calculate_imageStats.py

Usage: calculate_imageStats_spreadmodel [-q] [--debug] [-h] [-v] [-s LOC] [-p LOC] [--minfwhm FLOAT] [-o FILE] <fitsfiles>...

Options:
    -h, --help                  Show this screen
    -q, --quiet                 Quiet mode, suppress printout of parameters measured.
    --debug                     Show extra information and save extra files [default: False]
    -v, --verbose               Show extra information [default: False]     
    -s LOC, --sextractor LOC    Location of source extractor [default: /opt/local/bin/source-extractor] 
    -p LOC, --psfex LOC         Location of PSFEx [default: /opt/local/bin/psfex]
    --minfwhm FLOAT             Minimum FWHM_IMAGE for this instrument, used for star selection. Default set for DECam [default: 2.05]
    -o FILE, --out FILE         Save stats on each input image in this file 

Examples:
    python calculate_imageStats_spreadmodel.py *fits -o ./outfile.txt
    from datastats.calculate_imageStats_spreadmodel import calculate_imageStats_spreadmodel
        (catted_fitsfiles,
        FWHMs,FWHM_stds,ELLIPs,ELLIP_stds,
        N_SRCs,N_SRCs_stars,BGRs,BGR_stds) = calculate_imageStats_spreadmodel(fitsfiles,outfile='./stats.txt',
                                                                                sextractorloc='/opt/local/bin/source-extractor',
                                                                                psfexloc='/opt/local/bin/psfex',
                                                                                minfwhm=2.05,
                                                                                verbose=False,quietmode=False,
                                                                                debugmode=False)
"""

import docopt, os
import astropy.io.fits as fits
import subprocess
from astropy.io import ascii
import pandas as pd
import numpy as np
import ntpath
import time

# Jielai Modules
from datastats.measure_psf import measure_psf
from datastats.calculate_skyStats2 import calculate_skyStats

####################### Source Extractor Files Templates #######################
conv_name = "./temp_default_stats.conv"
nnw_name = "./temp_default_stats.nnw"
params_name = "./temp_params_stats.txt"
config_name = "./temp_default_stats.sex"

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
ELLIPTICITY
FLAGS
SPREAD_MODEL
MAG_MODEL'''

####################### Auxiliary Functions #######################
def create_temp_files_for_SourceExtractor(f_nnw,f_conv,f_params,
                                          nnw_name,conv_name,params_name,
                                          config_name,sextractorloc='/opt/local/bin/source-extractor'):

    fp = open(nnw_name, "w")
    fp.write(f_nnw)
    fp.close()
    fp = open(conv_name, "w")
    fp.write(f_conv)
    fp.close()
    fp = open(params_name, "w")
    fp.write(f_params)
    fp.close()

    # Create temporary default.sex file
    command = sextractorloc+' -d > '+ config_name 
    subprocess.call(command,shell=True)

    return None

def remove_temp_files(fs):
    for f in fs:
        os.remove(f)
    return None

def run_sourceExtractor(fitsfiles,
                        sextractorloc='/opt/local/bin/source-extractor',
                        psfexloc='/opt/local/bin/psfex',
                        verbose=False,quietmode=False,debugmode=False):

    # Create temporary files required for Source Extractor
    # These files are defined as global variables above
    create_temp_files_for_SourceExtractor(f_nnw,f_conv,f_params,nnw_name,conv_name,params_name,
                                          config_name,sextractorloc=sextractorloc)

    # Print useful information if verbose
    if verbose:
        print('#####################################')
        print('###### RUNNING SOURCE EXTRACTOR #####')
        print('#####################################')

    # Set Source Extractor verbose type
    if verbose:
        VERBOSE_TYPE = 'NORMAL'
    else:
        VERBOSE_TYPE = 'QUIET'

    # Run Source Extractor one image at a time
    catfiles = []
    psffiles = []
    catted_fitsfiles = []
    for f in fitsfiles:

        # Print separator
        if verbose:
            print('==============================================')

        try:
            if verbose:
                print('Currently measuring PSF using source extractor and psfex...')
            [f_psf] = measure_psf([f], outdir='./', savepsffits=False,
                                    sextractorloc=sextractorloc,
                                    psfexloc=psfexloc,
                                    verbose=verbose,quietmode=True)
        except:
            print(f'\nSKIPPED: PSF measurement unsuccessful for {f}')
            continue

        # Run Source Extractor on image
        try:
            if verbose:
                print('Currently running source extractor using .psf input...')
            catalog_name = f.replace('.fits','_imageStats.cat')
            command =   f'{sextractorloc} -c {config_name} '\
                        f'-CATALOG_NAME {catalog_name} '\
                        f'-CATALOG_TYPE ASCII_HEAD '\
                        f'-PARAMETERS_NAME {params_name} -FILTER_NAME {conv_name} '\
                        f'-STARNNW_NAME {nnw_name} -MAG_ZEROPOINT 25.0 '\
                        f'-PSF_NAME {f_psf} -PSF_NMAX 1 -PATTERN_TYPE GAUSS-LAGUERRE '\
                        f'-VERBOSE_TYPE {VERBOSE_TYPE} '\
                        f'{f}'
            if verbose:
                print('Executing command: %s\n' % command)
            rval = subprocess.run(command.split(), check=True)
            catfiles.append(catalog_name)
            catted_fitsfiles.append(f)
            os.remove(f_psf)
            if verbose:
                print('Success!')
        except subprocess.CalledProcessError as err:
            print('\nCould not run SExtractor with exit error %s\n'%err)
            print('Command used:\n%s\n'%command)
    
    # Remove temporary files required for Source Extractor
    remove_temp_files([nnw_name,conv_name,params_name,config_name])
    
    return catfiles, catted_fitsfiles # end run_sourceExtractor

def get_imageStats(catfiles,fitsfiles,minfwhm=2.05,verbose=False,quietmode=False,debugmode=False):

    FWHMs = []
    FWHM_stds = []
    ELLIPs = []
    ELLIP_stds = []
    N_SRCs = []
    N_SRCs_stars = []
    BGRs = []
    BGR_stds = []

    # Print useful information if verbose
    if verbose:
        print('#####################################################')
        print('###### READING .cat FILES and CALCULATING STATS #####')
        print('#####################################################')

    for catfile, fitsfile in zip(catfiles,fitsfiles):
        # Print separator
        if verbose:
            print('==============================================')
            print(f'Working on: {catfile}')

        # Read in catalog
        dat = ascii.read(catfile) 
        df  = pd.DataFrame(dat.as_array())
        if verbose:
            print('\nFor file                     : ',catfile)
            print('Number of sources extracted  : ',len(df['CLASS_STAR']))
            print('Average CLASS_STAR           : %.5f'%np.average(df['CLASS_STAR']))
            print('Average ELLIPTICITY          : %.5f'%np.average(df['ELLIPTICITY']))
            print('Average FWHM_IMAGE           : %.5f'%np.average(df['FWHM_IMAGE']))

        # Pick catalog entries most likely to be stars
        # CLASS_STAR  1=star
        # Ellipticity 0=round
        df_stars = df[  (df['ELLIPTICITY']<0.3) & 
                        (df['CLASS_STAR']>0.8) &
                        (df['SPREAD_MODEL']<0.002) &
                        (df['FLAGS'] == 0) &
                        (df['FWHM_IMAGE'] > minfwhm)
                     ]

        # Calculate FWHMs,FWHM_stds,ELLIPs,ELLIP_stds,N_SRCs,N_SRCs_stars and save
        fwhm_avg = np.average(df_stars['FWHM_IMAGE'])
        fwhm_std = np.std(df_stars['FWHM_IMAGE'])
        FWHMs.append(fwhm_avg)
        FWHM_stds.append(fwhm_std)
        ellip_avg = np.average(df_stars['ELLIPTICITY'])
        ellip_std = np.std(df_stars['ELLIPTICITY'])
        ELLIPs.append(ellip_avg)
        ELLIP_stds.append(ellip_std)
        nsrcs = len(df['ELLIPTICITY'])
        nsrcs_stars = len(df_stars['ELLIPTICITY'])
        N_SRCs.append(nsrcs)
        N_SRCs_stars.append(nsrcs_stars)


        # Calculate BGRs,BGR_stds and save
        sky, _, sky_pixel_std = calculate_skyStats(fitsfile,  
                                                   sextractorloc=sextractorloc, 
                                                   verbose=verbose, debugmode=debugmode, quietmode=True)
        BGRs.append(sky)
        BGR_stds.append(sky_pixel_std)

        if verbose:
            print(f'The Average FWHM (pixels) is  : {fwhm_avg:.5f}')  
            print(f'The FWHM std (pixels) is      : {fwhm_std:.5f}')
            print(f'The Average ELLIP is          : {ellip_avg:.5f}')
            print(f'The ELLIP std is              : {ellip_std:.5f}')
            print(f'The number of detected srcs   : {nsrcs}')
            print(f'The number of detected stars  : {nsrcs_stars}')
            print(f'The BGR is                    : {sky}')
            print(f'The BGR std is                : {sky_pixel_std}')

    # Print useful information if verbose
    if verbose:
        print('\nCLASS_STAR  1=star\nEllipticity 0=round')

    return FWHMs,FWHM_stds,ELLIPs,ELLIP_stds,N_SRCs,N_SRCs_stars,BGRs,BGR_stds





def write_imageStats_to_file(outfile,catted_fitsfiles,
                            FWHMs,FWHM_stds,ELLIPs,ELLIP_stds,
                            N_SRCs,N_SRCs_stars,BGRs,BGR_stds,
                            quietmode=False):
    # Format numbers as desired, add on "heading" per column
    s_FWHMs             = [format(x,'.4f') for x in FWHMs]
    s_FWHMs.insert(0,'FWHM')
    s_FWHM_stds         = [format(x,'.4f') for x in FWHM_stds]
    s_FWHM_stds.insert(0,'FWHM_ERR')
    s_ELLIPs            = [format(x,'.4f') for x in ELLIPs]
    s_ELLIPs.insert(0,'ELLIP')
    s_ELLIP_stds        = [format(x,'.4f') for x in ELLIP_stds]
    s_ELLIP_stds.insert(0,'ELLIP_ERR')
    s_N_SRCs            = [format(x,'d') for x in N_SRCs]
    s_N_SRCs.insert(0,'N_SRCS')
    s_N_SRCs_stars      = [format(x,'d') for x in N_SRCs_stars]
    s_N_SRCs_stars.insert(0,'N_SRCS_STAR')
    s_BGRs              = [format(x,'.4f') for x in BGRs]
    s_BGRs.insert(0,'BGR')
    s_BGR_stds          = [format(x,'.4f') for x in BGR_stds]
    s_BGR_stds.insert(0,'BGR_ERR')
    catted_basenames    = [ntpath.basename(x) for x in catted_fitsfiles]
    catted_basenames.insert(0,'FILENAME')
    catted_fitsfiles.insert(0,'FULLPATH')
    # Transpose for numpy.savetxt
    savetext            = np.transpose([catted_basenames,
                                        s_FWHMs,s_FWHM_stds,
                                        s_ELLIPs,s_ELLIP_stds,
                                        s_N_SRCs,s_N_SRCs_stars,
                                        s_BGRs,s_BGR_stds,
                                        catted_fitsfiles])
    h                   = '0. Filename\n'\
                          '1. FWHM (pixels)\n'\
                          '2. FWHM standard deviation (pixels)\n'\
                          '3. Ellipticity, between 0 and 1, 0=round\n'\
                          '4. Ellipticity standard deivtation\n'\
                          '5. Number of sources detected\n'\
                          '6. Number of sources determined to be stars\n'\
                          '7. Background (median of all pixels outside of grown source extractor segmentation map) (ADU)\n'\
                          '8. Background standard deviation (ADU)\n'\
                          '9. Full path'
    # Save
    np.savetxt(outfile, (savetext),fmt='%s',header=h)
    # Print info it not quietmode
    if not quietmode:
        print(f'\nSaved: {outfile}')
        print(f'This file can be read into python with t=astropy.io.ascii.read({outfile})')
        print("t['KEYS'] gets columns where KEYS are: FILENAME FWHM FWHM_ERR ELLIP ELLIP_ERR N_SRCS N_SRCS_STAR BGR BGR_ERR FULLPATH")

    return None

####################### MAIN function #######################

def calculate_imageStats_spreadmodel(fitsfiles,outfile=False,
                                    sextractorloc='/opt/local/bin/source-extractor',
                                    psfexloc='/opt/local/bin/psfex',
                                    minfwhm=2.05,
                                    verbose=False,quietmode=False,debugmode=False):

    # Run source extractor on all input fits files
    if not quietmode:
        print('During the running of calculate_imageStats_spreadmodel, temporary .psf _imageStats.cat files will be created in the input file directory(s)')
        print('During the running of calculate_imageStats_spreadmodel, temporary files will also be created within the current directory.')
        print('If the code does not successfully finish, check to remove these temporary files.\n\n')
        print(f'Running source extractor, psfex, and source extractor again to get catalogs on {len(fitsfiles)} files...')
        print(f'It should take about 30-60s per file ({len(fitsfiles)/2:0.2f} - {len(fitsfiles):0.2f} mins in total) depending on your computer.')
        tick = time.perf_counter()
    catfiles, catted_fitsfiles = run_sourceExtractor(fitsfiles,
                                                     sextractorloc=sextractorloc,
                                                     psfexloc=psfexloc,
                                                     verbose=verbose,quietmode=quietmode,debugmode=debugmode)
    if not quietmode:
        tock = time.perf_counter()
        print(f'Getting catalogs took {(tock-tick)/60.:0.2f} minutes')
        
    if not quietmode:
        print('Reading catalogs to get image statistics, this ...')
        tick = time.perf_counter()
    # Calculate average FWHM, average ELLIPTICITY, # sources, background of each input fits file
    (FWHMs,FWHM_stds,
    ELLIPs,ELLIP_stds,
    N_SRCs,N_SRCs_stars,
    BGRs,BGR_stds)          = get_imageStats(catfiles,catted_fitsfiles,
                                             minfwhm=minfwhm,
                                             verbose=verbose,quietmode=quietmode,debugmode=debugmode)
    if not quietmode:
        tock = time.perf_counter()
        print(f'Reading catalogs and calculating image stats took {(tock-tick)/60.:0.2f} minutes')

    # Save if required to
    if outfile:
        write_imageStats_to_file(outfile,catted_fitsfiles,
                                 FWHMs,FWHM_stds,ELLIPs,ELLIP_stds,
                                 N_SRCs,N_SRCs_stars,BGRs,BGR_stds,
                                 quietmode=quietmode)

    # Remove catalog files if not in debug mode
    if not debugmode:
        remove_temp_files(catfiles)

    return catted_fitsfiles,FWHMs,FWHM_stds,ELLIPs,ELLIP_stds,N_SRCs,N_SRCs_stars,BGRs,BGR_stds

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
    sextractorloc   = arguments['--sextractor']
    psfexloc        = arguments['--psfex']
    minfwhm         = float(arguments['--minfwhm']) 
    outfile         = arguments['--out']

    # Calculate
    _ = calculate_imageStats_spreadmodel(fitsfiles,outfile=outfile,
                                        sextractorloc=sextractorloc,
                                        psfexloc=psfexloc,
                                        minfwhm=minfwhm,
                                        verbose=verbose,quietmode=quietmode,debugmode=debugmode)
