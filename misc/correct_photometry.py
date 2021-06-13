#!/usr/bin/env python

""" correct_photometry.py -- Enter reference catalogue and catalogue to be photometry corrected. Output photometry corrected catalogue into output directory with name input_sci_cat_CORRECTED.cat. 

Usage: correct_photometry.py [-h] [-v] [--debug] [-q] [-o DIRECTORY] [--REF_RA_KEY STRING] [--REF_DEC_KEY STRING] [--REF_PHOTOM_KEY STRING] [--SCI_RA_KEY STRING] [--SCI_DEC_KEY STRING] [--SCI_PHOTOM_KEY STRING] [--RADIUS_THRESHOLD FLOAT] <reference_cat> <science_cat>

Options:
    -h, --help                              Show this screen
    -v, --verbose                           Show extra information [default: False]   
    -q, --quietmode                         Minimize print to screen. This is useful when this function is called in another function. [default: False]  
    --debug                                 Output more for debugging [default: False]
    -o DIRECTORY, --out DIRECTORY           Saved _CORRECTED.cat here. If not specified, same in same directory as input science catalogue.
    --REF_RA_KEY STRING                     Ref RA column name [default: RA]
    --REF_DEC_KEY STRING                    Ref DEC column name [default: DEC]
    --REF_PHOTOM_KEY STRING                 Ref magnitude column name [default: gmag]
    --SCI_RA_KEY STRING                     Sci RA column name [default: X_WORLD]
    --SCI_DEC_KEY STRING                    Sci DEC column name [default: Y_WORLD]
    --SCI_PHOTOM_KEY STRING                 Sci magnitude column name [default: MAG_AUTO]
    --RADIUS_THRESHOLD FLOAT                How close does ref source have to be to science source to be a match, arcsec [default: 1.0]

Examples:
"""
import docopt
import pandas as pd
import astropy.io.ascii as ascii
from astropy import units as u
import numpy as np

from datavis.ascii.match_catalogs import match_catalogs_df

__author__      = "Jielai Zhang"
__license__     = "MIT"
__version__     = "1.0.1"
__date__        = "2021-06-13"
__maintainer__  = "Jielai Zhang"
__email__       = "zhang.jielai@gmail.com"

#########################################
# ======= House Keeping Functions =======
#########################################

##############################################################
####################### Main Function ########################
##############################################################

def correct_photometry(f_ref_cat,f_sci_cat,savedir=None,
                            radius_threshold=1.0*u.arcsec,
                            REF_RA_KEY='RA',REF_DEC_KEY='DEC',REF_PHOTOM_KEY='gmag',
                            SCI_RA_KEY='X_WORLD',SCI_DEC_KEY='Y_WORLD',SCI_PHOTOM_KEY='MAG_AUTO',
                            verbose=False,debugmode=False,quietmode=False):

    # Read in Catalogues
    cat_ref = ascii.read(f_ref_cat)
    cat_sci = ascii.read(f_sci_cat)
    df_ref  = pd.DataFrame(cat_ref.as_array())
    df_sci  = pd.DataFrame(cat_sci.as_array())

    # Match catalogues by RA and DEC
    df_ref_matched, df_sci_matched = match_catalogs_df(df_ref, df_sci,radius_threshold=radius_threshold,
                                               ref_RA_KEY=REF_RA_KEY,ref_DEC_KEY=REF_DEC_KEY,
                                               sci_RA_KEY=SCI_RA_KEY,sci_DEC_KEY=SCI_DEC_KEY)
    if debugmode:
        print(f'DEBUG: The number of ref+sci matched sources = {len(df_ref_matched)}') 

    # Calculate ZP
    try:
        ref_mags = np.array([float(x) for x in df_ref_matched[REF_PHOTOM_KEY]])
    except:
        ref_mags_temp = np.array(df_ref_matched[REF_PHOTOM_KEY])
        for ii,v in enumerate(ref_mags_temp):
            if v == 'NA':
                ref_mags_temp[ii] = np.nan
        ref_mags = np.array([float(x) for x in ref_mags_temp])
    sci_mags = np.array([float(x) for x in df_sci_matched[SCI_PHOTOM_KEY]])
    zps      = ref_mags - sci_mags
    if debugmode:
        print(f'DEBUG: The average {SCI_PHOTOM_KEY} of matched sci sources is {np.nanmean(sci_mags)}')
        print(f'DEBUG: The average {REF_PHOTOM_KEY} of matched ref sources is {np.nanmean(ref_mags)}')
        print(f'DEBUG: The average of the difference is: {np.nanmean(zps)}')
        print(f'DEBUG: The median of the difference is: {np.nanmedian(zps)}')
        print('DEBUG: Use the median difference to correct photometry.')
    ZP = np.nanmedian(zps)

    # Calculate corrected photometry
    new_sci_mag = df_sci[SCI_PHOTOM_KEY]+ZP
    if debugmode:
        print(f'DEBUG: The median science magnitude is: {np.median(new_sci_mag)}')

    # Save new catalogue
    file_ext = f_sci_cat.split('.')[-1]
    if savedir==None:
        f_save = f_sci_cat.replace('.'+file_ext,'_CORRECTED.'+file_ext)
    else:   
        fname = ntpath.basename(f_sci_cat)
        f_save = savedir+os.sep+fname.replace('.'+file_ext,'_CORRECTED.'+file_ext)
    df_sci[SCI_PHOTOM_KEY] = new_sci_mag
    header_string = ''
    for k in np.array(df_sci.keys()):
        header_string=header_string+k+' '
    np.savetxt(f_save, df_sci.values, fmt='%f',header=header_string)
    if not quietmode:
        print(f'Saved: {f_save}')
    
    return f_save

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
    savedir         = arguments['--out']
    f_ref_cat       = arguments['<reference_cat>']
    f_sci_cat       = arguments['<science_cat>']

    radius_threshold    = float(arguments['--RADIUS_THRESHOLD'])*u.arcsec
    REF_RA_KEY          = arguments['--REF_RA_KEY']
    REF_DEC_KEY         = arguments['--REF_DEC_KEY']
    REF_PHOTOM_KEY      = arguments['--REF_PHOTOM_KEY']
    SCI_RA_KEY          = arguments['--SCI_RA_KEY']
    SCI_DEC_KEY         = arguments['--SCI_DEC_KEY']
    SCI_PHOTOM_KEY      = arguments['--SCI_PHOTOM_KEY']

    _ = correct_photometry(f_ref_cat,f_sci_cat,savedir=savedir,
                            radius_threshold=radius_threshold,
                            REF_RA_KEY=REF_RA_KEY,REF_DEC_KEY=REF_DEC_KEY,REF_PHOTOM_KEY=REF_PHOTOM_KEY,
                            SCI_RA_KEY=SCI_RA_KEY,SCI_DEC_KEY=SCI_DEC_KEY,SCI_PHOTOM_KEY=SCI_PHOTOM_KEY,
                            verbose=verbose,debugmode=debugmode,quietmode=quietmode)
