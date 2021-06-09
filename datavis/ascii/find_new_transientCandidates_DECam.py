#!/usr/bin/python

#!/usr/bin/env python

""" 
find_new_transientCandidates_DECam.py -- Find transient candidates in subtraction image catalogs. Transient candidates turn up in Science AND Subtraction Image, but not in Neg of Subtraction image and NOT in template image. 

Inputs 
    (1) subtraction image source extractor ASCII catalogs (without SPREAD_MODEL)
    (2) corresponding science source extractor ASCII catalogs (without SPREAD_MODEL)
    (3) corresponding -1*subtraction source extractor ASCII catalogs (without SPREAD_MODEL)
    {4} corresponding template image source extractor ASCII catalog (without SPREAD_MODEL)
    (4) Subtraction Image path    

Outputs three ds9 region files
    (1) Green circles: sources in subtraction and science image
    (2) Yellow pandas: sources in subtraction and science image, but not in inverted sub image
    (2) Magenta pandas: sources in subtraction and science image, but not in inverted sub image; Also not in template image

Caveats: 
!!! Currently parameters are optimised for DECam NOAO Community Pipeline Stacked images 

Usage: find_new_transientCandidates_DECam [-h] [-v] [--debug] [-q] [-o SAVEDIR] <sub_cat> <sci_cat> <neg_cat> <templ_cat>

Options:
    -h, --help                              Show this screen
    -v, --verbose                           Show extra information [default: False]   
    -q, --quietmode                         Minimize print to screen. This is useful when this function is called in another function. [default: False]  
    --debug                                 Output more for debugging [default: False]
    -o SAVEDIR, --out SAVEDIR               Saved output as [default: ./]

"""
import docopt
import astropy.io.fits as fits
import numpy as np
import sys, os
from astropy.io import ascii
import pandas as pd
from astropy.coordinates import SkyCoord
from astropy import units as u
import numpy as np
from pathlib import Path
import ntpath
import subprocess

from datavis.ascii.match_catalogs import match_catalogs_df

__author__      = "Jielai Zhang"
__license__     = "MIT"
__version__     = "1.0.1"
__date__        = "2021-06-06"
__maintainer__  = "Jielai Zhang"
__email__       = "zhang.jielai@gmail.com"

##############################################################
####################### Main Function ########################
##############################################################

def find_new_transientCandidates_DECam(f_sub_cat, f_sci_cat, f_neg_cat, f_templ_cat,
                                        savedir='./',
                                        verbose=False,debugmode=False,quietmode=False):

    if verbose:
        print('*'*10)
        print('find_new_transientCandidates_DECam.py running...')
        print('*'*10)

    # ..............................................................
    # If savedir doesn't exist, make it
    if not os.path.isdir(savedir):
        os.makedirs(savedir)

    # ..............................................................
    # Determine output file names
    reg_sciYes_subYes = f_sci_cat.replace('.cat','_sciYes_subYes.reg')
    reg_sciYes_subYes_negNo = f_sci_cat.replace('.cat','_sciYes_subYes_negNo.reg')
    reg_sciYes_subYes_negNo_templNo = f_sci_cat.replace('.cat','_sciYes_subYes_negNo_templNo.reg')

    # ..............................................................
    # Read in catalogs.
    dat_pos   = ascii.read(f_sub_cat)
    dat_neg   = ascii.read(f_neg_cat)
    dat_sci   = ascii.read(f_sci_cat)
    dat_templ = ascii.read(f_templ_cat)

    # ..............................................................
    # Turn it into pandas data frame. 
    df_neg_select    = pd.DataFrame(dat_neg.as_array())
    df_sci_select    = pd.DataFrame(dat_sci.as_array())
    df_pos_select    = pd.DataFrame(dat_pos.as_array())
    df_templ_select  = pd.DataFrame(dat_templ.as_array())

    if verbose:
        print(f'For {f_sci_cat}:')
        print(f'Sources in TEMPL: {len(df_templ_select)}')
        print(f'Sources in   SCI: {len(df_sci_select)}')
        print(f'Sources in   POS: {len(df_pos_select)}')
        print(f'Sources in   NEG: {len(df_neg_select)}')  

    # ..............................................................
    # Set candidates as a source in Science with corresponding 
    # source in Subtraction, but no nearby source in Negative.
    # XXX Hard coded section !!! (radius_threshold)

    # Track if there are candidates for each round before proceeding. 
    no_candidates = False
    no_sciYes_posYes_candidates = False
    no_sciYes_posYes_negNo_candidates = False
    no_sciYes_posYes_negNo_templNo_candidates = False

    # Round 1: some positive detections made
    if len(df_pos_select) == 0:
        print('No subtraction sources were detected.')
        no_candidates                     = True
        no_sciYes_posYes_candidates       = True
        no_sciYes_posYes_negNo_candidates = True
    else:
        # Get sci and sub matched
        radius_threshold = 1.1 * u.arcsec # Seeing ~0.8-1.2 at DECam
        df_sciYes_posYes, df_sciYes_posYes_sub = match_catalogs_df( df_sci_select, df_pos_select,
                                                                    radius_threshold=radius_threshold)

    # Round 2: for each science detection there is a positive detection associated
    if no_candidates == False:
        if len(df_sciYes_posYes) == 0:
            print('No subtraction sources were found to have a match in the science source list.')
            no_sciYes_posYes_candidates       = True
            no_sciYes_posYes_negNo_candidates = True
            no_sciYes_posYes_negNo_templNo_candidates = True
        else:        
            # Get sci/sub matched with no nearby neg
            radius_threshold = 4 * u.arcsec
            df_sciYes_posYes_negNo, df_sciYes_posYes_negNo_neg = match_catalogs_df( df_sciYes_posYes, df_neg_select,
                                                                                    radius_threshold=radius_threshold,find_close=False)

    # Round 3: for each sci/sub pair, there is no nearby negative detection
    if no_sciYes_posYes_candidates == False:
        if len(df_sciYes_posYes_negNo) == 0:
            print('No positive sources with matching subtraction sources with no nearby negative sources has been found.')
            no_sciYes_posYes_negNo_candidates = True    
            no_sciYes_posYes_negNo_templNo_candidates = True
        else:
            # Get sci/sub/no nearby neg with no templ
            radius_threshold = 1.1 * u.arcsec # Seeing ~0.8-1.2 at DECam
            df_sciYes_posYes_negNo_templNo, df_sciYes_posYes_negNo_neg_templNo_templ = match_catalogs_df( df_sciYes_posYes_negNo, df_templ_select,
                                                                                    radius_threshold=radius_threshold,find_close=False)
            
    # Round 4: for each sci/sub/no neg source, there is no associated source in template
    if no_sciYes_posYes_negNo_candidates == False:
        if len(df_sciYes_posYes_negNo_templNo) == 0:
            print('No pos sources with matching sub srcs with no nearby neg sources and no associated template source has been found.')
            no_sciYes_posYes_negNo_templNo_candidates = True
            
    if verbose:
        if no_sciYes_posYes_candidates == False:
            print(f'Sources in SCI matched in POS SUB: {len(df_sciYes_posYes)}')
        if no_sciYes_posYes_negNo_candidates == False:
            print(f'Sources in SCI matched in POS SUB, then not matched in NEG: {len(df_sciYes_posYes_negNo)}')
        if no_sciYes_posYes_negNo_candidates == False:
            print(f'Sources in SCI matched in POS SUB, then not matched in NEG, then not matched in template: {len(df_sciYes_posYes_negNo_templNo)}')

    # ..............................................................
    # Make Region files of sources matched in Sci and Sub

    if no_sciYes_posYes_candidates == False:
        # Save region file.
        Rs_sciYes_posYes = np.array(df_sciYes_posYes['X_WORLD'])
        Ds_sciYes_posYes = np.array(df_sciYes_posYes['Y_WORLD'])
        Sizes   = ["0.1'"]*len(Rs_sciYes_posYes)
        Shapes = ['circle']*len(Rs_sciYes_posYes)
        j2000s  = ['j2000;']*len(Rs_sciYes_posYes)
        Colors = ['#color=green']*len(Rs_sciYes_posYes)
        h       = 'Region file format: DS9 version 4.0'
        savetext = np.transpose([j2000s,Shapes,Rs_sciYes_posYes,Ds_sciYes_posYes,Sizes])
        np.savetxt(reg_sciYes_subYes, (savetext),fmt='%s',header=h)
        if verbose:
            print(f'\nSaved: {reg_sciYes_subYes}')
    else:
        # Save empty region file.
        h       = 'Region file format: DS9 version 4.0; Yes, this file is empty of regions, however, it should not happen!.'
        f_open= open(reg_sciYes_subYes,"w+")
        f_open.write(h)
        f_open.close()
        if verbose:
            print(f'Saved: {reg_sciYes_subYes}')


    # ..............................................................
    # Make Region files of sources matched in Sci and Sub and not in Neg

    if no_sciYes_posYes_negNo_candidates == False:
        # Save region file.
        Rs_sciYes_posYes_negNo = np.array(df_sciYes_posYes_negNo['X_WORLD'])
        Ds_sciYes_posYes_negNo = np.array(df_sciYes_posYes_negNo['Y_WORLD'])
        Shapes  = ['panda']*len(Rs_sciYes_posYes_negNo)
        j2000s  = ['j2000;']*len(Rs_sciYes_posYes_negNo)
        Colors  = ['#color=yellow']*len(Rs_sciYes_posYes_negNo)
        startangles=['0']*len(Rs_sciYes_posYes_negNo)
        stopangles =['360']*len(Rs_sciYes_posYes_negNo)
        nangles    =['4']*len(Rs_sciYes_posYes_negNo)
        inners     =['0.002']*len(Rs_sciYes_posYes_negNo)
        outers     =['0.0025']*len(Rs_sciYes_posYes_negNo)
        nradiuss   =['1']*len(Rs_sciYes_posYes_negNo)
        h       = 'Region file format: DS9 version 4.0'
        savetext = np.transpose([j2000s,Shapes,Rs_sciYes_posYes_negNo,Ds_sciYes_posYes_negNo,
                                 startangles,stopangles,nangles,inners,outers,nradiuss,Colors])
        np.savetxt(reg_sciYes_subYes_negNo, (savetext),fmt='%s',header=h)
        if verbose:
            print(f'Saved: {reg_sciYes_subYes_negNo}')
    else:
        # Save empty region file.
        h       = 'Region file format: DS9 version 4.0; Yes, this file is empty of regions, however, it is a bit surprising.'
        f_open= open(reg_sciYes_subYes_negNo,"w+")
        f_open.write(h)
        f_open.close()
        if verbose:
            print(f'Saved: {reg_sciYes_subYes_negNo}')


    # ..............................................................
    # Make Region files of sources matched in Sci and Sub and not in Neg

    if no_sciYes_posYes_negNo_templNo_candidates == False:
        # Save region file.
        Rs_sciYes_posYes_negNo_templNo = np.array(df_sciYes_posYes_negNo_templNo['X_WORLD'])
        Ds_sciYes_posYes_negNo_templNo = np.array(df_sciYes_posYes_negNo_templNo['Y_WORLD'])
        Shapes  = ['panda']*len(Rs_sciYes_posYes_negNo_templNo)
        j2000s  = ['j2000;']*len(Rs_sciYes_posYes_negNo_templNo)
        Colors  = ['#color=magenta']*len(Rs_sciYes_posYes_negNo_templNo)
        startangles=['0']*len(Rs_sciYes_posYes_negNo_templNo)
        stopangles =['360']*len(Rs_sciYes_posYes_negNo_templNo)
        nangles    =['4']*len(Rs_sciYes_posYes_negNo_templNo)
        inners     =['0.003']*len(Rs_sciYes_posYes_negNo_templNo)
        outers     =['0.0035']*len(Rs_sciYes_posYes_negNo_templNo)
        nradiuss   =['1']*len(Rs_sciYes_posYes_negNo_templNo)
        h       = 'Region file format: DS9 version 4.0'
        savetext = np.transpose([j2000s,Shapes,Rs_sciYes_posYes_negNo_templNo,Ds_sciYes_posYes_negNo_templNo,
                                 startangles,stopangles,nangles,inners,outers,nradiuss,Colors])
        np.savetxt(reg_sciYes_subYes_negNo_templNo, (savetext),fmt='%s',header=h)
        if verbose:
            print(f'Saved: {reg_sciYes_subYes_negNo_templNo}')
    else:
        # Save empty region file.
        h       = 'Region file format: DS9 version 4.0; Yes, this file is empty of regions.'
        f_open= open(reg_sciYes_subYes_negNo_templNo,"w+")
        f_open.write(h)
        f_open.close()
        if verbose:
            print(f'Saved: {reg_sciYes_subYes_negNo_templNo}')

        

    # ..............................................................
    # Helpful ds9 command hint:
    if verbose:
        print('\nHelpful ds9 command hint:')
        print(f'ds9 -zscale -lock frame wcs')

    # ..............................................................
    # Return appropriate outputs    
    return reg_sciYes_subYes, reg_sciYes_subYes_negNo, reg_sciYes_subYes_negNo_templNo

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
    f_sub_cat         = arguments['<sub_cat>']
    f_sci_cat         = arguments['<sci_cat>']
    f_neg_cat         = arguments['<neg_cat>']   
    f_templ_cat       = arguments['<templ_cat>']

    _,_,_ = find_new_transientCandidates_DECam(f_sub_cat, f_sci_cat, f_neg_cat, f_templ_cat, 
                                                savedir=savedir,
                                                verbose=verbose,debugmode=debugmode,quietmode=quietmode)
