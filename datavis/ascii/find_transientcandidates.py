#!/usr/bin/python

#!/usr/bin/env python

""" find_transientCandidates.py -- Given three input text files each containing a list of source extractor ASCII catalogs for (1) subtraction images (2) corresponding science images and (3) corresponding -1*subtraction images, output 
(1) a list of candidate transients as text files for each input subtraction image to a specified folder. This file also contains a ds9 command to look at the candidates in science and subtraction. 
(2) a region file for the found transient candidates to the same specified folder. 

Note: fits files are input here, but they are not read in at all. The path is only used to write ds9 commands. 

Caveats: 
!!! Currently source ID parameters are optimised for KMTNet data for DWF 
!!! The separation distance between sci/sub and sci/neg is set optimally for KTMNet right now. It is HARD CODED in! 

Usage: find_transientCandidates [-h] [-v] [--debug] [-q] [-o SAVEDIR] [--RAmin FLOAT] [--RAmax FLOAT] [--DECmin FLOAT] [--DECmax FLOAT] [--list] <sub_cat> <sci_cat> <neg_cat> <sub_fits> <sci_fits> <neg_fits>

Options:
    -h, --help                              Show this screen
    -v, --verbose                           Show extra information [default: False]   
    -q, --quietmode                         Minimize print to screen. This is useful when this function is called in another function. [default: False]  
    --debug                                 Output more for debugging [default: False]
    -o SAVEDIR, --out SAVEDIR               Saved output as [default: ./]
    --RAmax FLOAT                           If you don't want a list of cands in whole catalog, but just within RA DEC ranges. Input in degrees.
    --RAmin FLOAT                           If you don't want a list of cands in whole catalog, but just within RA DEC ranges. Input in degrees.
    --DECmax FLOAT                          If you don't want a list of cands in whole catalog, but just within RA DEC ranges. Input in degrees.
    --DECmin FLOAT                          If you don't want a list of cands in whole catalog, but just within RA DEC ranges. Input in degrees.
    --list                                  If specified, the input sub_cat, sci_cat, neg_cat should be lists [default: False]

Examples:
from datavis.ascii.find_transientCandidates import find_transientCandidates
     (Region_Files, 
        Candidate_List_Files, 
        Candidate_Lists) = find_transientCandidates(sub_cat, sci_cat, neg_cat, input_list=False, savedir='./',
                                                    verbose=False,debugmode=False,quietmode=False):
"""
import docopt
import astropy.io.fits as fits
import numpy as np
import sys, os
from astropy.io import ascii
import pandas as pd
from astropy.coordinates import SkyCoord
from astropy import units as u
import corner
import numpy as np
from pathlib import Path
import ntpath
import subprocess

from datavis.ascii.match_catalogs import match_catalogs_df

__author__      = "Jielai Zhang"
__license__     = "MIT"
__version__     = "1.0.1"
__date__        = "2020-11-03"
__maintainer__  = "Jielai Zhang"
__email__       = "zhang.jielai@gmail.com"

#################################
# ======= Other Functions =======
#################################
def find_transientCandidates_(sub_cat, sci_cat, neg_cat, savedir='./',
                              RAmax=None,RAmin=None,DECmax=None,DECmin=None,
                              verbose=False,debugmode=False,quietmode=False):

    # ..............................................................
    # Determine where reg file and transient cat is to be saved.
    fname                   = ntpath.basename(sub_cat)
    fname_stub              = Path(fname).stem
    candidate_list_file     = savedir + os.path.sep + fname_stub + '.list'
    ds9commands_file        = savedir + os.path.sep + fname_stub + '_ds9.sh'
    # This is sciYes_subYes_negNo region file
    region_file             = savedir + os.path.sep + fname_stub + '.reg' 
    # These are deleted after code is done.
    reg_pos_temp            = savedir + os.path.sep + '_pos_temp.reg'
    reg_sci_temp            = savedir + os.path.sep + '_sci_temp.reg'
    reg_neg_temp            = savedir + os.path.sep + '_neg_temp.reg'
    # These are not deleted after code is done, but their paths are not returned. 
    reg_all_selected_detections = savedir + os.path.sep + fname_stub + '_selected_detections.reg'
    reg_sciYes_subYes           = savedir + os.path.sep + fname_stub + '_sciYes_subYes.reg'
    
    # ..............................................................
    # Check that RA and DEC limits are all in, and not just partial.
    if RAmax or RAmin or DECmax or DECmin:
        if RAmax == None:
            sys.exit('ERROR: Please input RAmax if you input some RA DEC constraints')
        if RAmin == None:
            sys.exit('ERROR: Please input RAmin if you input some RA DEC constraints')
        if DECmax == None:
            sys.exit('ERROR: Please input DECmax if you input some RA DEC constraints')
        if DECmin == None:
            sys.exit('ERROR: Please input DECmin if you input some RA DEC constraints')

    # ..............................................................
    # Read in catalogs.
    dat_neg = ascii.read(neg_cat) 
    dat_sci = ascii.read(sci_cat)
    dat_pos = ascii.read(sub_cat) 

    # ..............................................................
    # Turn it into pandas data frame. 
    df_neg  = pd.DataFrame(dat_neg.as_array())
    df_sci  = pd.DataFrame(dat_sci.as_array())
    df_pos  = pd.DataFrame(dat_pos.as_array())

    if verbose:
        print(f'For {sub_cat}:')
        print(f'Sources in SCI: {len(df_sci)}')
        print(f'Sources in POS: {len(df_pos)}')
        print(f'Sources in NEG: {len(df_neg)}')

    # ..............................................................
    # Select sources that's within RA DEC Range
    if RAmax:
        df_sci_RADECselect = df_sci[(   (df_sci['X_WORLD']>RAmin) &
                                        (df_sci['X_WORLD']<RAmax) &
                                        (df_sci['Y_WORLD']>DECmin) &
                                        (df_sci['Y_WORLD']<DECmax) 
                                    )
                                    ]
        df_pos_RADECselect = df_pos[(   (df_pos['X_WORLD']>RAmin) &
                                        (df_pos['X_WORLD']<RAmax) &
                                        (df_pos['Y_WORLD']>DECmin) &
                                        (df_pos['Y_WORLD']<DECmax) 
                                    )
                                    ]
        df_neg_RADECselect = df_neg[(   (df_neg['X_WORLD']>RAmin) &
                                        (df_neg['X_WORLD']<RAmax) &
                                        (df_neg['Y_WORLD']>DECmin) &
                                        (df_neg['Y_WORLD']<DECmax) 
                                    )
                                    ]
        if verbose:
            print(f'Sources within specified RA, DEC range in SCI: {len(df_sci_RADECselect)}')
            print(f'Sources within specified RA, DEC range in POS: {len(df_pos_RADECselect)}')
            print(f'Sources within specified RA, DEC range in NEG: {len(df_neg_RADECselect)}')

    else:
        df_sci_RADECselect = df_sci
        df_pos_RADECselect = df_pos
        df_neg_RADECselect = df_neg

    # ..............................................................
    # Make cuts on sources detected
    df_sci_select = df_sci_RADECselect[(  (df_sci_RADECselect['FWHM_IMAGE']>3.0) &
                                          (df_sci_RADECselect['ELLIPTICITY']<0.77) &
                                          (df_sci_RADECselect['SPREAD_MODEL']>-0.025) &
                                          (df_sci_RADECselect['SPREAD_MODEL']<0.05) &
                                          (df_sci_RADECselect['FLAGS']==0) 
                                        )
                                       ]


    df_pos_select = df_pos_RADECselect[(  (df_pos_RADECselect['FWHM_IMAGE']<10) &
                                          (df_pos_RADECselect['FWHM_IMAGE']>2.0) &
                                          (df_pos_RADECselect['MAG_AUTO']!=99) &
                                          (df_pos_RADECselect['SPREAD_MODEL']>-0.025) &
                                          (df_pos_RADECselect['SPREAD_MODEL']<0.1) 
                                        )
                                       ]

    df_neg_select = df_neg_RADECselect[(  (df_neg_RADECselect['FWHM_IMAGE']<35) &
                                          (df_neg_RADECselect['FWHM_IMAGE']>2.0) &
                                          (df_neg_RADECselect['MAG_AUTO']!=99) 
                                        )
                                       ]
    if verbose:
        print(f'Sources specified RA, DEC range + satisfy some conditions in SCI: {len(df_sci_select)}')
        print(f'Sources specified RA, DEC range + satisfy some conditions in POS: {len(df_pos_select)}')
        print(f'Sources specified RA, DEC range + satisfy some conditions in NEG: {len(df_neg_select)}')    

    # ..............................................................
    # Set candidates as a source in Science with corresponding 
    # source in Subtraction, but no nearby source in Negative.

    radius_threshold = 2.5 * u.arcsec # 10 pixels for KMTNet
    df_sciYes_posYes, df_sciYes_posYes_sub = match_catalogs_df( df_sci_select, df_pos_select,
                                                                radius_threshold=radius_threshold)
    radius_threshold = 4 * u.arcsec
    df_sciYes_posYes_negNo, df_sciYes_posYes_negNo_neg = match_catalogs_df( df_sciYes_posYes, df_neg_select19,
                                                                            radius_threshold=radius_threshold,find_close=False)
    # Grab the corresponding sub too (this line doesn't change df_sciYes_posYes_negNo)
    df_sciYes_posYes_negNo, df_sciYes_posYes_negNo_sub = match_catalogs_df( df_sciYes_posYes_negNo, df_pos_select19,
                                                                            radius_threshold=radius_threshold,find_close=True)

    if verbose:
        print(f'Sources in SCI matched in POS SUB: {len(Rs_sciYes_posYes)}')
        print(f'Sources in SCI matched in POS SUB, then not matched in NEG: {len(Ds_sciYes_posYes_negNo)}')


    # ..............................................................
    # Open ds9commands file for overwriting
    f_ds9 = open(ds9commands_file, "w")
    f_ds9.write("# ====================================================")
    f_ds9.write(f'# {sub_cat}')

    # ..............................................................
    # Make Region files of SELECTED sources in Sci, Sub, Neg

    # Sub Detections Selected (Green squares)
    Rs = df_pos_select19['X_WORLD']
    Ds = df_pos_select19['Y_WORLD']
    Widths = ["0.09'"]*len(Rs)
    Heights = ["0.00'"]*len(Rs)
    Angles = ["0.0"]*len(Rs)
    Shapes = ['box']*len(Rs)
    j2000s = ['j2000;']*len(Rs)
    Colors = ['#color=green']*len(Rs)
    h = 'Region file format: DS9 version 4.0'
    savetext = np.transpose([j2000s,Shapes,Rs,Ds,Widths,Heights,Angles,Colors])
    np.savetxt(reg_pos_temp, (savetext),fmt='%s',header=h)

    # Neg Detections Selected (Red squares, bit bigger than Sub Squares)
    Rs = df_neg_select['X_WORLD']
    Ds = df_neg_select['Y_WORLD']
    Widths = ["0.11'"]*len(Rs)
    Heights = ["0.11'"]*len(Rs)
    Angles = ["0.0"]*len(Rs)
    Shapes = ['box']*len(Rs)
    j2000s = ['j2000;']*len(Rs)
    Colors = ['#color=red']*len(Rs)
    h = 'Region file format: DS9 version 4.0'
    savetext = np.transpose([j2000s,Shapes,Rs,Ds,Widths,Heights,Angles,Colors])
    np.savetxt(reg_neg_temp, (savetext),fmt='%s',header=h)

    # Sci Detection Selected (Blue squares, big bigger than Sub Squares, at 30 deg angle)
    Rs = df_sci_select19['X_WORLD']
    Ds = df_sci_select19['Y_WORLD']
    Widths = ["0.11'"]*len(Rs)
    Heights = ["0.11'"]*len(Rs)
    Angles = ["30.0"]*len(Rs)
    Shapes = ['box']*len(Rs)
    j2000s = ['j2000;']*len(Rs)
    Colors = ['#color=blue']*len(Rs)
    h = 'Region file format: DS9 version 4.0'
    savetext = np.transpose([j2000s,Shapes,Rs,Ds,Widths,Heights,Angles,Colors])
    np.savetxt(reg_sci_temp, (savetext),fmt='%s',header=h)

    # Concatenate to region file.
    command = f'tail -n +2 {reg_pos_temp} > {reg_all_selected_detections}'
    subprocess.call(command,shell=True)
    command = f'tail -n +2 {reg_neg_temp} >> {reg_all_selected_detections}'
    subprocess.call(command,shell=True)
    command = f'tail -n +2 {reg_sci_temp} >> {reg_all_selected_detections}'
    subprocess.call(command,shell=True)

    # Add ds9 command to ds9 command file
    f_ds9.write('# ---Sub, neg sub, sci, detected sources region file --- ')
    command =   f'ds9 -zscale -lock frame wcs {f_pos} {f_neg} {f_sci} '\
                f'-regions load all {regionfile_detections} &' 
print(command)

print('')
print('green = sub')
print('red = negative sub')
print('blue = sci')

    # This is sciYes_subYes_negNo region file
    region_file             = savedir + os.path.sep + fname_stub + '.reg' 
    # These are deleted after code is done.
    reg_pos_temp            = savedir + os.path.sep + '_pos_temp.reg'
    reg_sci_temp            = savedir + os.path.sep + '_sci_temp.reg'
    reg_neg_temp            = savedir + os.path.sep + '_neg_temp.reg'
    # These are not deleted after code is done, but their paths are not returned. 
    reg_all_selected_detections = savedir + os.path.sep + fname_stub + '_selected_detections.reg'
    reg_sciYes_subYes           = savedir + os.path.sep + fname_stub + '_sciYes_subYes.reg'


    # ..............................................................
    # Make Region files of sources matched in Sci and Sub
    Rs_sciYes_posYes = df_sciYes_posYes['X_WORLD']
    Ds_sciYes_posYes = df_sciYes_posYes['Y_WORLD']


    # ..............................................................
    # Make Region files of sources matched in Sci and Sub and not in Neg
    Rs_sciYes_posYes_negNo = np.array(df_sciYes_posYes_negNo['X_WORLD'])
    Ds_sciYes_posYes_negNo = np.array(df_sciYes_posYes_negNo['Y_WORLD'])

    
    # ..............................................................
    # Close ds9commands file 
    f_ds9.close()  

    return region_file, candidate_list_file, candidate_list, ds9commands_file

##############################################################
####################### Main Function ########################
##############################################################

def find_transientCandidates(sub_cat, sci_cat, neg_cat, input_list=False, savedir='./',
                            RAmax=None,RAmin=None,DECmax=None,DECmin=None,
                            verbose=False,debugmode=False,quietmode=False):

    Region_Files         = []
    Candidate_List_Files = []
    Candidate_Lists      = []
    ds9commands_Files    = []

    # Check that RA and DEC limits are all in, and not just partial
    if RAmax or RAmin or DECmax or DECmin:
        if RAmax == None:
            sys.exit('ERROR: Please input RAmax if you input some RA DEC constraints.')
        if RAmin == None:
            sys.exit('ERROR: Please input RAmin if you input some RA DEC constraints.')
        if DECmax == None:
            sys.exit('ERROR: Please input DECmax if you input some RA DEC constraints.')
        if DECmin == None:
            sys.exit('ERROR: Please input DECmin if you input some RA DEC constraints.')

    if input_list == True:
        
        # Read in lists

        # For each subtraction image, identify candidates
        for sub,sci,neg in zip():
            (region_file, 
            candidate_list_file, 
            candidate_list,
            ds9commands_file) = find_transientCandidates_(sub,sci,neg,savedir=savedir,
                                                        RAmax=RAmax,RAmin=RAmin,DECmax=DECmax,DECmin=DECmin,
                                                        verbose=verbose,debugmode=debugmode,quietmode=quietmode)
            Region_Files.append(region_file)
            Candidate_List_Files.append(candidate_list_file)
            Candidate_Lists.append(candidate_list)
            ds9commands_Files.append(ds9commands_file)

    elif input_list == False:

        # identify canddiates
        (region_file, 
        candidate_list_file, 
        candidate_list,
        ds9commands_list) = find_transientCandidates_(sub,sci,neg,savedir=savedir,
                                                    RAmax=RAmax,RAmin=RAmin,DECmax=DECmax,DECmin=DECmin,
                                                    verbose=verbose,debugmode=debugmode,quietmode=quietmode)    
    # Return appropriate outputs    
    if input_list == True:
        return Region_Files, Candidate_List_Files, Candidate_Lists, ds9commands_Files
    else:
        return region_file, candidate_list_file, candidate_list, ds9commands_list

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
    input_list      = arguments['--list']
    sub_cat         = arguments['<sub_cat>']
    sci_cat         = arguments['<sci_cat>']
    neg_cat         = arguments['<neg_cat>']
    RAmax           = arguments['--RAmax']
    RAmin           = arguments['--RAmin']
    DECmax          = arguments['--DECmax']
    DECmin          = arguments['--DECmin']

    _,_,_,_ = find_transientCandidates(sub_cat, sci_cat, neg_cat, input_list=input_list, savedir=savedir,
                                        RAmax=RAmax,RAmin=RAmin,DECmax=DECmax,DECmin=DECmin,
                                        verbose=verbose,debugmode=debugmode,quietmode=quietmode)
