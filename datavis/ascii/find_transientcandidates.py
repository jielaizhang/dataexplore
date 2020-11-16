#!/usr/bin/python

#!/usr/bin/env python

""" 
find_transientCandidates.py -- Find transient candidates in subtraction image catalogs. Transient candidates turn up in Science AND Subtraction Image, but not in Neg of Subtraction image. 

Inputs 
    (1) subtraction image source extractor ASCII catalogs (with SPREAD_MODEL)
    (2) corresponding science source extractor ASCII catalogs (with SPREAD_MODEL) 
    (3) corresponding -1*subtraction source extractor ASCII catalogs (without SPREAD_MODEL)
    (4) Subtraction Image path
    (5) Science Image path
    (6) Negative Image path
    Note: fits files are input here, but they are not read in at all. The path is only used to write ds9 commands.
    Note: Each input is a path to the file if --list is False
    Note: Each input is a path to a text file list if --list is True
    !!! Note: If --list, all lists needs to be MATCHED. I.e. e.g. 5th in each list corresponds to each other.    

Outputs for each input trio
    (1) Path to a text file of list of candidate transients.  
    (2) Path to a region file for detected transient candidates. 
    (3) Path to a text file containing ds9 commands. Commands include 
        (i) seeing all detected and selected sources in 3 input catalogs
        (ii) seeing sources in sci AND sub
        (iii) seeing sources in sci AND sub AND not neg
    (4) Array of candidate lists

 

Caveats: 
!!! Currently source ID parameters are optimised for KMTNet data for DWF 
!!! The separation distance between sci/sub and sci/neg is set optimally for KTMNet right now. It is HARD CODED in! 
!!! Error Region of FRB200914 is hard coded in right now. 


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
    --list                                  If specified, the required inputs should be lists [default: False]

Examples:
python ~/src/dataexplore/datavis/ascii/find_transientcandidates.py -v -o ~/Desktop/ --RAmin 325.72868375 --RAmax 327.6204304166667 --DECmax -80.9290338888889 --DECmin -81.22332805555556 sub.cat sci.cat neg.cat sub.fits sci.fits neg.fits
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

def find_transientCandidates_(sub_cat, sci_cat, neg_cat, 
                              sub_fits,sci_fits,neg_fits,
                              savedir='./',
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
        print('#-------------------------------')
        print(f'For {sub_cat}:')
        print(f'Sources in SCI: {len(df_sci)}')
        print(f'Sources in POS: {len(df_pos)}')
        print(f'Sources in NEG: {len(df_neg)}')

    # ..............................................................
    # Select sources that's within RA DEC Range
    if RAmax:
        RAmax = float(RAmax)
        RAmin = float(RAmin)
        DECmax = float(DECmax)
        DECmin = float(DECmin)

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
    # XXX Hard coded section !!! (FWHM, ELLIP, SPREAD_MODEL etc cuts)
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
    # XXX Hard coded section !!! (radius_threshold)

    radius_threshold = 2.5 * u.arcsec # 10 pixels for KMTNet
    df_sciYes_posYes, df_sciYes_posYes_sub = match_catalogs_df( df_sci_select, df_pos_select,
                                                                radius_threshold=radius_threshold)
    radius_threshold = 4 * u.arcsec
    df_sciYes_posYes_negNo, df_sciYes_posYes_negNo_neg = match_catalogs_df( df_sciYes_posYes, df_neg_select,
                                                                            radius_threshold=radius_threshold,find_close=False)
    # Grab the corresponding sub too (this line doesn't change df_sciYes_posYes_negNo)
    df_sciYes_posYes_negNo, df_sciYes_posYes_negNo_sub = match_catalogs_df( df_sciYes_posYes_negNo, df_pos_select,
                                                                            radius_threshold=radius_threshold,find_close=True)

    if verbose:
        print(f'Sources in SCI matched in POS SUB: {len(df_sciYes_posYes)}')
        print(f'Sources in SCI matched in POS SUB, then not matched in NEG: {len(df_sciYes_posYes_negNo)}')


    # ..............................................................
    # Open ds9commands file for overwriting
    f_ds9 = open(ds9commands_file, "w")
    f_ds9.write("# ====================================================")
    f_ds9.write(f'# {sub_cat}')

    # ..............................................................
    # Make Region files of SELECTED sources in Sci, Sub, Neg

    # Sub Detections Selected (Green squares)
    Rs = df_pos_select['X_WORLD']
    Ds = df_pos_select['Y_WORLD']
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
    Rs = df_sci_select['X_WORLD']
    Ds = df_sci_select['Y_WORLD']
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
    if verbose:
        print(f'Saved: {reg_all_selected_detections}')

    # Add ds9 command to ds9 commands file.
    f_ds9.write('# ---Sub, neg sub, sci, detected sources region file --- ')
    f_ds9.write('# green = sub')
    f_ds9.write('# red = neg sub')
    f_ds9.write('# blue 30 deg = sci')
    command =   f'ds9 -zscale -lock frame wcs {sub_fits} {neg_fits} {sci_fits} '\
                f'-regions load all {reg_all_selected_detections} &' 
    f_ds9.write(command)
    f_ds9.write('')


    # ..............................................................
    # Make Region files of sources matched in Sci and Sub

    # Save region file.
    Rs_sciYes_posYes = df_sciYes_posYes['X_WORLD']
    Ds_sciYes_posYes = df_sciYes_posYes['Y_WORLD']
    Sizes   = ["0.1'"]*len(Rs_sciYes_posYes)
    Shapes = ['circle']*len(Rs_sciYes_posYes)
    j2000s  = ['j2000;']*len(Rs_sciYes_posYes)
    Colors = ['#color=green']*len(Rs_sciYes_posYes)
    h       = 'Region file format: DS9 version 4.0'
    savetext = np.transpose([j2000s,Shapes,Rs_sciYes_posYes,Ds_sciYes_posYes,Sizes])
    np.savetxt(reg_sciYes_subYes, (savetext),fmt='%s',header=h)
    if verbose:
        print(f'Saved: {reg_sciYes_subYes}')

    # Add ds9 command to ds9 commands file.
    f_ds9.write('# ---Sub, sci, matched sources (sciYes_posYes) region file --- ')
    f_ds9.write('# green circle = sci/sub matched')
    command =   f'ds9 -zscale -lock frame wcs {sub_fits} {sci_fits} '\
                f'-regions load all {reg_sciYes_subYes} & \n'    
    f_ds9.write(command)
    f_ds9.write('')


    # ..............................................................
    # Make Region files of sources matched in Sci and Sub and not in Neg
    # XXX Hard coded section !!! (FRB200914 error region)

    # Save region file, include FRB200914 error region 
    Rs_sciYes_posYes_negNo = np.array(df_sciYes_posYes_negNo['X_WORLD'])
    Ds_sciYes_posYes_negNo = np.array(df_sciYes_posYes_negNo['Y_WORLD'])
    Shapes  = ['panda']*len(Rs_sciYes_posYes_negNo)
    j2000s  = ['j2000;']*len(Rs_sciYes_posYes_negNo)
    Colors  = ['#color=yellow']*len(Rs_sciYes_posYes_negNo)
    startangles=['0']*len(Rs_sciYes_posYes_negNo)
    stopangles =['360']*len(Rs_sciYes_posYes_negNo)
    nangles    =['4']*len(Rs_sciYes_posYes_negNo)
    inners     =['0.003']*len(Rs_sciYes_posYes_negNo)
    outers     =['0.0035']*len(Rs_sciYes_posYes_negNo)
    nradiuss   =['1']*len(Rs_sciYes_posYes_negNo)
    h       = 'Region file format: DS9 version 4.0'
    savetext = np.transpose([j2000s,Shapes,Rs_sciYes_posYes_negNo,Ds_sciYes_posYes_negNo,
                             startangles,stopangles,nangles,inners,outers,nradiuss,Colors])
    np.savetxt(region_file, (savetext),fmt='%s',header=h)
    fp = open(region_file, "a")
    fp.write("j2000; circle 326.65958333 -81.0772222 7.5' #color=yellow")
    fp.close()
    if verbose:
        print(f'Saved: {region_file}')

    # Add ds9 command to ds9 commands file- sub sci and neg
    f_ds9.write('# ---Sub, sci, neg; sciYes_posYes_negNo region file --- ')
    f_ds9.write('# yellow panda = sci/sub matched but not matched to neg')
    command =   f'ds9 -zscale -lock frame wcs {sub_fits} {sci_fits} {neg_fits} '\
                f'-regions load all {region_file} &'    
    f_ds9.write(command)
    f_ds9.write('')

    # Add ds9 command to ds9 commands file- sub sci only
    f_ds9.write('# ---Sub, sci; sciYes_posYes_negNo region file --- ')
    f_ds9.write('# yellow panda = sci/sub matched but not matched to neg')
    command =   f'ds9 -zscale -lock frame wcs {sub_fits} {sci_fits} '\
                f'-regions load all {region_file} &'    
    f_ds9.write(command)
    f_ds9.write('')


    # Candidate pandas table
    
    candidate_list = df_sciYes_posYes_negNo


#df_sciYes_posYes_negNo_sub[['X_WORLD','Y_WORLD','FLUX_RADIUS','FLUX_APER','ISOAREA_IMAGE','FLAGS','CLASS_STAR', 
#                            'MAG_AUTO','BACKGROUND', 'FWHM_IMAGE','ELLIPTICITY','SPREAD_MODEL']]
#df_sciYes_posYes_negNo[['X_WORLD','Y_WORLD','FLUX_RADIUS','FLUX_APER','ISOAREA_IMAGE','FLAGS','CLASS_STAR', 
#                            'MAG_AUTO','BACKGROUND', 'FWHM_IMAGE','ELLIPTICITY','SPREAD_MODEL']]
    
    # ..............................................................
    # Close ds9commands file 
    f_ds9.close()  
    if verbose:
        print(f'Saved: {ds9commands_file}')


    # ..............................................................
    return region_file, candidate_list_file, candidate_list, ds9commands_file

##############################################################
####################### Main Function ########################
##############################################################

def find_transientCandidates(sub_cat, sci_cat, neg_cat, sub_fits, sci_fits, neg_fits,
                            input_list=False, savedir='./',
                            RAmax=None,RAmin=None,DECmax=None,DECmin=None,
                            verbose=False,debugmode=False,quietmode=False):

    # If savedir doesn't exist, make it
    if not os.path.isdir(savedir):
        os.makedirs(savedir)

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

        # Initiate output arrays
        Region_Files         = []
        Candidate_List_Files = []
        Candidate_Lists      = []
        ds9commands_Files    = []
        
        # Read in lists
        sub_cat_list = [list(x) for x in read(sub_cat,format='no_header')]
        sub_cat_list = np.transpose(sub_cat_list)[0]
        sci_cat_list = [list(x) for x in read(sci_cat,format='no_header')]
        sci_cat_list = np.transpose(sci_cat_list)[0]
        neg_cat_list = [list(x) for x in read(neg_cat,format='no_header')]
        neg_cat_list = np.transpose(neg_cat_list)[0]
        sub_fits_list = [list(x) for x in read(sub_fits,format='no_header')]
        sub_fits_list = np.transpose(sub_fits_list)[0]
        sci_fits_list = [list(x) for x in read(sci_fits,format='no_header')]
        sci_fits_list = np.transpose(sci_fits_list)[0]
        neg_fits_list = [list(x) for x in read(neg_fits,format='no_header')]
        neg_fits_list = np.transpose(neg_fits_list)[0]
        

        # For each subtraction image, identify candidates
        for sub,sci,neg,_sub_fits,_sci_fits,_neg_fits in zip(sub_cat_list, sci_cat_list, neg_cat_list,
                                                             sub_fits_list,sci_fits_list,neg_fits_list):
            (region_file, 
            candidate_list_file, 
            candidate_list,
            ds9commands_file) = find_transientCandidates_(sub,sci,neg,_sub_fits,_sci_fits,_neg_fits,
                                                        savedir=savedir,
                                                        RAmax=RAmax,RAmin=RAmin,DECmax=DECmax,DECmin=DECmin,
                                                        verbose=verbose,debugmode=debugmode,quietmode=quietmode)
            Region_Files.append(region_file)
            Candidate_List_Files.append(candidate_list_file)
            Candidate_Lists.append(candidate_list)
            ds9commands_Files.append(ds9commands_file)

    elif input_list == False:

        # identify candiates
        (region_file, 
        candidate_list_file, 
        candidate_list,
        ds9commands_list) = find_transientCandidates_(sub_cat,sci_cat,neg_cat,sub_fits,sci_fits,neg_fits,
                                                    savedir=savedir,
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
    sub_fits        = arguments['<sub_fits>']
    sci_fits        = arguments['<sci_fits>']
    neg_fits        = arguments['<neg_fits>']    
    RAmax           = arguments['--RAmax']
    RAmin           = arguments['--RAmin']
    DECmax          = arguments['--DECmax']
    DECmin          = arguments['--DECmin']

    _,_,_,_ = find_transientCandidates(sub_cat, sci_cat, neg_cat, 
                                        sub_fits,sci_fits,neg_fits,
                                        input_list=input_list, savedir=savedir,
                                        RAmax=RAmax,RAmin=RAmin,DECmax=DECmax,DECmin=DECmin,
                                        verbose=verbose,debugmode=debugmode,quietmode=quietmode)
