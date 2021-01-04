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

def write_empty_region_file(region_file):
    # Write region file only with error region. 
    fp = open(region_file, "w")
    fp.write('# Region file format: DS9 version 4.0')
    fp.write('# No candidates detected.')
    fp.write("j2000; circle 326.65958333 -81.0772222 7.5' #color=yellow")
    fp.close()
    return region_file

def write_empty_candidate_list_file(candidate_list_file):
    # Write empty candidate list file
    fp = open(candidate_list_file, "w")
    fp.write('# No candidates detected.\n')
    fp.close()
    return candidate_list_file

def write_empty_ds9commands_file(ds9commands_file,sub_fits,sci_fits,neg_fits):

    # visually inspect sci, sub, neg
    f_ds9 = open(ds9commands_file, "w")
    f_ds9.write("# ====================================================\n")
    f_ds9.write(f'# {sub_cat}')
    command =   f'ds9 -zscale -lock frame wcs {sub_fits} {sci_fits} {neg_fits} &' 
    f_ds9.write(command)
    f_ds9.write('')
    f_ds9.close()

    return ds9commands_file

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
    # If sub catalog is empty, skip
    try:
        dat_pos = ascii.read(sub_cat)
    except:
        write_empty_region_file(region_file)
        write_empty_candidate_list_file(candidate_list_file)
        candidate_list = [] # create empty candidate list for return
        write_empty_ds9commands_file(ds9commands_file,sub_fits,sci_fits,neg_fits)
        print('Empty catalog for subtraction image.')
        return region_file, candidate_list_file, candidate_list, ds9commands_file
    dat_neg = ascii.read(neg_cat) 
    dat_sci = ascii.read(sci_cat) 

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

    no_candidates = False
    if len(df_pos_select) == 0:
        print('No subtraction sources has passed tests up to this point.')
        no_candidates = True
    else:
        # Get sci and sub matched
        radius_threshold = 2.5 * u.arcsec # 10 pixels for KMTNet
        df_sciYes_posYes, df_sciYes_posYes_sub = match_catalogs_df( df_sci_select, df_pos_select,
                                                                    radius_threshold=radius_threshold)
        # Get sci/sub matched with no nearby neg
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
    f_ds9.write("# ====================================================\n")
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

    # Remove temporary region files
    os.remove(reg_pos_temp)
    os.remove(reg_sci_temp)
    os.remove(reg_neg_temp)

    # Add ds9 command to ds9 commands file.
    f_ds9.write('# ---Sub, neg sub, sci, detected sources region file --- \n')
    f_ds9.write('# green = sub\n')
    f_ds9.write('# red = neg sub\n')
    f_ds9.write('# blue 30 deg = sci\n')
    if no_candidates == True:
        f_ds9.write('# No candidates found.\n')
    command =   f'ds9 -zscale -lock frame wcs {sub_fits} {neg_fits} {sci_fits} '\
                f'-regions load all {reg_all_selected_detections} &\n\n' 
    f_ds9.write(command)
    f_ds9.write('')


    # ..............................................................
    # Make Region files of sources matched in Sci and Sub

    if no_candidates == False:
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
        f_ds9.write('# ---Sub, sci, matched sources (sciYes_posYes) region file --- \n')
        f_ds9.write('# green circle = sci/sub matched\n')
        command =   f'ds9 -zscale -lock frame wcs {sub_fits} {sci_fits} '\
                    f'-regions load all {reg_sciYes_subYes} & \n\n'    
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
        f_ds9.write('# ---Sub, sci, neg; sciYes_posYes_negNo region file --- \n')
        f_ds9.write('# yellow panda = sci/sub matched but not matched to neg\n')
        command =   f'ds9 -zscale -lock frame wcs {sub_fits} {sci_fits} {neg_fits} '\
                    f'-regions load all {region_file} &\n\n'    
        f_ds9.write(command)
        f_ds9.write('')

        # Add ds9 command to ds9 commands file- sub sci only
        f_ds9.write('# ---Sub, sci; sciYes_posYes_negNo region file --- \n')
        f_ds9.write('# yellow panda = sci/sub matched but not matched to neg\n')
        command =   f'ds9 -zscale -lock frame wcs {sub_fits} {sci_fits} '\
                    f'-regions load all {region_file} &\n\n'    
        f_ds9.write(command)
        f_ds9.write('')

    elif no_candidates == True:

        # Write region file only with error region. 
        fp = open(region_file, "w")
        fp.write('# Region file format: DS9 version 4.0')
        fp.write('# No candidates detected.')
        fp.write("j2000; circle 326.65958333 -81.0772222 7.5' #color=yellow")
        fp.close()

    # Candidate pandas table
    h               =   'RA_sci DEC_sci RA_sub DEC_sub '\
                        'FLUX_RADIUS_sci FLUX_RADIUS_sub '\
                        'FLUX_APER_sci FLUXERR_APER_sci FLUX_APER_sub FLUXERR_APER_sub '\
                        'MAG_AUTO_sci MAGERR_AUTO_sci MAG_AUTO_sub MAGERR_AUTO_sub '\
                        'MAG_MODEL_sci MAGERR_MODEL_sci MAG_MODEL_sub MAGERR_MODEL_sub '\
                        'BACKGROUND_sci BACKGROUND_sub '\
                        'FWHM_IMAGE_sci FWHM_WORLD_sci FWHM_IMAGE_sub FWHM_WORLD_sub '\
                        'ISOAREAF_IMAGE_sci  ISOAREAF_IMAGE_sub '\
                        'ELLIPTICITY_sci ELLIPTICITY_sub '\
                        'FLAGS_sci FLAGS_sub '\
                        'CLASS_STAR_sci CLASS_STAR_sub '\
                        'SPREAD_MODEL_sci SPREAD_MODEL_sub'

    if no_candidates == False:
        RA_sci              = df_sciYes_posYes_negNo['X_WORLD']
        DEC_sci             = df_sciYes_posYes_negNo['Y_WORLD']
        FLUX_RADIUS_sci     = df_sciYes_posYes_negNo['FLUX_RADIUS']
        FLUX_APER_sci       = df_sciYes_posYes_negNo['FLUX_APER']
        FLUXERR_APER_sci    = df_sciYes_posYes_negNo['FLUXERR_APER']
        MAG_AUTO_sci        = df_sciYes_posYes_negNo['MAG_AUTO']
        MAGERR_AUTO_sci     = df_sciYes_posYes_negNo['MAGERR_AUTO']
        MAG_MODEL_sci       = df_sciYes_posYes_negNo['MAG_MODEL']
        MAGERR_MODEL_sci    = df_sciYes_posYes_negNo['MAGERR_MODEL']
        BACKGROUND_sci      = df_sciYes_posYes_negNo['BACKGROUND']
        FWHM_IMAGE_sci      = df_sciYes_posYes_negNo['FWHM_IMAGE']
        FWHM_WORLD_sci      = df_sciYes_posYes_negNo['FWHM_WORLD']
        ISOAREAF_IMAGE_sci   = df_sciYes_posYes_negNo['ISOAREAF_IMAGE']
        ELLIPTICITY_sci     = df_sciYes_posYes_negNo['ELLIPTICITY']
        FLAGS_sci           = df_sciYes_posYes_negNo['FLAGS']
        CLASS_STAR_sci      = df_sciYes_posYes_negNo['CLASS_STAR']
        SPREAD_MODEL_sci    = df_sciYes_posYes_negNo['SPREAD_MODEL']

        RA_sub              = df_sciYes_posYes_negNo_sub['X_WORLD']
        DEC_sub             = df_sciYes_posYes_negNo_sub['Y_WORLD']
        FLUX_RADIUS_sub     = df_sciYes_posYes_negNo_sub['FLUX_RADIUS']
        FLUX_APER_sub       = df_sciYes_posYes_negNo_sub['FLUX_APER']
        FLUXERR_APER_sub    = df_sciYes_posYes_negNo_sub['FLUXERR_APER']
        MAG_AUTO_sub        = df_sciYes_posYes_negNo_sub['MAG_AUTO']
        MAGERR_AUTO_sub     = df_sciYes_posYes_negNo_sub['MAGERR_AUTO']
        MAG_MODEL_sub       = df_sciYes_posYes_negNo_sub['MAG_MODEL']
        MAGERR_MODEL_sub    = df_sciYes_posYes_negNo_sub['MAGERR_MODEL']
        BACKGROUND_sub      = df_sciYes_posYes_negNo_sub['BACKGROUND']
        FWHM_IMAGE_sub      = df_sciYes_posYes_negNo_sub['FWHM_IMAGE']
        FWHM_WORLD_sub      = df_sciYes_posYes_negNo_sub['FWHM_WORLD']
        ISOAREAF_IMAGE_sub   = df_sciYes_posYes_negNo_sub['ISOAREAF_IMAGE']
        ELLIPTICITY_sub     = df_sciYes_posYes_negNo_sub['ELLIPTICITY']
        FLAGS_sub           = df_sciYes_posYes_negNo_sub['FLAGS']
        CLASS_STAR_sub      = df_sciYes_posYes_negNo_sub['CLASS_STAR']
        SPREAD_MODEL_sub    = df_sciYes_posYes_negNo_sub['SPREAD_MODEL']

        candidate_list  = np.transpose([RA_sci,DEC_sci,RA_sub,DEC_sub,
                                        FLUX_RADIUS_sci, FLUX_RADIUS_sub,
                                        FLUX_APER_sci,FLUXERR_APER_sci, FLUX_APER_sub,FLUXERR_APER_sub,
                                        MAG_AUTO_sci,MAGERR_AUTO_sci,   MAG_AUTO_sub,MAGERR_AUTO_sub,
                                        MAG_MODEL_sci,MAGERR_MODEL_sci, MAG_MODEL_sub,MAGERR_MODEL_sub,
                                        BACKGROUND_sci, BACKGROUND_sub, 
                                        FWHM_IMAGE_sci,FWHM_WORLD_sci,  FWHM_IMAGE_sub,FWHM_WORLD_sub,
                                        ISOAREAF_IMAGE_sci,  ISOAREAF_IMAGE_sub,
                                        ELLIPTICITY_sci, ELLIPTICITY_sub,
                                        FLAGS_sci, FLAGS_sub,
                                        CLASS_STAR_sci,CLASS_STAR_sub,
                                        SPREAD_MODEL_sci, SPREAD_MODEL_sub])
        np.savetxt(candidate_list_file, (candidate_list),fmt='%s',header=h)

    elif no_candidates == True:
        # Write empty candidate list file
        fp = open(candidate_list_file, "w")
        fp.write('# No candidates detected.\n')
        fp.write(h)
        fp.close()

        # Make empty candidate list for output
        candidate_list = []

    if verbose:
        print(f'Saved: {candidate_list_file}') 
    
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
        sub_cat_list = [list(x) for x in ascii.read(sub_cat,format='no_header')]
        sub_cat_list = np.transpose(sub_cat_list)[0]
        sci_cat_list = [list(x) for x in ascii.read(sci_cat,format='no_header')]
        sci_cat_list = np.transpose(sci_cat_list)[0]
        neg_cat_list = [list(x) for x in ascii.read(neg_cat,format='no_header')]
        neg_cat_list = np.transpose(neg_cat_list)[0]
        sub_fits_list = [list(x) for x in ascii.read(sub_fits,format='no_header')]
        sub_fits_list = np.transpose(sub_fits_list)[0]
        sci_fits_list = [list(x) for x in ascii.read(sci_fits,format='no_header')]
        sci_fits_list = np.transpose(sci_fits_list)[0]
        neg_fits_list = [list(x) for x in ascii.read(neg_fits,format='no_header')]
        neg_fits_list = np.transpose(neg_fits_list)[0]
        
        # Count how many
        if verbose:
            print(f'There are {len(sub_cat_list)} files to go through.')

        # For each subtraction image, identify candidates
        upto_number = 0
        for sub,sci,neg,_sub_fits,_sci_fits,_neg_fits in zip(sub_cat_list, sci_cat_list, neg_cat_list,
                                                                                    sub_fits_list,sci_fits_list,neg_fits_list):
            upto_number += 1
            if verbose:
                print('#-------------------------------')
                print(f'...Processing file {upto_number} out of {len(sub_cat_list)}: {sci.split("/")[-1]}.')
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
