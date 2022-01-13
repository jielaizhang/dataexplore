#!/usr/bin/env python

""" -- 
Usage: radec.py [-h] [-v] [--debug] [-q] (sexa2deci | deci2sexa ) <RADEC>

Arguments:
    RA,DEC (string)
        E.g. 20:11:01,-3:11:01 or 1.23234,3.21
Options:
    -h, --help                              Show this screen
    -v, --verbose                           Show extra information [default: False]   
    -q, --quietmode                         Minimize print to screen. This is useful when this function is called in another function. [default: False]  
    --debug                                 Output more for debugging [default: False]

Examples: 
Bash:
    python /Users/jielaizhang/src/dataexplore/misc/radec.py deci2sexa 495.7592,-3.0339
    python /Users/jielaizhang/src/dataexplore/misc/radec.py sexa2deci 33:3:2.2,-3:2:2.2
Python:
    from dataexplore.misc.radec import sexagesimal2deg
    sexagesimal2deg('20:11:01', '-70:21:12')
    from dataexplore.misc.radec import deg2sexagecimal
    sexagesimal2deg(150.1, 70.5)
    from dataexplore.misc.radec import sexagesimal2degs
    sexagesimal2degs([['20:11:01', '-70:21:12'],
                      ['20:11:01', '-70:21:12']])
    from dataexplore.misc.radec import deg2sexagecimals
    sexagesimal2degs([[150.1, 70.5],
                      [150.1, 70.5]])

    from misc.radec import hms2deg, dms2deg
    hms2deg(['20:10','4:34:23'])
    dms2deg(['-70:10','-4:34:23'])
"""
import docopt
import astropy.io.fits as fits
import numpy as np
import sys, os, ntpath
from pathlib import Path

__author__      = "Jielai Zhang"
__license__     = "MIT"
__version__     = "1.0.1"
__date__        = "2020-08-06"
__maintainer__  = "Jielai Zhang"
__email__       = "zhang.jielai@gmail.com"

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

##############################################################
#################### other RADEC Functions ###################
##############################################################

def minussignornot(x):
    if x >= 0:
        sign = ''
    if x < 0:
        sign = '-'
    return sign

def hms2deg(ras):
    '''hms format: e.g. 12:23:34.32'''
    ras_deg = []
    for ra in ras:
        try:
            rh, rm, rs = [float(r) for r in ra.split(':')]
        except:
            rh, rm     = [float(r) for r in ra.split(':')]
            rs         = 0.0
        ra_deg     = rh*15 + rm/4 + rs/240
        ras_deg.append(ra_deg)
    return ras_deg

def dms2deg(decs):
    '''dms format: e.g. 12:23:34.32 or -12:23:34.32 '''
    decs_deg = []
    for dec in decs:
        try:
            dd, dm, ds = [float(d) for d in dec.split(':')]
        except: 
            dd, dm     = [float(d) for d in dec.split(':')]
            ds         = 0.0
        if dd < 0:
            sign   = -1
        else:
            sign   = 1
        dec_deg    = dd + sign*dm/60 + sign*ds/3600
        decs_deg.append(dec_deg)
    return decs_deg

def deg2dms(decs):
    decs_dms = []
    for dec in decs:
        dec    = float(dec)
        pmsign = minussignornot(dec)
        dec  = abs(dec)
        d  = np.floor(dec)
        m  = np.floor((dec-d)*60.)
        s  = ((dec-d)*60.-m)*60.
        dms = f'{pmsign}{int(d)}:{int(m)}:{s}'
        decs_dms.append(dms)
    return decs_dms

def deg2hms(ras):
    ras_hms = []
    for ra in ras:
        ra = np.float(ra)
        h  = np.floor(ra/15)
        m  = np.floor( (ra/15 - h)*60. )
        s  = ((ra/15 - h)*60. -m)*60. 
        hms = f'{int(h)}:{int(m)}:{s}'
        ras_hms.append(hms)
    return ras_hms

##############################################################
####################### Main Functions #######################
##############################################################

def sexagesimal2deg(RA,DEC,verbose=False,debugmode=False,quietmode=False):
    [ra]  = hms2deg([RA])
    [dec] = dms2deg([DEC])
    print(ra,dec)
    return None

def deg2sexagecimal(RA,DEC,verbose=False,debugmode=False,quietmode=False):
    [ra]  = deg2hms([RA])
    [dec] = deg2dms([DEC])
    print(ra,dec)
    return None

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
    input_deg          = arguments['deci2sexa']
    input_sexagecimal  = arguments['sexa2deci']
    RADEC              = arguments['<RADEC>']
    RA, DEC = RADEC.split(',')

    if input_sexagecimal:
        _ = sexagesimal2deg(RA,DEC,verbose=verbose,debugmode=debugmode,quietmode=quietmode)

    if input_deg:
        _ = deg2sexagecimal(RA,DEC,verbose=verbose,debugmode=debugmode,quietmode=quietmode)
