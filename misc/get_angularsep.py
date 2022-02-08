#!/usr/bin/env python

""" get_angularsep.py -- given 2 RA and DECs, determine and print angular separation. Can input in hms,dms or hms,deg or deg,dms or deg,deg. 
alias getsep='python /Users/jielaizhang/src/dataexplore/misc/get_angularsep.py '

Usage: get_angularsep.py [-h] [-v] [--debug] [-q] <RADEC1> <RADEC2>

Arguments:
    RADEC1 (string)
        E.g. 20:11:01,-3:11:01 or 1.23234,3.21 or 20:11:01,3.21 or 1.23234,-3:11:01
    RADEC2 (string)
        E.g. 20:11:01,-3:11:01 or 1.23234,3.21 or 20:11:01,3.21 or 1.23234,-3:11:01

Options:
    -h, --help                              Show this screen
    -v, --verbose                           Show extra information [default: False]   
    -q, --quietmode                         Minimize print to screen. This is useful when this function is called in another function. [default: False]  
    --debug                                 Output more for debugging [default: False]

Examples: 
Bash:
    python /Users/jielaizhang/src/dataexplore/misc/get_angularsep.py 495.7592,-3.0339
Python:
    from misc.get_angularsep import get_angularsep
"""
import docopt
import astropy.units as u
from astropy.coordinates import SkyCoord

# dataexplore modules
from misc.radec import hms2deg, dms2deg, deg2hms, deg2dms

__author__      = "Jielai Zhang"
__license__     = "MIT"
__version__     = "1.0.1"
__date__        = "2021-01-14"
__maintainer__  = "Jielai Zhang"
__email__       = "zhang.jielai@gmail.com"


##############################################################
####################### Main Functions #######################
##############################################################

def get_angularsep(RA1,DEC1,RA2,DEC2,verbose=False,debugmode=False,quietmode=False):

    # check if input format is ok
    # do this later

    # convert all to decimal 
    if ':' in 'RA1':
        ra1 = hms2deg(RA1)
    else:
        ra1 = float(RA1)
    if ':' in 'RA2':
        ra2 = hms2deg(RA2)
    else:
        ra2 = float(RA2)    
    if ':' in 'DEC1':
        dec1 = dms2deg(DEC1)
    else:
        dec1 = float(DEC1)
    if ':' in 'DEC2':
        dec2 = dms2deg(DEC2)
    else:
        dec2 = float(DEC2)
    
    # Calculate angular separation
    c1 = SkyCoord(ra1*u.deg, dec1*u.deg, frame='icrs')
    c2 = SkyCoord(ra2*u.deg, dec2*u.deg, frame='icrs')
    angular_sep = c1.separation(c2)
    print(f'{angular_sep.deg:.4f} deg |', angular_sep, f'| {angular_sep.arcsec:.3f} arcsec')
    return angular_sep

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
    RADEC1          = arguments['<RADEC1>']
    RA1, DEC1       = RADEC1.split(',')
    RADEC2          = arguments['<RADEC2>']
    RA2, DEC2       = RADEC2.split(',')

    _ = get_angularsep(RA1,DEC1,RA2,DEC2,verbose=verbose,debugmode=debugmode,quietmode=quietmode)
