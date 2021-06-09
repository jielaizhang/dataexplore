#!/usr/bin/python

#!/usr/bin/env python

""" -- 
Usage: XXX [-h] [-v] [--debug] [-q] [-o SAVELOC] [-b STRING] (add | subtract) <blah>

Arguments:
    blah (string)
        Info on blah

Options:
    -h, --help                              Show this screen
    -v, --verbose                           Show extra information [default: False]   
    -q, --quietmode                         Minimize print to screen. This is useful when this function is called in another function. [default: False]  
    --debug                                 Output more for debugging [default: False]
    -o SAVELOC, --out SAVELOC               Saved output as [default: ./profile.png]
    -b STRING, --boundonaxis STRING         xmin,xmax,ymin,ymax e.g. 10,50,0,100; e.g. 10,nan,0,100

Examples:
from datastats.measure_psf import measure_psf
from misc.convert_wavelength_frequency_energy import eVs_to_Hzs
"""
import docopt
import astropy.io.fits as fits
import numpy as np
import sys, os
import ntpath
from pathlib import Path
from datetime import datetime
import string
import random

import copy

__author__      = "Jielai Zhang"
__license__     = "MIT"
__version__     = "1.0.1"
__date__        = "2020-05-07"
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
####################### Main Function ########################
##############################################################

def XXX(verbose=False,debugmode=False,quietmode=False):
    # Determine save directory
    savedir = os.sep.join(saveloc.split(os.sep)[0:-1]) + os.sep

    # Determine saveloc for subtracted image based on input savedir and input file name
    fname       = ntpath.basename(fitspath)
    fname_stub  = Path(fname).stem
    stub        = Path(ntpath.basename(fitspath)).stem
    saveloc     = savedir + os.path.sep + fname_stub + '_sub.fits'

    # Convert to dt
    dt = datetime.strptime(sunset_CT.iso, '%Y-%m-%d %H:%M:%S.%f') #2020-10-23 16:01:22.058 
    
    # Random string
    S=5 # number of characters in random string
    ran = ''.join(random.choices(string.ascii_uppercase + string.digits, k = S))

    # Make temporoary random directory
    temp_dir = './temporary_'+ran
    os.makedirs(temp_dir)

    # Saving fits file
    fits.writeto(output,stellar_mask) 
    printme = f'SAVED  : {output}'
    print(printme)

    # Plot
    fig	= plt.figure()
    ax	= fig.add_subplot(1,1,1)
    # plot something here
    if xyaxislims:
        if ~np.isfinite(xmin):
            xmin = ax.get_xlim()[0]
        if ~np.isfinite(xmax):
            xmax = ax.get_xlim()[1]
        if ~np.isfinite(ymin):
            ymin = ax.get_ylim()[0]
        if ~np.isfinite(ymax):
            ymax = ax.get_ylim()[1]
        plt.ylim([ymin,ymax])
        plt.xlim([xmin,xmax])

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
    saveloc         = arguments['--out']
    input_add       = arguments['add']
    input_subtract  = arguments['subtract']
    blah            = arguments['<blah>']
    xyaxislims      = arguments['--boundonaxis']
    if xyaxislims:
        xmin,xmax,ymin,ymax = np.array(xyaxislims.split(',')).astype(np.float)

    _ = XXX(verbose=verbose,debugmode=debugmode,quietmode=quietmode)
