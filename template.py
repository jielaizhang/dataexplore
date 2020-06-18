#!/usr/bin/python

#!/usr/bin/env python

""" -- 
Usage: XXX [-h] [-v] [--debug] [-o SAVELOC] [-b STRING] (add | subtract) <blah>

Arguments:
    blah (string)
        Info on blah

Options:
    -h, --help                              Show this screen
    -v, --verbose                           Show extra information [default: False]     
    --debug                                 Output more for debugging [default: False]
    -o SAVELOC, --out SAVELOC               Saved output as [default: ./profile.png]
    -b STRING, --boundonaxis STRING         xmin,xmax,ymin,ymax e.g. 10,50,0,100; e.g. 10,nan,0,100

Examples:
"""
import docopt
import astropy.io.fits as fits
import numpy as np
import sys, os
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
        else:
            print(f"VERBOSE: {printme}",file=sys.stdout)

def print_debug_string(printme,debugmode=False,underscores=False):
    if debugmode:
        if underscores:
            print("-" * len(f"DEBUG  : {printme}"),file=sys.stdout)
        else:
            print(f"DEBUG  : {printme}",file=sys.stdout)

'''These functions help organise or delete files'''
def clearit(fname):
    if os.path.isfile(fname):
        os.remove(fname)
    return None

############################################################################
####################### BODY OF PROGRAM STARTS HERE ########################
############################################################################

if __name__ == "__main__":

    # Read in input arguments
    arguments       = docopt.docopt(__doc__)
    verbose         = arguments['--verbose']
    debugmode       = arguments['--debug']
    if debugmode:
        print(arguments)   
    saveloc         = arguments['--out']
    input_add       = arguments['add']
    input_subtract  = arguments['subtract']
    blah            = arguments['<blah>']
    xyaxislims      = arguments['--boundonaxis']
    if xyaxislims:
        xmin,xmax,ymin,ymax = np.array(xyaxislims.split(',')).astype(np.float)

    # Determine save directory
    savedir = os.sep.join(saveloc.split(os.sep)[0:-1]) + os.sep

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
