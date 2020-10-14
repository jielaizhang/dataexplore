#!/usr/bin/env python

# This script is not done, as it doesn't serve my purpose right now. Instead, I will write makeScript_modheadRelativeZPs.py

""" embed_relativezps.py -- Input a reference fits file and a series of other fits files with significant FOV overlap, calculate the relative ZP of the other fits files. 

Usage: embed_relativezps.py [-h] [-v] [--debug] [-q] [-u] [--catdir STRING] <ref> <others>...

Arguments:
    ref (string, fits file)
    others (strings, fits files)

Options:
    -h, --help                              Show this screen
    -v, --verbose                           Show extra information [default: False]   
    -q, --quietmode                         Minimize print to screen. This is useful when this function is called in another function. [default: False]  
    --debug                                 Output more for debugging [default: False]
    -u, --update                            Embed or update into fits header the ref file name and relative zp [default: False]
    --catdir STRING                         Directory with catalogs of the ref and other fits files, and the add on stubb. E.g. 

Examples:
"""
import docopt
import numpy as np
import os
from glob import glob
import pandas as pd
import corner
import matplotlib.pyplot as plt
import pylab
import astropy.io.ascii as ascii
from datetime import datetime
import ntpath

# from Jielai modules
from datavis.ascii.match_catalogs import match_catalogs

__author__      = "Jielai Zhang"
__license__     = "MIT"
__version__     = "1.0.1"
__date__        = "2020-10-13"
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

#########################################
# =========== Other Functions ===========
#########################################
def do_header_update(f,key,value,comment,quietmode=False):
    d,h = fits.getdata(f,header=True)
    if not quietmode:
        if key in h:
            oldval = h[key]
            print(f'{f}: {key}: {oldval}-->{value}')
        else:
            print(f'{f}: new key- {key}: {value}')
    h[key] = (value, comment)
    fits.writeto(f,d,header=h)
    return None

def embed_relativezp(ref_fits_file,other_fits_file,update_header=False,
                    verbose=False,debugmode=False,quietmode=False):
    # See Jupyter Notebook KNe_02_JielaiZhang_RelativePhotometry. 
    # See makeScript_modheadRelativeZPs.py

    return 

##############################################################
####################### Main Function ########################
##############################################################

def embed_relativezps(ref_fits_file,other_fits_files,update_header=False,
                    verbose=False,debugmode=False,quietmode=False):
    keyzp           = 'ZPREL'
    short_ref_name  = ref_fits_file[-60:]
    comment         = f'ZP referenced to: {short_ref_name}'

    # Embed in fits header ZPREL = 0 for ref file
    if update_header:
        do_header_update(ref_fits_file,keyzp,0.0,comment,quietmode=quietmode) 

    # For each input other_fits_files, calculate relative photometry and update header
    for f in other_fits_files:
        

    return None

############################################################################
####################### BODY OF PROGRAM STARTS HERE ########################
############################################################################

if __name__ == "__main__":

    # Read in input arguments
    arguments           = docopt.docopt(__doc__)
    verbose             = arguments['--verbose']
    debugmode           = arguments['--debug']
    quietmode           = arguments['--quietmode']
    if debugmode:
        print(arguments)   
    ref_fits_file       = arguments['<ref>']
    other_fits_files    = arguments['<others>']
    update_header       = arguments['--update']

    _ = embed_relativezps(ref_fits_file,other_fits_files,update_header=update_header,
                         verbose=verbose,debugmode=debugmode,quietmode=quietmode)
