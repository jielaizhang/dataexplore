#!/usr/bin/python

#!/usr/bin/env python

""" plot_frameoverlap.py -- Input fits files, and plot their footprint. Save as png.

Usage: plot_frameoverlap.py [-h] [-v] [--debug] [-q] [-p SAVELOC] [--RADECS STRING] <fitsfiles>... 

Arguments:
    fitsfiles (string)        

Options:
    -h, --help                              Show this screen
    -v, --verbose                           Show extra information [default: False]   
    -q, --quietmode                         Minimize print to screen. This is useful when this function is called in another function. [default: False]
    --debug                                 Output more for debugging [default: False]
    -p SAVELOC, --plot SAVELOC              Saved output as [default: ./frameoverlap.png]
    --RADECS STRING                         E.g. hh:mm:ss,+dd:mm:ss;hh:mm:ss,+dd:mm:ss

Examples:
    Bash: python plot_frameoverlap.py f1.fits f2.fits f3.fits
    Python: from datavis.fits.plot_frameoverlap import plot_frameoverlap
    plot_path = plot_frameoverlap( [f1,f2,f3],verbose=False,debugmode=False,quietmode=False)
"""
import docopt
import astropy.io.fits as fits
import numpy as np
from astropy.wcs import wcs
import sys, os
from shapely.geometry.polygon import Polygon
import matplotlib.pyplot as plt

# Jielai modules
from misc.radec import hms2deg, dms2deg

__author__      = "Jielai Zhang"
__license__     = "MIT"
__version__     = "1.0.1"
__date__        = "2020-06-25"
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

#################################
# ======= Other Functions =======
#################################
def get_cornerRADECs(f):
        d,h       = fits.getdata(f,header=True)
        w         = wcs.WCS(h)
        lenx,leny = np.shape(d)
        four_corners_pixels = [[0,0],[0,lenx-1],[leny-1,lenx-1],[leny-1,0]]
        four_corners_RADEC  = w.all_pix2world(four_corners_pixels,0)
        return four_corners_RADEC

def get_imagesize(corner_RADECs,w):
    # Get four corner coordinates in pixel coordinates
    four_corners_pixels = w.all_world2pix(corner_RADECs,0)
    # Get max and min x,y
    max_x               = np.max(four_corners_pixels[:,0])
    min_x               = np.min(four_corners_pixels[:,0])
    XLEN                = max_x-min_x
    max_y               = np.max(four_corners_pixels[:,1])
    min_y               = np.min(four_corners_pixels[:,1])
    YLEN                = max_y-min_y
    return XLEN, YLEN

def plot_polygonOverlap(polys,plotsave,RADECS=None,quietmode=False):

    # Initiate figure
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_title('Image Footprints')

    # Plot polys
    for poly in polys:
        x1,y1 = poly.exterior.xy
        ax.plot(x1, y1, color='blue', alpha=0.7,
                linewidth=0.5, solid_capstyle='round', zorder=2)

    if RADECS:
        RADEC_array = RADECS.split(';')
        for RADEC in RADEC_array:
            RA,DEC=RADEC.split(',')
            RA  = hms2deg([RA])
            DEC = dms2deg([DEC])
            ax.plot(RA,DEC,color='blue',marker='x')
        
    # Save and print if not quietmode 
    plt.savefig(plotsave)
    if not quietmode:
        print(f'SAVED  : {plotsave}')

    return None

##############################################################
####################### Main Function ########################
##############################################################
def plot_frameoverlap( fitsfiles,plotsave='./frameoverlap.png',RADECS=None,
                       verbose=False,debugmode=False,quietmode=False):
    ''' Determine the max number of pixels in AXIS1 and AXIS2 of the input fits images, print to screen unless suppressed. Plot up the footprint if -p option is on. Right now, only takes in 2 fits images (YOLO).'''
    
    polys = []
    for f in fitsfiles:
        d,h = fits.getdata(f,header=True)
        four_corners_RADECs = get_cornerRADECs(f)
        poly =  Polygon(four_corners_RADECs)
        polys.append(poly)
        
    plot_polygonOverlap(polys,plotsave,RADECS=RADECS,quietmode=quietmode)

    return plotsave

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
    fitsfiles       = arguments['<fitsfiles>']
    plotsave        = arguments['--plot']
    RADECS          = arguments['--RADECS']


    _ = plot_frameoverlap(fitsfiles,plotsave=plotsave,RADECS=RADECS,
                          verbose=verbose,debugmode=debugmode,quietmode=quietmode)

