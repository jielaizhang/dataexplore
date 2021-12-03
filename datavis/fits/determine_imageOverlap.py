#!/usr/bin/python

#!/usr/bin/env python

""" determine_imageOverlap.py -- Determine the max number of overlapping pixels in AXIS1 and AXIS2 of the input fits images, print to screen unless suppressed. Plot up the footprint if -p option is on. Right now, only takes in 2 fits images (YOLO).

Known error: if the WCS of the two input images are "flipped" in CDX_X, then the overlap area will not be calculated correctly. 

Usage: determine_imageOverlap.py [-h] [-v] [--debug] [-q] [-p SAVELOC] <fitsfile1> <fitsfile2>

Arguments:
    fitsfile1 (string)
    fitsfile2 (string)
        

Options:
    -h, --help                              Show this screen
    -v, --verbose                           Show extra information [default: False]   
    -q, --quietmode                         Minimize print to screen. This is useful when this function is called in another function. [default: False]
    --debug                                 Output more for debugging [default: False]
    -p SAVELOC, --plot SAVELOC              Saved output as.

Examples:
    Bash: python determine_imageOverlap.py -p ./imageoverlap.png f1.fits f2.fits
    Python: from datavis.fits.determine_imageOverlap import determine_imageOverlap
    XLEN,YLEN = determine_imageOverlap( f1,f2,plotsave=False,
                                        verbose=False,debugmode=False,quietmode=False)
"""
import docopt
import astropy.io.fits as fits
from astropy.wcs import wcs
from astropy.wcs import utils
import numpy as np
import sys, os
from shapely.geometry.polygon import Polygon
import matplotlib.pyplot as plt

import copy

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

def perform_pixelscaleChecks(w1,w2,verbose=False,debugmode=False):
    # Check that the pixel scale in 2 directions on the image are similar, if not exit, as this script isn't coded for it (YOLO)
    pixel_size1 = utils.proj_plane_pixel_scales(w1)
    pixel_size2 = utils.proj_plane_pixel_scales(w2)  
    if (pixel_size1[0]-pixel_size1[1])/pixel_size1[0] > 0.01:
        sys.exit(f'ERROR  : input to determine_imageOverlap, {f1}, has pixels that are not very square, this code is not set up to deal with that YOLO')
    if (pixel_size2[0]-pixel_size2[1])/pixel_size2[0] > 0.01:
        sys.exit(f'ERROR  : input to determine_imageOverlap, {f2}, has pixels that are not very square, this code is not set up to deal with that YOLO')

    # Check the pixel scale of two input images are similar, if not, throw a warning
    if abs(pixel_size1[0]-pixel_size2[0])/pixel_size1[0] < 0.01:
        printme=f'The two image pixel scales are very similar, so this should not be an issue: {pixel_size1},{pixel_size2}'
        print_debug_string(printme,debugmode=debugmode)
    elif abs(pixel_size1[0]-pixel_size2[0])/pixel_size1[0] < 0.1:
        printme=f'The two image pixel scales are within 10%, so this should not be an issue: {pixel_size1},{pixel_size2}'
        print_debug_string(printme,debugmode=debugmode)
    elif abs(pixel_size1[0]-pixel_size2[0])/pixel_size1[0] < 2.0:
        printme=f'The two image pixel scales have a diff > 10% but smaller than a factor of 2, could be an issue: {pixel_size1},{pixel_size2}'
        print_debug_string(printme,debugmode=debugmode) 
    else:
        sys.exit('ERROR  : Input images to determine_imageOverlap differ in pixel scale by more than factor of 2, check if this is ok, if ok, then time to update determine_imageOverlap to allow this.')

    # Print some extra information if in verbose mode
    printme=f'If SWarp is used with PIXELSCALE_TYPE default of MEDIAN, the output image pix scale is the median of pixel scales at centre of input images.'
    print_verbose_string(printme,verbose=verbose)
    return None

def plot_polygonOverlap(poly1,poly2,quietmode=False):

    # Initiate figure
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_title('Image Footprints')

    # Plot poly1 and poly2
    x1,y1 = poly1.exterior.xy
    ax.plot(x1, y1, color='blue', alpha=0.7,
            linewidth=1, solid_capstyle='round', zorder=2)
    x2,y2 = poly2.exterior.xy
    ax.plot(x2, y2, color='green', alpha=0.7,
            linewidth=1, solid_capstyle='round', zorder=2)  

    # Get intersection polygon
    poly_intersect = poly1.intersection(poly2)

    # Get coordinates of corners of the intersection (overlapping) polygon
    x_intersect,y2_intersect = poly_intersect.exterior.xy

    # Plot intersecting polygon
    ax.plot(x_intersect,y2_intersect, color='red', alpha=0.4,
            linewidth=4, solid_capstyle='round', zorder=2) 

    # Save and print if not quietmode 
    plt.savefig(plotsave)
    if not quietmode:
        print(f'SAVED  : {plotsave}')

    return None

##############################################################
####################### Main Function ########################
##############################################################
def determine_imageOverlap( f1,f2,plotsave=False,
                            verbose=False,debugmode=False,quietmode=False):
    ''' Determine the max number of pixels in AXIS1 and AXIS2 of the input fits images, print to screen unless suppressed. Plot up the footprint if -p option is on. Right now, only takes in 2 fits images (YOLO).'''
    
    # Read in f1,f2
    d1,h1 = fits.getdata(f1,header=True)
    d2,h2 = fits.getdata(f2,header=True)
    w1    = wcs.WCS(h1)
    w2    = wcs.WCS(h2)

    # Get the RA and DEC coordinates of the four corners of each image
    four_corners_RADECs1 = get_cornerRADECs(f1)
    four_corners_RADECs2 = get_cornerRADECs(f2)

    # Get polygon of each image corners
    poly1 = Polygon(four_corners_RADECs1)
    poly2 = Polygon(four_corners_RADECs2)

    # Get intersection polygon
    poly_intersect = poly1.intersection(poly2)

    # Get RA DEC of corners of the intersection (overlapping) polygon
    x_intersect,y2_intersect = poly_intersect.exterior.xy

    # Get intersection RA and DEC. 
    # Note that the last coordinate is a repeat of the first, so remove it.
    corner_RADECs_intersect = np.transpose(np.array([x_intersect[0:-1],y2_intersect[0:-1]]))

    # Do pixel scale checks: are pixels square, are pixel scales of two images similar enough
    # Exits if input images found to be outside what this is coded for; Prints warnings if in debugmode but not exiting.
    perform_pixelscaleChecks(w1,w2,verbose=verbose,debugmode=debugmode)

    # Get XLEN and YLEN- Diff pix scales give diff XLEN, YLEN for same RA,DEC range
    # Get the average XLEN and YLEN for the two input image pixscales
    XLEN1,YLEN1 = get_imagesize(corner_RADECs_intersect,w1)
    XLEN2,YLEN2 = get_imagesize(corner_RADECs_intersect,w2)
    if verbose:
        print(f'XLEN1,YLEN1: {XLEN1},{YLEN1}')
        print(f'XLEN2,YLEN2: {XLEN2},{YLEN2}')
    #XLEN        = int(np.average([XLEN1,XLEN2])) # auto floor
    #YLEN        = int(np.average([YLEN1,YLEN2])) # auto floor
    XLEN        = min([XLEN1,XLEN2]) # Take the smallest
    YLEN        = min([YLEN1,YLEN2]) # Take the smallest

    #  Print information if not quiet mode, or if verbose
    if not quietmode and not verbose:
        print(f'-IMAGE_SIZE option in SWarp can be set to: {XLEN},{YLEN}')
    printme=f'SWarp Command to align two images: swarp {f1} {f2}  -SUBTRACT_BACK N -RESAMEPLE Y -COMBINE N -CENTER_TYPE MOST -IMAGE_SIZE {XLEN},{YLEN} -RESAMPLE_DIR ./'
    print_verbose_string(printme,verbose=verbose)

    # If plotsave, plot and save.
    if plotsave:
        plot_polygonOverlap(poly1,poly2,quietmode=quietmode)

    return XLEN,YLEN

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
    fitsfile1       = arguments['<fitsfile1>']
    fitsfile2       = arguments['<fitsfile2>']
    plotsave        = arguments['--plot']


    _ = determine_imageOverlap( fitsfile1,fitsfile2,plotsave=plotsave,
                                verbose=verbose,debugmode=debugmode,quietmode=quietmode)
