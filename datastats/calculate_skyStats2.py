#!/usr/bin/python

"""calculate_skyStats.py -- Mask stars aggressively, determine sky background stats in remaining pixels. Stats are
    Average
    Median
    Standard deviation of multiple sampling of the sky value  

    Three options for sampling the sky are:
    (1) input --box s -n n: calculate sky in n boxes of size sxs randomly placed around the image
    (2) input --annulus 'xc1,yc1,a1,b1,ang1,xc2,yc2,a2,b2,ang2' -n n: calculate the sky in n elliptical annuli placed randomly. xc1,yc1,a1,b1,ang1 specifies the inner ellipse, and xc2,yc2,a2,b2,ang2 specifies the outer ellipse of the annuli. 
    (3) input --annulusallover 'xc1,yc1,a1,b1,ang1,xc2,yc2,a2,b2,ang2' -n n: calculate the sky in n elliptical annuli placed randomly, but letting it move around a little more than --annulus option
    If no option is set, then the average, median and standard deviation of the unmasked pixels in the image is calculated. 

    There is the option to input a mask to be combined with the mask produced by this script to mask stars.

    Output are XXX

Required input 
    fitsimage - image for which background stats are desired.

Usage:
    calculate_skystats.py [-h] [-v] [-b STRING] [-a STRING] [--annulusallover STRING] [-n INT] [-m FILE] [-s SEXLOC] <fitsimage>

Options:
    -h, --help                          Print this screen.
    -v, --verbose                       Print extra information [default: False]
    --debug                             Print extra extra information and save extra files [default: False]
    -b STRING, --box STRING             Input box parameters: size of box (pixels), number of random boxes to be placed. E.g. 10,100.
    -a STRING, --annulus STRING         Select annulus with params 'xc1,yc1,a1,b1,ang1,xc2,yc2,a2,b2,ang2', angle in degrees, counter clockwise rotation; place random annuli around galaxy.
    --annulusallover STRING             Select annulus with params 'xc1,yc1,a1,b1,ang1,xc2,yc2,a2,b2,ang2', angle in degrees, counter clockwise rotation; place annulus around galaxy, but let it move around a little more than above option.
    -n INT, --niterations INT           Input number of random annuli to be placed. (not used for boxes). [default: 100]
    -m FILE, --mask FILE                Required. Input mask to be combined with grown sextractor mask. Important to input a mask that masks out the galaxy and its extended low surface brightness features!

Example:
    python calculate_skyStats.py -v -b 10 -n 100 fitsimage.fits
    python calculate_skyStats.py -v -b 10 -n 100 -m extra_mask.fits fitsimage.fits
    python calculate_skyStats.py -v --annulus 500,500,150,150,500,500,200,200 -n 100 fitsimage.fits
    python calculate_skyStats.py -v --annulusallover 500,500,150,150,500,500,200,200 -n 100 fitsimage.fits
"""

import docopt
import numpy as np
import astropy.io.fits as fits
from scipy import ndimage
import subprocess
import os, sys, copy

import matplotlib.pyplot as plt
import aplpy
from astropy.visualization import ZScaleInterval

# modules by Jielai Zhang
from datavis.fits import create_stellarMask


###########################################
# ======= Sky Calculation Functions =======
###########################################

def calculate_sky_box(fitsimage,boxsize_pix,n_iterations,extra_mask_file = extra_mask_file,verbose=verbose)
#calculate_sky_box(fitsimage,image,total_mask,boxsize_pix,nboxes):
    '''Place n_iterations number of boxsize_pix sized boxes randomly in image with total_mask, 
    calculate average, median in each box and standard deviation of averages. 
    Also calculate the average, median and std of all unmasked pixels.
    Output 
    sky         = average of median of boxes 
    sky_error   = std of median of boxes
    sky_pixel_std = std of all unmasked pixels.
    
    if verbose, output
    figure to show where boxes were placed'''
    
    sky_counts  = []
    pix_counts  = []
    n_counter   = 0
    n_notfinite = 0

    # Read in image and header
    f,h = fits.getdata(fitsimage, header=True)

    # Set boxes to be placed not too near edge
    xmin = 1.5*boxsize_pix
    ymin = 1.5*boxsize_pix
    xmax = float(h['NAXIS1'])-1.5*boxsize_pix
    ymax = float(h['NAXIS2'])-1.5*boxsize_pix

    # Start figure to plot up box locations
    fig = plt.figure(figsize=(48, 36))
    f1 = aplpy.FITSFigure(fitsimage,figure=fig)
    f1.ticks.hide()
    f1.tick_labels.hide_x()
    f1.tick_labels.hide_y()
    f1.axis_labels.hide()
    interval = ZScaleInterval()
    vmin,vmax = interval.get_limits(image)
    f1.show_grayscale(invert=True, stretch='linear', vmin=vmin, vmax=vmax)    

    xlen = int(h['NAXIS2'])
    ylen = int(h['NAXIS1'])
    xtomesh = np.arange(0, ylen, 1)
    ytomesh = np.arange(0, xlen, 1)
    X, Y    = np.meshgrid(xtomesh, ytomesh)
    
    while n_counter <= nboxes:
        
        # Choose a random spot
        row = np.random.randint(low=ymin,high=ymax)
        col = np.random.randint(low=xmin,high=xmax)

        # Make a box
        image_box = image[row-int(boxsize_pix/2):row+int(boxsize_pix/2)+1,col-int(boxsize_pix/2):col+int(boxsize_pix/2)+1]
        mask_box  = total_mask[row-int(boxsize_pix/2):row+int(boxsize_pix/2)+1,col-int(boxsize_pix/2):col+int(boxsize_pix/2)+1]

        # Plot up location of box for display using show_contour
        display_mask = np.zeros((xlen,ylen))
        display_mask[row-int(boxsize_pix/2):row+int(boxsize_pix/2)+1,col-int(boxsize_pix/2):col+int(boxsize_pix/2)+1] = 1.0
        CS = plt.contour(X, Y, display_mask,linewidths=1.0,alpha=0.1,colors='red')

        # Measure average counts in this masked box
        counts = np.ma.mean(np.ma.masked_array(image_box,mask=mask_box))
        # Measure number of pixels not masked in this masked box
        no_pixels_notmasked = np.sum(mask_box)

        # Add average to sky_counts if finite
        # Also increment box count 
        # Else increment n_notfinite
        if np.isfinite(counts):
            sky_counts.append(counts)
            pix_counts.append(no_pixels_notmasked)
            n_counter += 1
        else:
            n_notfinite += 1

    # Save figure to of annuli locations
    outname = './skyregionlocs.png'
    f1.save(outname)
    print(' ')
    print('***OUTPUT: Box location plot saved here: ',outname)

        printme=f'Number of attempts where average sky count in box was not finite: {n_notfinite}'
        print_verbose_string(printme,verbose=verbose)

    
    return sky, sky_error, sky_pixel_std


###############################
# ======= Main function =======
###############################

def calculate_skyStats.py(  fitsimage, 
                            place_boxes          = place_boxes, 
                            place_annuli         = place_annuli, 
                            place_annuli_allover = place_annuli_allover,
                            n_iterations         = n_iterations,
                            maskfile             = maskfile,
                            sextractorloc        = sextractorloc,
                            verbose = verbose, debugmode = debugmode):
    '''calculate_skyStats.py -- Mask stars aggressively, determine sky background stats in remaining pixels. Stats are
    Average
    Median
    Standard deviation of multiple sampling of the sky value  

    Three options for sampling the sky are:
    (1) input --box s -n n: calculate sky in n boxes of size sxs randomly placed around the image
    (2) input --annulus 'xc1,yc1,a1,b1,ang1,xc2,yc2,a2,b2,ang2' -n n: calculate the sky in n elliptical annuli placed randomly. xc1,yc1,a1,b1,ang1 specifies the inner ellipse, and xc2,yc2,a2,b2,ang2 specifies the outer ellipse of the annuli. 
    (3) input --annulusallover 'xc1,yc1,a1,b1,ang1,xc2,yc2,a2,b2,ang2' -n n: calculate the sky in n elliptical annuli placed randomly, but letting it move around a little more than --annulus option
    If no option is set, then the average, median and standard deviation of the unmasked pixels in the image is calculated. 

    There is the option to input a mask to be combined with the mask produced by this script to mask stars.

    Output are XXX
    '''

    # Check which option(s) were set. If more than 1, ask user to re-set, if zero, set place_no_shapes to True.
    place_no_shapes = False
    options_true = [place_boxes, place_annuli, place_annuli_allover]
    if sum(options_true) > 1:
        error_message = 'ERROR: Please only select one of --boxes, --annulus, and --annulusallover'
        sys.exit(error_message)
    if sum(options_true) == 0:
        place_no_shapes = True
        printme='No option for placing shapes for calculating sky was set, so the avg, median, std of all unmasked pixels in the image will be calculated.'
        print_verbose_string(printme,verbose=verbose)
        
    # ===== --box =====
    # Place nboxes random boxsize_pix pixel sized boxes on total_masked image
    # Calculate average in each box and
    # Standard deviation of these averages
    if place_boxes:
        boxsize_pix     = int(place_boxes)
        sky, sky_error, sky_pixel_variance = calculate_sky_box(fitsimage,boxsize_pix,n_iterations,extra_mask_file = extra_mask_file,verbose=verbose)

    # ===== --annulus =====

    # ===== --annulusallover =====

    # ===== no option set, just calculate stats for all unmasked pixels =====


    return None



###########################################################
###################### Start of main ######################
###########################################################
if __name__=='__main__':

    # Import arguments
    arguments = docopt.docopt(__doc__)
    
    fitsimage           = arguments['<fitsimage>']
    verbose             = arguments['--verbose']
    debugmode           = arguments['--debug']
    place_boxes         = arguments['--box']
    place_annuli        = arguments['--annulus']    
    place_annuli_allover = arguments['--annulusallover']
    n_iterations        = int(arguments['--niterations'])
    extra_mask_file     = arguments['--mask']
    sextractorloc       = arguments['--sextractor']

    if debugmode:
        print(arguments)

    calculate_skyStats.py(  fitsimage, 
                            place_boxes          = place_boxes, 
                            place_annuli         = place_annuli, 
                            place_annuli_allover = place_annuli_allover,
                            n_iterations         = n_iterations,
                            maskfile             = maskfile,
                            sextractorloc        = sextractorloc,
                            verbose = verbose, debugmode = debugmode)

