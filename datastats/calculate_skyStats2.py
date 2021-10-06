#!/usr/bin/python

"""calculate_skyStats.py -- Mask stars aggressively, determine sky background stats in remaining pixels. Stats are
    sky:            Average of many sky box/annuli medians (which are each a sample of the sky value) 
    sky_error:      Standard deviation of many sky box/annuli medians (which are each a sample of the sky value) 
    sky_pixel_std:  Standard deviation of all non-masked pixels  

    Three options for sampling the sky are:
    (1) input --box s -n n: calculate sky in n boxes of size sxs randomly placed around the image
    (2) input --annulus 'xc1,yc1,a1,b1,ang1,xc2,yc2,a2,b2,ang2' -n n: calculate the sky in n elliptical annuli placed randomly. xc1,yc1,a1,b1,ang1 specifies the inner ellipse, and xc2,yc2,a2,b2,ang2 specifies the outer ellipse of the annuli. 
    (3) input --annulusallover 'xc1,yc1,a1,b1,ang1,xc2,yc2,a2,b2,ang2' -n n: as above, but letting it move around a little more than --annulus option
    If no option is set, then the average, median and standard deviation of the unmasked pixels in the image is calculated. 

    Option to input additional mask (mask should be zeros and ones; zero pixels are not masked, one pixels are masked). Important to input a mask that masks out the galaxy and its extended low surface brightness features!

    Option to output checkims: originalfilepath_skymask_contour.png and originalfilepath_skyregionlocs.png (overwrites)

Required input 
    fitsimage - image for which background stats are desired.

Usage:
    calculate_skystats.py [-h] [-q] [-v] [--debug] [--sextractorloc LOC] [-b STRING] [-a STRING] [--annulusallover STRING] [-n INT] [-m FILE] [--checkims] <fitsimage>

Options:
    -h, --help                          Print this screen.
    -q, --quietmode                     Do not print calculated sky values [default: False]
    -v, --verbose                       Print extra information [default: False]
    --debug                             Print extra extra information and save extra files [default: False]
    --sextractorloc LOC                 Source-extractor path [default: /opt/local/bin/source-extractor]
    -b STRING, --box STRING             Input box parameters: size of box (pixels). E.g. 10.
    -a STRING, --annulus STRING         Select annulus with params 'xc1,yc1,a1,b1,ang1,xc2,yc2,a2,b2,ang2', angle in degrees, counter clockwise rotation; place random annuli around inital specs.
    --annulusallover STRING             As above, but when placing annulus around galaxy, let it move around a little more than above option.
    -n INT, --niterations INT           Input number of random annuli to be placed. (also used for boxes). [default: 100]
    -m FILE, --mask FILE                Input mask to be combined with program calculated source mask. 
    -c, --checkims                      Output two check images for masking and sky regions used. [default: False]

Example:
    Bash:
    python calculate_skyStats.py -v -b 10 -n 100 fitsimage.fits
    python calculate_skyStats.py -v -b 10 -n 100 -m extra_mask.fits fitsimage.fits
    python calculate_skyStats.py -v --annulus 500,500,150,150,500,500,200,200 -n 100 fitsimage.fits
    python calculate_skyStats.py -v --annulusallover 500,500,150,150,500,500,200,200 -n 100 fitsimage.fits
    Python:
    from datastats.calculate_skyStats2 import calculate_skyStats
    sky, sky_error, sky_pixel_std = calculate_skyStats( fitsimage, 
                                                        place_boxes          = False|boxsize, 
                                                        place_annuli         = False|'xc1,yc1,a1,b1,ang1,xc2,yc2,a2,b2,ang2', 
                                                        place_annuli_allover = False|'xc1,yc1,a1,b1,ang1,xc2,yc2,a2,b2,ang2',
                                                        n_iterations         = 100 (or some other number),
                                                        input_mask_file      = False|'./mask.fits',
                                                        verbose = False|True, debugmode = False|True)
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
from datavis.fits.create_sourceMask import create_sourceMask

#########################################
# ======= House Keeping Functions =======
#########################################
'''These functions standardise verbose, debug printing'''
def print_verbose_string(printme,verbose=False,underscores=False):
    if verbose:
        if underscores:
            print("-" * len(f"VERBOSE: {printme}"),file=sys.stdout)
            print(f"VERBOSE: {printme}",file=sys.stdout)
        else:
            print(f"VERBOSE: {printme}",file=sys.stdout)

def print_debug_string(printme,debugmode=False,underscores=False):
    if debugmode:
        if underscores:
            print("-" * len(f"DEBUG  : {printme}"),file=sys.stdout)
            print(f"DEBUG  : {printme}",file=sys.stdout)
        else:
            print(f"DEBUG  : {printme}",file=sys.stdout)

'''These functions help organise or delete files'''
def clearit(fname):
    if os.path.isfile(fname):
        os.remove(fname)
    return None

##########################################
# ======= Annuli Related Functions =======
##########################################

def read_annulusparams(annulusparams):
    '''Read out annulus parameters of form xc1,yc1,a1,b1,ang1,xc2,yc2,a2,b2,ang2'''
    params = annulusparams.split(',')
    xc1,yc1,a1,b1,ang1,xc2,yc2,a2,b2,ang2 = params
    return float(xc1),float(yc1),float(a1),float(b1),float(ang1),float(xc2),float(yc2),float(a2),float(b2),float(ang2)

def make_annulus_mask(xlen,ylen,xc1,yc1,a1,b1,ang1,xc2,yc2,a2,b2,ang2):
    '''Read in annulus parameters and create grabber of annulus (1 inside and 0 outside)'''
    ang1_rad = (ang1/360.)*2*np.pi
    ang2_rad = (ang2/360.)*2*np.pi

    # Ellipse 1
    mask1 = np.zeros((xlen,ylen))
    xv,yv = np.meshgrid(np.linspace(0,xlen-1,xlen),np.linspace(0,ylen-1,ylen))
    
    A   = ( (xv-xc1)*np.cos(ang1_rad) + (yv-yc1)*np.sin(ang1_rad) )**2 / a1**2
    B   = ( (xv-xc1)*np.sin(ang1_rad) - (yv-yc1)*np.cos(ang1_rad) )**2 / b1**2
    xi,yi    = np.where( A+B < 1.0 )
    mask1[xi,yi] = 1

    # Ellipse 2
    mask2 = np.zeros((xlen,ylen))
    
    A   = ( (xv-xc2)*np.cos(ang2_rad) + (yv-yc2)*np.sin(ang2_rad) )**2 / a2**2
    B   = ( (xv-xc2)*np.sin(ang2_rad) - (yv-yc2)*np.cos(ang2_rad) )**2 / b2**2
    xi,yi    = np.where( A+B < 1.0 )
    mask2[xi,yi] = 1

    # Combine Ellipse 1 and 2 --> annulus
    mask3   = np.ones((xlen,ylen)).astype(int)
    tmp     = mask1+mask2
    xi,yi   = np.where(tmp == 1.0)
    mask3[xi,yi]    = 0

    return mask3.astype(bool)

##################################################################################
# ======= Other Functions (used for both box/annulus sky sampling options) =======
##################################################################################

def combine_masks(source_mask,input_mask_file,verbose=False):

    # If extra mask file, combine masks
    if input_mask_file:
        input_mask                  = fits.getdata(input_mask_file)
        total_mask                  = source_mask + input_mask
        total_mask                  = total_mask[total_mask!=0] == 1.
        total_mask                  = total_mask == 1 # Change to boolean, True where 1
        printme = 'Input mask is combined with program created source mask.'
        print_verbose_string(printme,verbose=verbose)
    else:
        total_mask                  = source_mask
        total_mask                  = total_mask == 1 # Change to boolean, True where 1
        printme = 'No input mask file supplied, mask used is the program created source mask.'
        print_verbose_string(printme,verbose=verbose)

    return total_mask

def plot_contourOnImage(fitsfile,total_mask_bool,verbose=False):

    # Read in image
    image,h = fits.getdata(fitsfile,header=True)

    # Create header with wcs
    contour_fits        = fits.PrimaryHDU()
    contour_fits.data   = total_mask_bool.astype('int')
    contour_fits.header['CTYPE1']   = h['CTYPE1']
    contour_fits.header['CRPIX1']   = h['CRPIX1']
    contour_fits.header['CRVAL1']   = h['CRVAL1']
    contour_fits.header['CTYPE2']   = h['CTYPE2']
    contour_fits.header['CRPIX2']   = h['CRPIX2']
    contour_fits.header['CRVAL2']   = h['CRVAL2']
    contour_fits.header['CD1_1']    = h['CD1_1']
    contour_fits.header['CD1_2']    = h['CD1_2']
    contour_fits.header['CD2_1']    = h['CD2_1']
    contour_fits.header['CD2_2']    = h['CD2_2']
    try:
        contour_fits.header['EQUINOX']  = h['EQUINOX']
    except:
        print('IMPORTANT NOTE!!!! Equinox of input image assumed to be 2000.0')
        print('                   This is just for plotting checkim purposes')
        contour_fits.header['EQUINOX']  = 2000.0

    # Save contour_image to file, with fitsfile WCS
    total_mask_fitsWithWCS = './contour.fits'
    contour_fits.writeto(total_mask_fitsWithWCS)
    printme = f'SAVED  : {total_mask_fitsWithWCS}'
    print_verbose_string(printme,verbose=verbose)
    
    # Plot total_mask as contour on fits image
    fig = plt.figure(figsize=(48, 36))
    f2 = aplpy.FITSFigure(fitsfile,figure=fig)
    f2.ticks.hide()
    f2.tick_labels.hide_x()
    f2.tick_labels.hide_y()
    f2.axis_labels.hide()
    interval = ZScaleInterval()
    vmin,vmax = interval.get_limits(image)
    f2.show_grayscale(invert=True, stretch='linear', vmin=vmin, vmax=vmax)
    f2.show_contour(data=total_mask_fitsWithWCS,linewidths=3.0,colors='MediumPurple')
    cont_name = fitsfile.replace('.fits','_skymask_contour.png')
    f2.save(cont_name)
    print(f'SAVED  : {cont_name}')

    # Remove contour_image fits file
    clearit(total_mask_fitsWithWCS)
    printme = f'REMOVED: {total_mask_fitsWithWCS}'
    print_verbose_string(printme,verbose=verbose)

    return None

def calculate_stats_andPrint(image,mask,sky_counts,sky_counts_avg,pix_counts,verbose=False, quietmode=False):

    image_masked = np.ma.masked_array(image,mask=mask)
    
    sky = np.average(np.array(sky_counts))
    sky_error = np.std(np.array(sky_counts))
    sky_pixel_std = np.nanstd(image_masked)

    if not quietmode:
        printme = '------------------'
        print(printme)
        printme = 'PRINTING SKY STATS'
        print(printme)
        printme = '------------------'
        print(printme)

        printme = f'\n# SKY: Average of sky box/annuli medians      : {sky:.4f}'
        print(printme)
        printme = f'# SKY_ERROR: STD of sky box/annuli medians    : {sky_error:.4f}'
        print(printme)
        printme = f'# SKY_PIXEL_STD: STD of all non-masked pixels : {sky_pixel_std:.4f}\n'
        print(printme) 

    # Calculate other things for verbose printing
    if verbose:

        # Signpost what's about to be printed
        printme = '------------------------'
        print_verbose_string(printme,verbose=verbose)
        printme = 'PRINTING EXTRA SKY STATS'
        print_verbose_string(printme,verbose=verbose)
        printme = '------------------------'
        print_verbose_string(printme,verbose=verbose)

        # Calculate More Stats: stats on median of boxes/ annuli
        sky_avgOfAvg = np.average(np.array(sky_counts_avg))
        sky_medOfAvg = np.median(np.array(sky_counts_avg))
        sky_stdOfAvg = np.std(np.array(sky_counts_avg))

        # Calculate More Stats: stats on average of boxes/ annuli
        sky_avgOfMed = np.average(np.array(sky_counts)) # sky above
        sky_medOfMed = np.median(np.array(sky_counts))
        sky_stdOfMed = np.std(np.array(sky_counts)) # sky_error above

        # Print More Stats: stats on median of boxes/ annuli
        printme = f'# Average of sky box/annuli AVERAGES : {sky_avgOfAvg:.4f}'
        print_verbose_string(printme,verbose=verbose)
        printme = f'# Median  of sky box/annuli AVERAGES : {sky_medOfAvg:.4f}'
        print_verbose_string(printme,verbose=verbose)
        printme = f'# STD     of sky box/annuli AVERAGES : {sky_stdOfAvg:.4f}'
        print_verbose_string(printme,verbose=verbose)

        # Print More Stats: stats on average of boxes/ annuli
        printme = f'# Average of sky box/annuli MEDIANS  : {sky_avgOfMed:.4f}'
        print_verbose_string(printme,verbose=verbose)
        printme = f'# Median  of sky box/annuli MEDIANS  : {sky_medOfMed:.4f}'
        print_verbose_string(printme,verbose=verbose)
        printme = f'# STD     of sky box/annuli MEDIANS  : {sky_stdOfMed:.4f}'
        print_verbose_string(printme,verbose=verbose)

        # Print More Stats: Average number of non-masked pixels in boxes/annuli:
        printme = f'# Avg # unmasked pix in boxes/annuli : {np.average(pix_counts)}'
        print_verbose_string(printme,verbose=verbose)

        # Global Stats
        print(' ')
        printme = f'# Average of all non-masked pixels   : {np.nanmean(image_masked):.4f}'
        print_verbose_string(printme,verbose=verbose)
        printme = f'# Median  of all non-masked pixels   : {np.nanmedian(image_masked):.4f}'
        print_verbose_string(printme,verbose=verbose)
        printme = f'# Number of unmasked pixels          : {np.ma.sum(image_masked.mask)}'
        print_verbose_string(printme,verbose=verbose)
        printme = f'# Total number of pixels             : {np.size(image)}'
        print_verbose_string(printme,verbose=verbose)

    return sky,sky_error,sky_pixel_std


###########################################
# ======= Sky Calculation Functions =======
###########################################

def calculate_sky_box(fitsimage, boxsize_pix, n_iterations, input_mask_file=False, checkims=False, sextractorloc='/opt/local/bin/source-extractor',verbose=False, quietmode=False):
    '''Place n_iterations number of boxsize_pix sized boxes randomly in image with total_mask, 
    Output 
    sky         = average of median of boxes 
    sky_error   = std of median of boxes
    sky_pixel_std = std of all unmasked pixels.
    
    if verbose, output
    figure to show where boxes were placed'''
    
    sky_counts     = [] # median
    sky_counts_avg = []
    pix_counts     = []
    n_counter      = 0
    n_notfinite    = 0

    # Read in image and header
    image,h = fits.getdata(fitsimage, header=True)

    # Make source mask
    source_mask                    = create_sourceMask(fitsimage,sextractorloc=sextractorloc)

    # Combine with input_mask if input_mask_file supplied
    total_mask_bool = combine_masks(source_mask,input_mask_file,verbose=verbose)  

    # Plot total_mask as a contour on fits image
    if checkims:
        plot_contourOnImage(fitsimage,total_mask_bool,verbose=verbose)

    # Set boxes to be placed not too near edge
    xmin = 1.5*boxsize_pix
    ymin = 1.5*boxsize_pix
    xmax = float(h['NAXIS1'])-1.5*boxsize_pix
    ymax = float(h['NAXIS2'])-1.5*boxsize_pix

    # Start figure to plot up box locations
    if checkims:
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
    
    while n_counter <= n_iterations:
        
        # Choose a random spot
        row = np.random.randint(low=ymin,high=ymax)
        col = np.random.randint(low=xmin,high=xmax)

        # Make a box
        image_box = image[row-int(boxsize_pix/2):row+int(boxsize_pix/2)+1,col-int(boxsize_pix/2):col+int(boxsize_pix/2)+1]
        mask_box  = total_mask_bool[row-int(boxsize_pix/2):row+int(boxsize_pix/2)+1,col-int(boxsize_pix/2):col+int(boxsize_pix/2)+1]

        # Plot up location of box for display using show_contour
        if checkims:
            display_mask = np.zeros((xlen,ylen))
            display_mask[row-int(boxsize_pix/2):row+int(boxsize_pix/2)+1,col-int(boxsize_pix/2):col+int(boxsize_pix/2)+1] = 1.0
            CS = plt.contour(X, Y, display_mask,linewidths=1.0,alpha=0.1,colors='red')

        # Measure median counts in this masked box
        counts     = np.ma.median(np.ma.masked_array(image_box,mask=mask_box))
        counts_avg = np.ma.average(np.ma.masked_array(image_box,mask=mask_box))
        # Measure number of pixels not masked in this masked box
        no_pixels_notmasked = np.sum(mask_box)

        # Add average to sky_counts if finite
        # Also increment box count 
        # Else increment n_notfinite
        if np.isfinite(counts):
            sky_counts.append(counts)
            sky_counts_avg.append(counts_avg)
            pix_counts.append(no_pixels_notmasked)
            n_counter += 1
        else:
            n_notfinite += 1

    # Save figure to of box locations
    if checkims:
        outname = fitsimage.replace('.fits','_skyregionlocs.png')
        f1.save(outname)
        print(f'\nSAVED  : Box location plot saved: {outname}')

    printme=f'Number of attempts where average sky count in box was not finite: {n_notfinite}'
    print_verbose_string(printme,verbose=verbose)

    sky,sky_error,sky_pixel_std = calculate_stats_andPrint(image,total_mask_bool,sky_counts,sky_counts_avg,pix_counts,verbose=verbose,quietmode=quietmode)

    return sky, sky_error, sky_pixel_std # end calculate_sky_box

def calculate_sky_annuli(fitsimage,annulusparams,n_iterations,input_mask_file = False,checkims=False,verbose=False, quietmode=False):
    '''Place n_iterations number of elliptical annuli randomly in image with total_mask.
    Output 
    sky         = average of median of annuli 
    sky_error   = std of median of annuli
    sky_pixel_std = std of all unmasked pixels.'''
    exitmsg = 'ERROR: sky calculations using sky annuli not implemented yet.'
    sys.exit(exitmsg)

    # Read in image and header
    image,h = fits.getdata(fitsimage, header=True)

    # Make source mask
    source_mask                    = create_sourceMask(fitsimage)

    # Combine with input_mask if input_mask_file supplied
    total_mask_bool = combine_masks(source_mask,input_mask_file,verbose=verbose)  

    # Plot total_mask as a contour on fits image
    if checkims:
        plot_contourOnImage(fitsimage,total_mask_bool,verbose=verbose)

    # Calculate sky in input annulus
    xc1,yc1,a1,b1,ang1,xc2,yc2,a2,b2,ang2 = read_annulusparams(annulusparams)
    h    = fits.getheader(fitsimage)
    xlen = int(h['NAXIS2'])
    ylen = int(h['NAXIS1'])
    mask = make_annulus_mask(xlen,ylen,xc1,yc1,a1,b1,ang1,xc2,yc2,a2,b2,ang2)
    initial_annuli_mask_data = mask.copy()
    image_annuli        = copy.copy(image)
    image_annuli[mask]  = float('nan')     
    image_annuli[total_mask_bool] = float('nan')
    initial_annuli_name = 'annuli_input.fits'
    fits.writeto(initial_annuli_name,image_annuli)
    printme = 'SAVED  : {initial_annuli_name} (temporary)'
    print_verbose_string(printme)
    print('Average in input sky annulus is: ',np.nanmean(image_annuli))
    print('Median in input sky annulus is : ',np.nanmedian(image_annuli))
    print('Std in input sky annulus is    : ',np.nanstd(image_annuli))
    print('Number of finite non masked pixels in input sky annulus: ',np.sum(np.isfinite(image_annuli)))

    # Plonk some random annuli, calculate average of averages and std of averages
    # Vary xc,yc within width of annuli randomly (move xc2,yc2 by same amount)
    # AND vary a1 randomly while keeping a1-a2 constant, varations up to width of annuli
    annuli_thickness = abs(a1-a2)/2.

    # Start figure to plot up annuli locations
    fig = plt.figure(figsize=(48, 36))
    f1 = aplpy.FITSFigure(fitsimage,figure=fig)
    f1.ticks.hide()
    f1.tick_labels.hide_x()
    f1.tick_labels.hide_y()
    f1.axis_labels.hide()
    interval = ZScaleInterval()
    vmin,vmax = interval.get_limits(image)
    f1.show_grayscale(invert=True, stretch='linear', vmin=vmin, vmax=vmax)
    
    sky_counts      = []
    sky_counts_avg  = []
    pix_counts      = []
    n_counter       = 0
    n_notfinite     = 0
    xtomesh = np.arange(0, ylen, 1)
    ytomesh = np.arange(0, xlen, 1)
    X, Y    = np.meshgrid(xtomesh, ytomesh)
    while n_counter < n_iterations:

        # Choose X random values for xc,yc and a1
        xc_shift = np.random.randint(low=-annuli_thickness,high=annuli_thickness)
        yc_shift = np.random.randint(low=-annuli_thickness,high=annuli_thickness)
        a1_shift = np.random.randint(low=-annuli_thickness,high=annuli_thickness)
        new_xc1 = xc1+xc_shift
        new_xc2 = xc2+xc_shift
        new_yc1 = yc1+yc_shift
        new_yc2 = yc2+yc_shift
        new_a1  = a1+a1_shift
        new_a2  = a2+a1_shift
        new_b1  = (b1/a1)*(new_a1)
        new_b2  = (b2/a2)*(new_a2) 

        # Make mask for new annuli 
        mask = make_annulus_mask(xlen,ylen,
                                    new_xc1,new_yc1,new_a1,new_b1,ang1,
                                    new_xc2,new_yc2,new_a2,new_b2,ang2)
        image_annuli        = copy.copy(image)
        image_annuli[mask]  = float('nan')     
        image_annuli[total_mask_bool] = float('nan')

        # Plot up location annulus for display using show_contour
        CS = plt.contour(X, Y, mask,linewidths=1.0,alpha=0.1,colors='red')

        # Calculate average and number of pixels in average to array
        #counts = 3.*np.nanmedian(image_annuli) - 2.*np.nanmean(image_annuli)
        counts_avg  = np.nanmean(image_annuli)
        counts      = np.nanmedian(image_annuli)

        # Add average to sky_counts if finite
        # Also increment n_counter
        # Else increment n_notfinite
        if np.isfinite(counts):
            sky_counts.append(counts)
            sky_counts_avg.append(counts_avg)
            pix_counts.append(np.sum(np.isfinite(image_annuli)))
            n_counter += 1
        else:
            n_notfinite += 1
        
        # Increment counter
        n_counter += 1

    # Plot initial sky ellipse
    # Copy wcs to total_mask_name, and show initial ellipse contour
    CS = plt.contour(X, Y, initial_annuli_mask_data,linewidths=6.0,colors='green')

    # Save figure to of annuli locations
    outname = './skyregionlocs.png'
    f1.save(outname)
    printme=f'SAVED  : {outname}'
    print(printme)

    # Clear temporary files
    clearit(initial_annuli_name)
 
    # Print useful information
    print_verbose_string(f'Number of annuli placed randomly is: {n_counter}',verbose=verbose)
    print_verbose_string(f'#Number of attempts where average sky count in box/annuli was not finite: {str(n_notfinite)}',verbose=verbose)

    sky,sky_error,sky_pixel_std = calculate_stats_andPrint(image,total_mask_bool,sky_counts,sky_counts_avg,pix_counts,verbose=verbose, quietmode=quietmode)

    return sky, sky_error, sky_pixel_std # end calculate_sky_annuli

def calculate_sky_annuli_allover(fitsimage,annulusparams,n_iterations,input_mask_file = False,checkims=False,sextractorloc='/opt/local/bin/source-extractor',verbose=False, quietmode=False):
    '''Place n_iterations number of elliptical annuli randomly in image with total_mask.
    Annuli placed is allowed to move around more than the --annuli option.
    Output 
    sky         = average of median of annuli 
    sky_error   = std of median of annuli
    sky_pixel_std = std of all unmasked pixels.'''
    exitmsg = 'ERROR: sky calculations using sky annuli not implemented yet.'
    sys.exit(exitmsg)

    if verbose:
        print('#Number of attempts where average sky count in box/annuli was not finite: ',str(n_notfinite))

def calculate_sky_allunmaked(fitsimage,input_mask_file=False,checkims=False, sextractorloc='/opt/local/bin/source-extractor',verbose=False, quietmode=False):
    # Read in image and header
    image,h = fits.getdata(fitsimage, header=True)

    # Make source mask
    source_mask                    = create_sourceMask(fitsimage,sextractorloc=sextractorloc)

    # Combine with input_mask if input_mask_file supplied
    total_mask_bool = combine_masks(source_mask,input_mask_file,verbose=verbose)  

    # Plot total_mask as a contour on fits image
    if checkims:
        plot_contourOnImage(fitsimage,total_mask_bool,verbose=verbose)

    # Calculate sky stats
    masked_image  = np.ma.masked_array(image,mask=total_mask_bool)
    sky           = np.ma.median(masked_image)
    sky_error     = np.ma.std(masked_image)
    sky_pixel_std = sky_error

    if not quietmode:
        printme = '------------------'
        print(printme)
        printme = 'PRINTING SKY STATS'
        print(printme)
        printme = '------------------'
        print(printme)

        printme = f'\n# SKY: Average of sky box/annuli medians      : {sky:.4f}'
        print(printme)
        printme = f'# SKY_ERROR: STD of sky box/annuli medians    : {sky_error:.4f}'
        print(printme)
        printme = f'# SKY_PIXEL_STD: STD of all non-masked pixels : {sky_pixel_std:.4f}\n'
        print(printme) 

    return sky, sky_error, sky_pixel_std

###############################
# ======= Main function =======
###############################

def calculate_skyStats( fitsimage, 
                        sextractorloc        = '/opt/local/bin/source-extractor',
                        place_boxes          = False, 
                        place_annuli         = False, 
                        place_annuli_allover = False,
                        n_iterations         = 100,
                        input_mask_file      = False,
                        checkims             = False,
                        quietmode = False, verbose = False, debugmode = False):
    '''calculate_skyStats.py -- Mask stars aggressively, determine sky background stats in remaining pixels. Stats are
    sky:            Average of many sky box/annuli medians (which are each a sample of the sky value) 
    sky_error:      Standard deviation of many sky box/annuli medians (which are each a sample of the sky value) 
    sky_pixel_std:  Standard deviation of all non-masked pixels  

    Three options for sampling the sky are:
    (1) input --box s -n n: calculate sky in n boxes of size sxs randomly placed around the image
    (2) input --annulus 'xc1,yc1,a1,b1,ang1,xc2,yc2,a2,b2,ang2' -n n: calculate the sky in n elliptical annuli placed randomly. xc1,yc1,a1,b1,ang1 specifies the inner ellipse, and xc2,yc2,a2,b2,ang2 specifies the outer ellipse of the annuli. 
    (3) input --annulusallover 'xc1,yc1,a1,b1,ang1,xc2,yc2,a2,b2,ang2' -n n: as above, but letting it move around a little more than --annulus option
    If no option is set, then the average, median and standard deviation of the unmasked pixels in the image is calculated. 

    Option to input additional mask (mask should be zeros and ones; zero pixels are not masked, one pixels are masked). Important to input a mask that masks out the galaxy and its extended low surface brightness features!

    Output are XXX
    '''

    print_debug_string(f'sextractor loc being used by calculate_skyStats is: {sextractorloc}',debugmode=debugmode)


    # Check which option(s) were set. If more than 1, ask user to re-set, if zero, set place_no_shapes to True.
    place_no_shapes = False
    options_true = [place_boxes, place_annuli, place_annuli_allover]
    options_true = [bool(x) for x in options_true]
    if sum(options_true) > 1:
        error_message = 'ERROR: Please only select one of --boxes, --annulus, and --annulusallover'
        sys.exit(error_message)
    if sum(options_true) == 0:
        place_no_shapes = True
        printme='No option for placing shapes for calculating sky was set, so the avg, median, std of all unmasked pixels in the image will be calculated.'
        print_verbose_string(printme,verbose=verbose)
        
    # ===== --box =====
    # Place n_iterations random boxsize_pix pixel sized boxes on total_masked image
    # Calculate median in each box, average of these medians and
    # Standard deviation of these medians
    if place_boxes:
        boxsize_pix     = int(place_boxes)
        sky, sky_error, sky_pixel_std = calculate_sky_box(fitsimage,boxsize_pix, n_iterations, input_mask_file = input_mask_file, checkims=checkims, sextractorloc=sextractorloc, verbose=verbose, quietmode=quietmode)

    # ===== --annulus =====
    # Place X random X pixel elliptical annuli 
    # Calculate median in each annuli, average of these medians and
    # Standard deviation of these medians
    if place_annuli:
        annulusparams = place_annuli
        sky, sky_error, sky_pixel_std = calculate_sky_annuli(fitsimage,annulusparams,n_iterations,input_mask_file=input_mask_file,checkims=checkims,sextractorloc=sextractorloc, verbose=verbose, quietmode=quietmode)

    # ===== --annulusallover =====
    # Place X random X pixel elliptical annuli 
    # Calculate median in each annuli, average of these medians and
    # Standard deviation of these medians
    if place_annuli_allover:
        annulusparams = place_annuli_allover
        sky, sky_error, sky_pixel_std = calculate_sky_annuli_allover(fitsimage,annulusparams,n_iterations,input_mask_file=input_mask_file,checkims=checkims,sextractorloc=sextractorloc,verbose=verbose, quietmode=quietmode)

    # ===== no option set, just calculate stats for all unmasked pixels =====
    if (not  place_boxes) and (not place_annuli) and (not place_annuli_allover):        
        sky, sky_error, sky_pixel_std = calculate_sky_allunmaked(fitsimage,input_mask_file=input_mask_file,checkims=checkims,sextractorloc=sextractorloc,verbose=verbose, quietmode=quietmode)

    return sky, sky_error, sky_pixel_std



###########################################################
###################### Start of main ######################
###########################################################
if __name__=='__main__':

    # Import arguments
    arguments = docopt.docopt(__doc__)
    
    fitsimage            = arguments['<fitsimage>']
    quietmode            = arguments['--quietmode']
    verbose              = arguments['--verbose']
    debugmode            = arguments['--debug']
    sextractorloc        = arguments['--sextractorloc']
    place_boxes          = arguments['--box']
    place_annuli         = arguments['--annulus']    
    place_annuli_allover = arguments['--annulusallover']
    n_iterations         = int(arguments['--niterations'])
    input_mask_file      = arguments['--mask']
    checkims             = arguments['--checkims']

    if debugmode:
        print(arguments)

    calculate_skyStats(  fitsimage, 
                            sextractorloc        = sextractorloc,
                            place_boxes          = place_boxes, 
                            place_annuli         = place_annuli, 
                            place_annuli_allover = place_annuli_allover,
                            n_iterations         = n_iterations,
                            input_mask_file      = input_mask_file,
                            checkims             = checkims,
                            quietmode=quietmode, verbose=verbose, debugmode=debugmode)

