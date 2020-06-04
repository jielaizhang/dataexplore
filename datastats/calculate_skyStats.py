#!/usr/bin/python

"""calculate_skystats.py -- Mask stars aggressively, together with input mask, determine sky background stats in remaining pixels. Stats are
    Average
    Median
    Standard deviation of x by x pixels's average values 
        (y boxes are placed randomly)

Note that for using aplpy to plot location of randomly placed sky areas, the vmin and vmax for fits image is hard coded in right now. So if you get black png files, that's why. For randomly placed boxes, take in vmin, vmax command line options; but not for annuli! 

Required input 
    fitsimage - image for which background stats are desired.

Usage:
    calculate_skystats.py [-h] [-v] [-b STRING] [-a STRING] [--annulusallover STRING] [-n INT] [-m FILE] [-s SEXLOC] [--vmin FLOAT] [--vmax FLOAT] <fitsimage>

Options:
    -h, --help                      Print this screen.
    -v, --verbose                   Print extra information [default: False]
    -b STRING, --box STRING         Input box parameters: size of box (pixels), number of random boxes to be placed. E.g. 8,1000.
    -a STRING, --annulus STRING     Select annulus with params 'xc1,yc1,a1,b1,ang1,xc2,yc2,a2,b2,ang2', angle in degrees, counter clockwise rotation; place random annuli around galaxy.
    --annulusallover STRING         Select annulus with params 'xc1,yc1,a1,b1,ang1,xc2,yc2,a2,b2,ang2', angle in degrees, counter clockwise rotation; place annulus around galaxy, but let it move around a little more than above option.
    -n INT, --niterations INT       Input number of random annuli to be placed. (not used for boxes). [default: 100]
    -m FILE, --mask FILE            Required. Input mask to be combined with grown sextractor mask. Important to input a mask that masks out the galaxy and its extended low surface brightness features!
    -s SEXLOC, --sex SEXLOC         SExtractor location [default: /opt/local/bin/sex]   
    --vmin FLOAT                    For plotting up box/ annuli locations [default: 4.25]
    --vmax FLOAT                    For plotting up box/ annuli locations [default: 4.32]

Example:
    python calculate_skystats.py -v fitsimage.fits
"""

import docopt
import numpy as np
import astropy.io.fits as fits
from scipy import ndimage
import subprocess
import os, sys, copy

import matplotlib.pyplot as plt
import aplpy

sextractor_params = """NUMBER
FLUX_AUTO
FLUXERR_AUTO
FLUX_APER
FLUXERR_APER
X_IMAGE
Y_IMAGE
X_WORLD
Y_WORLD
FLUX_RADIUS
FLAGS
CLASS_STAR
BACKGROUND
ELLIPTICITY
FWHM_IMAGE
"""

sextractor_config = """
    ANALYSIS_THRESH 3
        BACK_FILTERSIZE 3
        BACKPHOTO_TYPE LOCAL
        BACK_SIZE 32
        CATALOG_NAME test.cat
        CATALOG_TYPE ASCII_HEAD
        CHECKIMAGE_TYPE SEGMENTATION
        CHECKIMAGE_NAME {check_name}
        CLEAN Y
        CLEAN_PARAM 1.
        DEBLEND_MINCONT 0.001
        DEBLEND_NTHRESH 32
        DETECT_MINAREA 5
        DETECT_THRESH 3
        DETECT_TYPE CCD
        FILTER Y
        FILTER_NAME {filter_name}
        FLAG_IMAGE flag.fits
        GAIN 1.0
        MAG_GAMMA 4.
        MAG_ZEROPOINT 0.0
        MASK_TYPE CORRECT
        MEMORY_BUFSIZE 1024
        MEMORY_OBJSTACK 3000
        MEMORY_PIXSTACK 300000
        PARAMETERS_NAME {parameters_name}
        PHOT_APERTURES 5
        PHOT_AUTOPARAMS 2.5, 3.5
        PIXEL_SCALE 2.85
        SATUR_LEVEL 50000.
        SEEING_FWHM 2.5
        STARNNW_NAME {starnnw_name}
        VERBOSE_TYPE {verbose_type}
"""

default_conv = """CONV NORM
# 3x3 ``all-ground'' convolution mask with FWHM = 2 pixels.
1 2 1
2 4 2
1 2 1
"""

default_nnw = """NNW
# Neural Network Weights for the SExtractor star/galaxy classifier (V1.3)
# inputs:       9 for profile parameters + 1 for seeing.
# outputs:      ``Stellarity index'' (0.0 to 1.0)
# Seeing FWHM range: from 0.025 to 5.5'' (images must have 1.5 < FWHM < 5 pixels)
# Optimized for Moffat profiles with 2<= beta <= 4.

 3 10 10  1

-1.56604e+00 -2.48265e+00 -1.44564e+00 -1.24675e+00 -9.44913e-01 -5.22453e-01  4.61342e-02  8.31957e-01  2.15505e+00  2.64769e-01
 3.03477e+00  2.69561e+00  3.16188e+00  3.34497e+00  3.51885e+00  3.65570e+00  3.74856e+00  3.84541e+00  4.22811e+00  3.27734e+00

-3.22480e-01 -2.12804e+00  6.50750e-01 -1.11242e+00 -1.40683e+00 -1.55944e+00 -1.84558e+00 -1.18946e-01  5.52395e-01 -4.36564e-01 -5.30052e+00
 4.62594e-01 -3.29127e+00  1.10950e+00 -6.01857e-01  1.29492e-01  1.42290e+00  2.90741e+00  2.44058e+00 -9.19118e-01  8.42851e-01 -4.69824e+00
-2.57424e+00  8.96469e-01  8.34775e-01  2.18845e+00  2.46526e+00  8.60878e-02 -6.88080e-01 -1.33623e-02  9.30403e-02  1.64942e+00 -1.01231e+00
 4.81041e+00  1.53747e+00 -1.12216e+00 -3.16008e+00 -1.67404e+00 -1.75767e+00 -1.29310e+00  5.59549e-01  8.08468e-01 -1.01592e-02 -7.54052e+00
 1.01933e+01 -2.09484e+01 -1.07426e+00  9.87912e-01  6.05210e-01 -6.04535e-02 -5.87826e-01 -7.94117e-01 -4.89190e-01 -8.12710e-02 -2.07067e+01
-5.31793e+00  7.94240e+00 -4.64165e+00 -4.37436e+00 -1.55417e+00  7.54368e-01  1.09608e+00  1.45967e+00  1.62946e+00 -1.01301e+00  1.13514e-01
 2.20336e-01  1.70056e+00 -5.20105e-01 -4.28330e-01  1.57258e-03 -3.36502e-01 -8.18568e-02 -7.16163e+00  8.23195e+00 -1.71561e-02 -1.13749e+01
 3.75075e+00  7.25399e+00 -1.75325e+00 -2.68814e+00 -3.71128e+00 -4.62933e+00 -2.13747e+00 -1.89186e-01  1.29122e+00 -7.49380e-01  6.71712e-01
-8.41923e-01  4.64997e+00  5.65808e-01 -3.08277e-01 -1.01687e+00  1.73127e-01 -8.92130e-01  1.89044e+00 -2.75543e-01 -7.72828e-01  5.36745e-01
-3.65598e+00  7.56997e+00 -3.76373e+00 -1.74542e+00 -1.37540e-01 -5.55400e-01 -1.59195e-01  1.27910e-01  1.91906e+00  1.42119e+00 -4.35502e+00

-1.70059e+00 -3.65695e+00  1.22367e+00 -5.74367e-01 -3.29571e+00  2.46316e+00  5.22353e+00  2.42038e+00  1.22919e+00 -9.22250e-01 -2.32028e+00


 0.00000e+00
 1.00000e+00
"""

def run_SExtractor(image_name):

    'Create temporary directory'
    if not os.path.exists('tmpSkystats'):
        os.makedirs('tmpSkystats')    
    else:
        print('./tmpSkystats directory existed already')

    'Names of required config files'
    sextractor_config_name = './tmpSkystats/scamp.sex'
    params_name = './tmpSkystats/scamp.param'
    conv_name = './tmpSkystats/default.conv'
    nnw_name = './tmpSkystats/default.nnw'

    catalog_name = image_name.split('.fits')[0]+'_bkg.cat'
    check_name = image_name.split('.fits')[0]+'_bkg_segmap.fits'

    if verbose:
        verbose_type = 'NORMAL'
    else:
        verbose_type = 'QUIET'

    'Stick content in config files'
    configs = zip([sextractor_config_name,params_name,conv_name,nnw_name],[sextractor_config,sextractor_params,default_conv,default_nnw])

    for fname,fcontent in configs:
        fout = open(fname,'w')

        if 'scamp.sex' in fname:
            fout.write(fcontent.format(filter_name=conv_name,
            parameters_name=params_name,starnnw_name=nnw_name,
            verbose_type=verbose_type,check_name=check_name))

        else:
            fout.write(fcontent)

        fout.close()

    if verbose:
        print('SExtracting...')


    'SExtractor command'
    command = sexloc + ' -c {config} -CATALOG_NAME {catalog} {image}'.format(config=sextractor_config_name,catalog=catalog_name,image=image_name)

    if verbose:
        print('Running this command:')
        print(command+'\n')

    'Run SExtractor'
    subprocess.call(command,shell=True)

    'Clear unnecessary files'
    for fname in [sextractor_config_name,params_name,conv_name,nnw_name]:
        clearit(fname)
    
    'Remove temp directory if its not empty'
    try:
        os.rmdir('tmpSkystats')
    except OSError as ex:
        if ex.errno == errno.ENOTEMPTY:
            print("directory not empty")

    return check_name

def clearit(fname):
    if os.path.isfile(fname):
        os.remove(fname)
    return None

def writeFITS(im,saveAs,header=None):
    if header != None:
        hdu = fits.PrimaryHDU(data=im,header=header)
    else:
        hdu = fits.PrimaryHDU(data=im)

    hdulist = fits.HDUList([hdu])
    hdulist.writeto(saveAs,overwrite=True)
    hdulist.close()
    return None

def calculate_sky_box(fitsimage,image,total_mask,boxsize_pix,nboxes,vmin=4.25,vmax=4.32):
    '''Place nboxes boxsize_pix sized boxes randomly in image with total_mask, calculate average in each box and standard deviation of averages'''
    
    sky_counts  = []
    pix_counts  = []
    n_counter   = 0
    n_notfinite = 0

    # Read in image size and set boxes 
    # to be placed not too near edge
    h     = fits.getheader(fitsimage)
    xmin = 1.5*boxsize_pix
    ymin = 1.5*boxsize_pix
    xmax = float(h['NAXIS1'])-1.5*boxsize_pix
    ymax = float(h['NAXIS2'])-1.5*boxsize_pix

    # Start figure to plot up box locations
    fig = plt.figure(figsize=(48, 36))
    f1 = aplpy.FITSFigure(fitsimage,figure=fig)
    f1.set_tick_labels_font(size='xx-small')
    f1.ticks.hide()
    f1.tick_labels.hide_x()
    f1.tick_labels.hide_y()
    f1.axis_labels.hide()
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
        CS = plt.contour(X, Y, display_mask,thickness=1.0,alpha=0.1,colors='red')

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
    
    return sky_counts, pix_counts, n_notfinite, h

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

def calculate_sky_annuli(fitsimage,image,total_mask,annulusparams,n_iterations):
    '''Save sky count averages in n_iteration annuli, also plot up where random n_iterations of annuli were placed on fits image.'''
    
    # Calculate sky in input annulus
    xc1,yc1,a1,b1,ang1,xc2,yc2,a2,b2,ang2 = read_annulusparams(annulusparams)
    h    = fits.getheader(fitsimage)
    xlen = int(h['NAXIS2'])
    ylen = int(h['NAXIS1'])
    mask = make_annulus_mask(xlen,ylen,xc1,yc1,a1,b1,ang1,xc2,yc2,a2,b2,ang2)
    initial_annuli_mask_data = mask.copy()
    image_annuli        = copy.copy(image)
    image_annuli[mask]  = float('nan')     
    image_annuli[total_mask] = float('nan')
    initial_annuli_name = 'annuli_input.fits'
    writeFITS(image_annuli,initial_annuli_name)
    print(' ')
    print('***OUTPUT: Sky calculation annulus saved here: ',initial_annuli_name)
    print(' ')
    print('Average in input sky annulus is: ',np.nanmean(image_annuli))
    print('Median in input sky annulus is: ',np.nanmedian(image_annuli))
    print('Std in input sky annulus is: ',np.nanstd(image_annuli))
    print('Number of finite non masked pixels in input sky annulus: ',np.sum(np.isfinite(image_annuli)))

    # Plonk some random annuli, calculate average of averages and std of averages
    # Vary xc,yc within width of annuli randomly (move xc2,yc2 by same amount)
    # AND vary a1 randomly while keeping a1-a2 constant, varations up to width of annuli
    annuli_thickness = abs(a1-a2)/2.

    # Start figure to plot up annuli locations
    fig = plt.figure(figsize=(48, 36))
    f1 = aplpy.FITSFigure(fitsimage,figure=fig)
    f1.set_tick_labels_font(size='xx-small')
    f1.ticks.hide()
    f1.tick_labels.hide_x()
    f1.tick_labels.hide_y()
    f1.axis_labels.hide()
    f1.show_grayscale(invert=True, stretch='linear', vmin=4.25, vmax=4.32)
    # for g-band ngc 2841: vmin=2.38, vmax=2.42
    
    sky_counts  = []
    pix_counts  = []
    n_counter   = 0
    n_notfinite = 0
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
        image_annuli[total_mask] = float('nan')

        # Plot up location annulus for display using show_contour
        CS = plt.contour(X, Y, mask,thickness=1.0,alpha=0.1,colors='red')

        # Calculate average and number of pixels in average to array
        #counts = 3.*np.nanmedian(image_annuli) - 2.*np.nanmean(image_annuli)
        counts = np.nanmean(image_annuli)

        # Add average to sky_counts if finite
        # Also increment n_counter
        # Else increment n_notfinite
        if np.isfinite(counts):
            sky_counts.append(counts)
            pix_counts.append(np.sum(np.isfinite(image_annuli)))
            n_counter += 1
        else:
            n_notfinite += 1
        
        # Increment counter
        n_counter += 1

    # Plot initial sky ellipse
    # Copy wcs to total_mask_name, and show initial ellipse contour
    CS = plt.contour(X, Y, initial_annuli_mask_data,thickness=6.0,colors='green')

    # Save figure to of annuli locations
    outname = './skyregionlocs.png'
    f1.save(outname)
    print(' ')
    print('***OUTPUT: Annuli location plot saved here: ',outname)
    
    if verbose:
        print('Number of annuli placed randomly is: ',n_counter)

    return sky_counts, pix_counts, n_notfinite, h


def copy_wcs(fits_withwcs,fits_withoutwcs):
    h = fits.getheader(fits_withwcs)
    f           = fits.open(fits_withoutwcs)
    newf        = fits.PrimaryHDU()
    newf.header = f[0].header
    newf.data   = f[0].data
    newf.header['CTYPE1']   = h['CTYPE1']
    newf.header['CRPIX1']   = h['CRPIX1']
    newf.header['CRVAL1']   = h['CRVAL1']
    newf.header['CTYPE2']   = h['CTYPE2']
    newf.header['CRPIX2']   = h['CRPIX2']
    newf.header['CRVAL2']   = h['CRVAL2']
    newf.header['CD1_1']    = h['CD1_1']
    newf.header['CD1_2']    = h['CD1_2']
    newf.header['CD2_1']    = h['CD2_1']
    newf.header['CD2_2']    = h['CD2_2']
    #newf.header['RADECSYS'] = h['RADECSYS']
    newf.header['EQUINOX']  = h['EQUINOX']
    saveloc     = fits_withoutwcs.split('.')[0]+'_wcs.fits'
    newf.writeto(saveloc, overwrite=True)     
    return saveloc

def calculate_sky_annuli_alloverim(fitsimage,image,total_mask,annulusparams,n_iterations):
    '''Save sky count averages in n_iteration annuli, also plot up where random n_iterations of annuli were placed on fits image.'''
    
    # Calculate sky in input annulus
    xc1,yc1,a1,b1,ang1,xc2,yc2,a2,b2,ang2 = read_annulusparams(annulusparams)
    h    = fits.getheader(fitsimage)
    xlen = int(h['NAXIS2'])
    ylen = int(h['NAXIS1'])
    mask = make_annulus_mask(xlen,ylen,xc1,yc1,a1,b1,ang1,xc2,yc2,a2,b2,ang2)
    initial_annuli_mask_data = mask.copy()
    image_annuli        = copy.copy(image)
    image_annuli[mask]  = float('nan')     
    image_annuli[total_mask] = float('nan')
    initial_annuli_name = 'annuli_input.fits'
    writeFITS(image_annuli,initial_annuli_name)
    print(' ')
    print('***OUTPUT: Sky calculation annulus saved here: ',initial_annuli_name)
    print(' ')
    print('Average in input sky annulus is: ',np.nanmean(image_annuli))
    print('Median in input sky annulus is: ',np.nanmedian(image_annuli))
    print('Std in input sky annulus is: ',np.nanstd(image_annuli))
    print('Number of finite non masked pixels in input sky annulus: ',np.sum(np.isfinite(image_annuli)))

    # Plonk some random annuli, calculate average of averages and std of averages
    # Vary xc,yc within width of annuli randomly (move xc2,yc2 by same amount)
    # AND vary a1 randomly while keeping a1-a2 constant, varations up to width of annuli
    annuli_thickness = abs(a1-a2)/2.

    # Start figure to plot up annuli locations
    fig = plt.figure(figsize=(48, 36))
    f1 = aplpy.FITSFigure(fitsimage,figure=fig)
    f1.set_tick_labels_font(size='xx-small')
    f1.ticks.hide()
    f1.tick_labels.hide_x()
    f1.tick_labels.hide_y()
    f1.axis_labels.hide()
    f1.show_grayscale(invert=True, stretch='linear', vmin=4.25, vmax=4.32)
    # g-band ngc 2841 vmin=2.38, vmax=2.42

    sky_counts  = []
    pix_counts  = []
    n_counter   = 0
    n_notfinite = 0
    xtomesh = np.arange(0, ylen, 1)
    ytomesh = np.arange(0, xlen, 1)
    X, Y    = np.meshgrid(xtomesh, ytomesh)
    while n_counter < n_iterations:

        # Choose X random values for xc,yc and a1
        xc_shift = np.random.randint(low=-a1/3.,high=a1/3.)
        yc_shift = np.random.randint(low=-a1/3.,high=a1/3.)
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
        image_annuli[total_mask] = float('nan')

        # Plot up location annulus for display using show_contour
        CS = plt.contour(X, Y, mask,thickness=1.0,alpha=0.1,colors='red')

        # Calculate average and number of pixels in average to array
        #counts = 3.*np.nanmedian(image_annuli) - 2.*np.nanmean(image_annuli)
        counts = np.nanmean(image_annuli)

        # Add average to sky_counts if finite
        # Also increment n_counter
        # Else increment n_notfinite
        if np.isfinite(counts):
            sky_counts.append(counts)
            pix_counts.append(np.sum(np.isfinite(image_annuli)))
            n_counter += 1
        else:
            n_notfinite += 1
        
        # Increment counter
        n_counter += 1

    # Plot initial sky ellipse
    # Copy wcs to total_mask_name, and show initial ellipse contour
    CS = plt.contour(X, Y, initial_annuli_mask_data,thickness=6.0,colors='green')

    # Save figure to of annuli locations
    outname = './skyregionlocs.png'
    f1.save(outname)
    print(' ')
    print('***OUTPUT: Annuli location plot saved here: ',outname)
    
    if verbose:
        print('Number of annuli placed randomly is: ',n_counter)

    return sky_counts, pix_counts, n_notfinite, h

###################### Start of main ######################
if __name__=='__main__':

    # Import arguments
    arguments = docopt.docopt(__doc__)
    
    fitsimage           = arguments['<fitsimage>']
    verbose             = arguments['--verbose']
    place_boxes         = arguments['--box']
    place_annuli        = arguments['--annulus']    
    place_annuli_allover = arguments['--annulusallover']
    n_iterations        = int(arguments['--niterations'])
    maskfile            = arguments['--mask']
    sexloc              = arguments['--sex']

    vmax                = float(arguments['--vmax'])
    vmin                = float(arguments['--vmin'])

    if verbose:
        print('Measuring background sky stats in ',fitsimage)

    # Read in image
    image = fits.getdata(fitsimage)
    # Read in mask
    if maskfile:
        mask_input = fits.getdata(maskfile)
    else:
        sys.exit('sys.exit Error: Please input a mask file to mask out the central galaxy/ other galaxies in the input field')

    # Run SExtractor to produce segmentation map
    segname = run_SExtractor(fitsimage)
    print('***OUTPUT: SExtractor seg map saved here: ',segname)
    segmap  = fits.getdata(segname)
    
    # Grow mask so that sources are definitely covered
    grownsegmap = ndimage.grey_dilation(segmap,size=(15,15))
    
    # Save mask for reference
    grownsegmapname = segname.split('segmap.fits')[0]+'grownsegmap.fits'
    writeFITS(grownsegmap,grownsegmapname)
    print(' ')
    print('***OUTPUT: Grown seg map saved here: ',grownsegmapname)

    # Create total mask file
    seg_mask = copy.copy(grownsegmap)
    seg_mask[np.where(seg_mask>0)] = 1
    total_mask = (seg_mask | mask_input.astype('int'))
    total_mask = total_mask.astype('Bool')
    
    # Save total mask for reference
    total_mask_name = grownsegmapname = segname.split('segmap.fits')[0]+'total_mask.fits'
    writeFITS(total_mask.astype('int'),total_mask_name)
    print(' ')
    print('***OUTPUT: Mask used to block out non-sky pixels saved here: ',total_mask_name)

    # Copy wcs to total_mask_name
    total_mask_name_withwcs = copy_wcs(fitsimage,total_mask_name)

    # Plot up total_mask_name as contour on fits image
    fig = plt.figure(figsize=(48, 36))
    f2 = aplpy.FITSFigure(fitsimage,figure=fig)
    f2.set_tick_labels_font(size='xx-small')
    f2.ticks.hide()
    f2.tick_labels.hide_x()
    f2.tick_labels.hide_y()
    f2.axis_labels.hide()
    f2.show_grayscale(invert=True, stretch='linear', vmin=4.25, vmax=4.32)
    # g-band ngc 2841: vmin=2.38, vmax=2.42
    f2.show_contour(data=total_mask_name_withwcs,thickness=3.0,colors='MediumPurple')
    cont_name = './mask_contour.png'
    f2.save(cont_name)
    print(' ')
    print('***OUTPUT: Plot showing masked areas saved here: ',cont_name)

# ====== Boxes ======
    # Place nboxes random boxsize_pix pixel sized boxes on total_masked image
    # Calculate average in each box and
    # Standard deviation of these averages
    if place_boxes:
        boxsize_pix = int(place_boxes.split(',')[0])
        nboxes      = int(place_boxes.split(',')[1])
        sky_counts, pix_counts, n_notfinite, h = calculate_sky_box(fitsimage,image,total_mask,boxsize_pix,nboxes)

        if verbose:
            print('#Number of attempts where average sky count in box/annuli was not finite: ',str(n_notfinite))

# ====== Annuli ======
    # Place X random X pixel elliptical annuli 
    # Calculate average in each annuli and
    # Standard deviation of these averages
    if place_annuli:
        annulusparams = place_annuli
        sky_counts, pix_counts, n_notfinite, h = calculate_sky_annuli(fitsimage,image,total_mask,annulusparams,n_iterations)

        if verbose:
            print('#Number of attempts where average sky count in box/annuli was not finite: ',str(n_notfinite))

# ====== Annuli anywhere on image ======
    # Place X random X pixel elliptical annuli 
    # Calculate average in each annuli and
    # Standard deviation of these averages
    if place_annuli_allover:
        annulusparams = place_annuli_allover
        sky_counts, pix_counts, n_notfinite, h = calculate_sky_annuli_alloverim(fitsimage,image,total_mask,annulusparams,n_iterations)

        if verbose:
            print('#Number of attempts where average sky count in box/annuli was not finite: ',str(n_notfinite))

# ====== Calculate sky stats ======
    # Create Masked away to calculate other sky stats
    image_masked = np.ma.masked_where(grownsegmap>0,image)

    # Print out results
    print(' ')
    if place_annuli or place_boxes or place_annuli_allover:
        print(' ')
        print('#STD of sky box averages: ',np.std(np.array(sky_counts)))
        print('#Average of sky box averages: ',np.average(np.array(sky_counts)))
        print('#Median of sky box averages: ',np.median(np.array(sky_counts)))
        print('#Average number of pixels averaged per box: ',np.ma.average(pix_counts))
    print(' ')
    print('#Average of all non masked pix: ',np.ma.average(image_masked))
    print('#Median of all non masked pix: ',np.ma.median(image_masked))
    print('#STD of all non masked pix: ',np.ma.std(image_masked))
    print('#Number of unmasked pixels: ',np.ma.sum(image_masked.mask))
    print('#Total number of pixels: ',float(h['NAXIS1'])*float(h['NAXIS2']))


    
        

