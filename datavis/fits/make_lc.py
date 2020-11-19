#!/usr/bin/env python

""" make_lc.py -- input RA DEC file and string of fits files, make lc for each RA and DEC. 
 
Usage: make_lc [-h] [-v] [--debug] [-q] [-o SAVEDIR] <RADECfile> <fitsfiles>...
!!! Photometry apertures currently hard coded in. 

Arguments:
    RADECfile (string)
    fitsfiles (string)

Options:
    -h, --help                              Show this screen
    -v, --verbose                           Show extra information [default: False]   
    -q, --quietmode                         Minimize print to screen. This is useful when this function is called in another function. [default: False]  
    --debug                                 Output more for debugging [default: False]
    -o SAVEDIR, --out SAVEDIR               Saved output as [default: ./]

Examples:
"""
import docopt
import astropy.io.fits as fits
import numpy as np
import sys, os, ntpath
from pathlib import Path
from datetime import datetime
from astropy import units as u
from astropy.coordinates import SkyCoord
from photutils import SkyCircularAperture, SkyCircularAnnulus, CircularAperture
from photutils import aperture_photometry
from astropy.wcs import WCS
from astropy.stats import sigma_clipped_stats
import astropy.io.ascii as ascii
import matplotlib.pyplot as plt

import copy

__author__      = "Jielai Zhang"
__license__     = "MIT"
__version__     = "1.0.1"
__date__        = "2020-11-19"
__maintainer__  = "Jielai Zhang"
__email__       = "zhang.jielai@gmail.com"

#########################################
# ======= House Keeping Functions =======
#########################################
'''These functions help organise or delete files'''
def clearit(fnames,debugmode=False):
    for fname in fnames:
        if os.path.isfile(fname):
            os.remove(fname)
            printme = f'Temporary file deleted: {fname}.'
            print_debug_string(printme,debugmode=debugmode)
    return None

#########################################
# ======= Other Functions =======
#########################################
# Given RA, DEC, 
# get x_centre,y_centre,aperture_sum as an array, aperture area as an array, background_median as an array
def get_photometry_ingredients(fits_file,RAs,DECs,aperture_radii,annuli_inner_radii,annuli_outer_radii):
    
    # Get data, wcs
    data,h = fits.getdata(fits_file,header=True)
    w      = WCS(h)    
    
    # Define positions where aperture photometry will be performed
    positions = SkyCoord(ra=RAs, dec=DECs,unit='deg')
    
    # Create apertures at positions for photometry
    aperture_sky = [SkyCircularAperture(positions, r=r) for r in radii]
    
    # Turn apertures defined in sky coordinates into pixel coordinates
    aperture_pixel = [a_sky.to_pixel(wcs=w) for a_sky in aperture_sky]
    
    # Perform aperture photometry
    phot = aperture_photometry(data, aperture_pixel)
    
    # Get out Array of aperture sums and areas
    aperture_sums = []
    aperture_areas = []
    for ii in range(len(aperture_radii)):
        aperture_sums.append(phot['aperture_sum_'+str(ii)])
        aperture_areas.append(aperture_pixel[ii].area)

    background_medians = []
    # Get background medians
    for r_inner, r_outer in zip(annuli_inner_radii,annuli_outer_radii):
        
        # Define annulus aperture
        annulus_aperture_sky = SkyCircularAnnulus(positions, r_in=r_inner, r_out=r_outer)
        
        # Turn annulus defined in sky coordinates to pixel coordinates
        annulus_masks_pixel = annulus_aperture_sky.to_pixel(wcs=w)
        
        # Define annuli masks to get annuli values from data:
        annulus_masks = annulus_masks_pixel.to_mask(method='center')
        
        # Get median background value around each apperture position
        bkg_median = []
        for mask in annulus_masks:
            annulus_data = mask.multiply(data)
            annulus_data_1d = annulus_data[mask.data > 0]
            _, median_sigclip, _ = sigma_clipped_stats(annulus_data_1d)
            bkg_median.append(median_sigclip)
        bkg_median = np.array(bkg_median)
        background_medians.append(bkg_median)
        
    return np.array(aperture_sums), np.array(aperture_areas), np.array(background_medians)


# Given aperture sum array, aperture area array, background median array, 
# get photometry for each aperture and background combination
def do_photometry_cals(aperture_sums,aperture_areas,bgrs_per_pixel,
                       radii,annuli_inner_radii,annuli_outer_radii):
    
    # Set up output size
    phot = np.zeros([len(aperture_areas),len(bgrs_per_pixel),np.shape(aperture_sums)[1]])
    phot_explainer = np.array([[' '*20]*len(bgrs_per_pixel)]*len(aperture_areas))
    
    # For each aperture sum, get photometry with each background annulus size.
    for ia,(ap_sum,area) in enumerate(zip(aperture_sums,aperture_areas)):
        for ib,bgr_per_pixel in enumerate(bgrs_per_pixel):
            # Calculate photometry
            bgr = bgr_per_pixel*area
            phot[ia,ib]=ap_sum-bgr
            # Explain what elements of output are
            explainer = str(str(radii[ia].value)+
                        '_'+str(annuli_inner_radii[ib].value)+
                        '_'+str(annuli_outer_radii[ib].value)  
                        )
            phot_explainer[ia,ib]=explainer

    return phot, phot_explainer

def get_photometry(f,RAs,DECs,
                   radii,annulus_inners,annulus_outers):
    (aperture_sums, 
     aperture_areas, 
     bgrs_per_pixel) = get_photometry_ingredients(f,RAs,DECs,
                                                  radii,annulus_inners,annulus_outers)
    phot,key = do_photometry_cals(aperture_sums, aperture_areas, bgrs_per_pixel,
                                 radii,annulus_inners,annulus_outers)
    return phot, key

##############################################################
####################### Main Function ########################
##############################################################

def make_lc(RADECfile, fitsfiles, savedir=savedir, verbose=False,debugmode=False,quietmode=False):

    # Create output directory if it doesn't exist.
    if not os.path.isdir(savedir):
        os.makedirs(savedir)

    # Read in RADECfile
    cand_list = ascii.read(candlist_file)
    RAs = regions['RA_sci']
    DECs = regions['DEC_sci']
    

    lc_raw = np.zeros([len(fitsfiles),len(RAs)])

    for ii,f in enumerate(fitsfiles):

        # Determine saveloc for subtracted image based on input savedir and input file name
        fname        = ntpath.basename(fitspath)
        fname_stub   = Path(fname).stem
        saveloc      = savedir + os.path.sep + fname_stub + '_sub.fits'
        textout_path = 
        pngout_path  = 

        # get dt 
        dt = datetime.strptime(sunset_CT.iso, '%Y-%m-%d %H:%M:%S.%f') #2020-10-23 16:01:22.058


        # get mjd

        # Get fits file data
        d,h = fits.getdata(f,header=True)
        w   = WCS(h)

        # get photometry for all RA and DECs
        radii = [5.*u.arcsec]
        annulus_inners = [11. * u.arcsec]
        annulus_outers = [16. * u.arcsec]
        phot,_ = get_photometry(f,RAs,DECs,radii,annulus_inners,annulus_outers)
        lc_raw[ii] = phot[0,0,:]

    # save lc text file


    # save lc png file
    

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
    savedir         = arguments['--out']
    RADECfile       = arguments['<RADECfile>']
    fitsfiles       = arguments['<fitsfiles>']

    _ = make_lc(RADECfile, fitsfiles, savedir=savedir, verbose=verbose,debugmode=debugmode,quietmode=quietmode)
