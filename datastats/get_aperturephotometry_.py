#from datastats.get_aperturephotometry_ import get_aperturephotometry
#phot_use,key_use                  = get_aperturephotometry(f,unique_cand_RAs_1,unique_cand_DECs_1,radii,
#                                                           annulus_inners,annulus_outers)

from astropy import units as u
from astropy.coordinates import SkyCoord
from photutils import SkyCircularAperture, SkyCircularAnnulus, CircularAperture
import numpy as np
import astropy.io.fits as fits
from photutils import aperture_photometry
from astropy.wcs import WCS
from astropy.stats import sigma_clipped_stats

# Given RA, DEC, 
# get x_centre,y_centre,aperture_sum as an array, aperture area as an array, background_median as an array
def get_photometry_ingredients(fits_file,RAs,DECs,aperture_radii,annuli_inner_radii,annuli_outer_radii):
    
    # Get data, wcs
    data,h = fits.getdata(fits_file,header=True)
    w      = WCS(h)    
    
    # Define positions where aperture photometry will be performed
    positions = SkyCoord(ra=RAs, dec=DECs,unit='deg')
    
    # Create apertures at positions for photometry
    aperture_sky = [SkyCircularAperture(positions, r=r) for r in aperture_radii]
    
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

def get_aperturephotometry(f,RAs,DECs,
                            radii,annulus_inners,annulus_outers):
    (aperture_sums, 
     aperture_areas, 
     bgrs_per_pixel) = get_photometry_ingredients(f,RAs,DECs,
                                                  radii,annulus_inners,annulus_outers)
    phot,key = do_photometry_cals(aperture_sums, aperture_areas, bgrs_per_pixel,
                                 radii,annulus_inners,annulus_outers)
    return phot, key
