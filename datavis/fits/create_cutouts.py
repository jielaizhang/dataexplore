#!/usr/bin/env python

"""create_cutouts.py - Create cutouts based on max, min RA and max, min dec or Centre RA/DEC and radius. Cutout fits file will have the right wcs information. 

Usage:
    create_cutouts [-h] [-v] [-d] [--corners STRING] [-centre STRING] [-radius STRING] [-o DIRECTORY] <image>...

Options:
    -h, --help                  Show this screen [default: False]
    -v, --verbose               Show extra information [default: False]
    -d, --debug                 Show debugging information [default: False]
    --corners STRING            Cutout's min RA, max RA, min DEC, max DEC e.g. 20.0,31.0,-28,-24
    --centre STRING             Cutout's central RA, DEC e.g. 34.3,-30.0
    --radius STRING             Cutout's radius if -C. e.g. 3,arcmin or 30,arcsec or 1,deg
    -o DIRECTORY, --outdir DIRECTORY    Output directory where cutout images are saved [default: .]
"""

import docopt
import astropy.io.fits as fits
from astropy import wcs
import numpy
import os

def sexi2deci(RA,dec):
    # Not used
    # Convert sexigesimal string coords to decimal floats
    #RAmax_deci,decmax_deci  = sexi2deci(RAmax,decmax)
    #RAmin_deci,decmin_deci  = sexi2deci(RAmin,decmin)

    # Calculate RA in decimal format
    RA_hr   = float(RA.split(':')[0])
    RA_min  = float(RA.split(':')[1])
    RA_sec  = float(RA.split(':')[2])
    RA_deci = (RA_hr + RA_min/60.0 + RA_sec/60.0/60.0)*15 

    # Calculate dec in decimal format
    dec_deg = float(dec.split(':')[0])
    dec_min = float(dec.split(':')[1])
    dec_sec = float(dec.split(':')[2])

    # Take into account if dec is neg or pos
    dec_sign = numpy.sign(dec_deg)
    if dec_sign > 0:
        dec_deci = dec_deg + dec_min/60.0 + dec_sec/60.0/60.0
    else:
        dec_deci = dec_deg - dec_min/60.0 - dec_sec/60.0/60.0

    return RA_deci, dec_deci

def mkdirp(dirpath,,debug=False):
    if not os.path.isdir(dirpath):
        os.makedirs(dirpath)
        if debug:
            print(f'DEBUG: created directory {dirpath}')
    return None

def create_cutouts_centre(fitsfile,RADECcentre,radius,outdir,verbose=False,debug=False):
    return cutout, cutout_file

def create_cutouts_corners(filesfile,RADECcorners,outdir,verbose=False,debug=False):

    # Read in RA DEC min maxes
    RAmin,RAmax,DECmin,DECmax = [float(x) for x in RADECcorners.split(',')]

    # read in image and WCS
    d,h = fits.getdata(image,header=True)
    w = wcs.WCS(h)

    # Convert decimial RA, dec to pixel values
    [[x1,y1],[x2,y2]]=w.wcs_world2pix([ [RAmax_deci,decmax_deci],
                                        [RAmin_deci,decmin_deci] ],1)

    # Order x, y pixel values to max and min
    xmax = int(max(x1,x2))
    xmin = int(min(x1,x2))
    ymax = int(max(y1,y2))
    ymin = int(min(y1,y2))

    # print out pixel cutout values if verbose
    if debug:
        print('For image ' + image + ' the cutout limits in pixel values are:')
        print('xmax: ',xmax)
        print('xmin: ',xmin)
        print('ymax: ',ymax)
        print('ymin: ',ymin)
    
    # Grab cutout of larger image, WCS does not need to update. 
    # Rotation, scaling should remain the same with the reference pixel changed.
    new_d = d[ymin:ymax,xmin:xmax]
    new_h = h
    # Don't use this because it changes CD to PC matrices
    # new wcs = w[ymin:ymax,xmin:xmax].to_header()
    new_h['CRPIX2']=h['CRPIX2']-ymin
    new_h['CRPIX1']=h['CRPIX1']-xmin

    # Save cutout image with correct wcs
    basename    = image.split('/')[-1]
    stemname    = basename.split('.')[0]
    saveloc     = outdir+'/'+stemname+'_cutout.fits'
    fits.writeto(saveloc, new_d, new_h, overwrite=True)    
    
    return new_d, new_h, saveloc


def create_cutouts(in_images,RADECcorners,RADECcentre,radius,outdir,verbose=False,debug=False):
    mkdirp(outdir)

    for fitsfile in in_images:
        

    
    return cutout, cutout_file

#################### Main part of the program ####################

if __name__=='__main__':

    # Read in arguments
    arguments   = docopt.docopt(__doc__)
    if debug:
        print(arguments)
    in_images   = arguments['<image>']
    outdir      = arguments['--outdir']
    debug       = arguments['--debug']

    # Make cutouts for each image and write out correct wcs
    for image in in_images:


