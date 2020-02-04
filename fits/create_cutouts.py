#!/usr/bin/env python

"""create_cutouts.py - Create cutouts based on max, min RA and max, min dec. Cutout fits file will have the right wcs information. 

Usage:
    create_cutouts [-h] [-v] [-R NUMBER] [-r NUMBER] [-D NUMBER] [-d NUMBER] [-o DIRECTORY] <image>...

Options:
    -h, --help                          Show this screen [default: False]
    -v, --verbose                       Show extra information [default: False]
    -R NUMBER, --RAmax NUMBER           Cutout's max RA in sexigesimal notation e.g. 10:20:0
    -r NUMBER, --ramin NUMBER           Cutout's min RA in sexigesimal notation e.g. 9:1:40.3
    -D NUMBER, --Decmax NUMBER          Cutout's max dec in sexigesimal notation e.g. 2:1:50.5
    -d NUMBER, --decmin NUMBER          Cutout's min dec in sexigesimal notation e.g. -1:18:1
    -o DIRECTORY, --outdir DIRECTORY    Output directory where cutout images are saved [default: .]
"""

import docopt
import astropy.io.fits as fits
from astropy import wcs
import numpy
import os

def sexi2deci(RA,dec):

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

def mkdirp(dirpath):
    if not os.path.isdir(dirpath):
        os.makedirs(dirpath)
    return None

#################### Main part of the program ####################

if __name__=='__main__':

    # Read in arguments
    arguments   = docopt.docopt(__doc__)

    in_images   = arguments['<image>']
    out_dir     = arguments['--outdir']
    mkdirp(out_dir)

    RAmax       = arguments['--RAmax']
    RAmin       = arguments['--ramin']
    decmax      = arguments['--Decmax']
    decmin      = arguments['--decmin']

    verbose     = arguments['--verbose']
    if verbose:
        print(arguments)

    # Convert sexigesimal string coords to decimal floats
    RAmax_deci,decmax_deci  = sexi2deci(RAmax,decmax)
    RAmin_deci,decmin_deci  = sexi2deci(RAmin,decmin)

    # Make cutouts for each image and write out correct wcs
    for image in in_images:

        # read in images
        d,h = fits.getdata(image,header=True)
        w = wcs.WCS(h)

        # Convert decimial RA, dec to pixel values
        [[x1,y1],[x2,y2]]=w.wcs_world2pix(
                                    [ [RAmax_deci,decmax_deci],
                                      [RAmin_deci,decmin_deci] ]
                                                    ,1)

        # Order x, y pixel values to max and min
        xmax = int(max(x1,x2))
        xmin = int(min(x1,x2))
        ymax = int(max(y1,y2))
        ymin = int(min(y1,y2))

        # print out pixel cutout values if verbose
        if verbose:
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
        saveloc     = out_dir+'/'+stemname+'_cutout.fits'
        fits.writeto(saveloc, new_d, new_h, overwrite=True)
