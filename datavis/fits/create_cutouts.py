#!/usr/bin/env python

"""create_cutouts.py - Create cutouts based on max, min RA and max, min dec or Centre RA/DEC and radius. Cutout fits file will have the right wcs information. 

Usage:
    create_cutouts [-h] [-v] [-d] [--corners STRING] [--centre STRING] [--radius STRING] [-o DIRECTORY] <image>...

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
import os
import ntpath
from pathlib import Path
from astropy.nddata import Cutout2D
from astropy import units as u
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord

def mkdirp(dirpath,debug=False):
    if not os.path.isdir(dirpath):
        os.makedirs(dirpath)
        if debug:
            print(f'DEBUG: created directory {dirpath}')
    return None

def create_cutout_centre(fitsfile,RA,DEC,image_size,verbose=False,debug=False):
    # Read data
    d,h = fits.getdata(fitsfile,header=True)
    w = WCS(h)
    # Get position in SkyCoords
    pos = SkyCoord(RA, DEC, unit=u.deg)
    # Create cutout
    cutout = Cutout2D(d, size=image_size, position=pos, wcs=w)
    # Create new cutout header with right WCS
    h.update(cutout.wcs.to_header())
    return cutout.data, h

def create_cutout_corners(filesfile,RAmin,RAmax,DECmin,DECmax,verbose=False,debug=False):

    # read in image and WCS
    d,h = fits.getdata(filesfile,header=True)
    w = WCS(h)

    # Convert decimial RA, dec to pixel values
    [[x1,y1],[x2,y2]]=w.wcs_world2pix([ [RAmax,DECmax],
                                        [RAmin,DECmin] ],1)

    # Order x, y pixel values to max and min
    xmax = int(max(x1,x2))
    xmin = int(min(x1,x2))
    ymax = int(max(y1,y2))
    ymin = int(min(y1,y2))

    # print out pixel cutout values if debug
    if debug:
        print('DEBUG: For image ' + image + ' the cutout limits in pixel values are:')
        print('DEBUG: xmax: ',xmax)
        print('DEBUG: xmin: ',xmin)
        print('DEBUG: ymax: ',ymax)
        print('DEBUG: ymin: ',ymin)
    
    # Grab cutout of larger image, WCS does not need to update. 
    # Rotation, scaling should remain the same with the reference pixel changed.
    new_d = d[ymin:ymax,xmin:xmax]
    new_h = h
    # Don't use this because it changes CD to PC matrices
    # new wcs = w[ymin:ymax,xmin:xmax].to_header()
    new_h['CRPIX2']=h['CRPIX2']-ymin
    new_h['CRPIX1']=h['CRPIX1']-xmin  
    
    return new_d, new_h


def create_cutouts(in_images,RADECcorners,RADECcentre,radius,outdir,verbose=False,debug=False):
    '''this function is mainly for use when doing stuff in terminal'''
    mkdirp(outdir)

    if debug:
        print('DEBUG: debugging...')

    cutout_files = []

    if RADECcentre != None:
        if debug:
            print(f'DEBUG: RADECcentre = {RADECcentre}')
        # Get right image size
        radius,unit = radius.split(',')
        radius = float(radius)
        if unit == 'arcsec':
            image_size = radius*u.arcsec
        elif unit == 'arcmin':
            image_size = radius*u.arcmin
        elif unit == 'deg':
            image_size = radius*u.deg
        else:
            print('Warning: assuming radius is entered in ARCSEC units.')
            image_size = radius*u.arcsec    
        # Get image centre
        RA,DEC = [float(x) for x in RADECcentre.split(',')]

    elif RADECcorners != None:
        if debug:
            print(f'DEBUG: RADECcorners = {RADECcorners}')
        # Read in RA DEC min maxes
        RAmin,RAmax,DECmin,DECmax = [float(x) for x in RADECcorners.split(',')]


    for fitsfile in in_images:
        # Create cutout
        if RADECcentre != None:
            cutout,cutout_header = create_cutout_centre(fitsfile,
                                                        RA,DEC,
                                                        image_size,
                                                        outdir)
        elif RADECcorners != None:
            cutout,cutout_header = create_cutout_corners(fitsfile,
                                                            RAmin,RAmax,DECmin,DECmax,
                                                            verbose=verbose,debug=debug)
        
        # Save cutout image with correct wcs
        stub        = Path(ntpath.basename(fitsfile)).stem
        saveloc     = outdir + os.path.sep + stub + '_cutout.fits'
        fits.writeto(saveloc, cutout, cutout_header, overwrite=True) 
        if verbose:
            print(f'VERBOSE: saved {saveloc}') 
        cutout_files.append(saveloc)

    return cutout_files

#################### Main part of the program ####################

if __name__=='__main__':

    # Read in arguments
    arguments   = docopt.docopt(__doc__)
    debug       = arguments['--debug']
    if debug:
        print(arguments)
    verbose      = arguments['--verbose']
    in_images    = arguments['<image>']
    outdir       = arguments['--outdir']
    RADECcorners = arguments['--corners']
    RADECcentre  = arguments['--centre']
    radius       = arguments['--radius']

    _ = create_cutouts(in_images,RADECcorners,RADECcentre,radius,outdir,verbose=verbose,debug=debug)
