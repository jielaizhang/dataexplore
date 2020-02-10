#!/usr/bin/env python

"""create_cutouts.py - Create cutouts based on max, min x-coords and max, min y-coords. Cutout fits file will have the right wcs information. 

Usage:
    create_cutouts [-h] [-v] [-X NUMBER] [-x NUMBER] [-Y NUMBER] [-y NUMBER] [-o DIRECTORY] [-w] <image>

Options:
    -h, --help                          Show this screen [default: False]
    -v, --verbose                       Show extra information [default: False]
    -X NUMBER, --Xmax NUMBER            Cutout's max X
    -x NUMBER, --xmin NUMBER            Cutout's min x
    -Y NUMBER, --Ymax NUMBER            Cutout's max Y
    -y NUMBER, --ymin NUMBER            Cutout's min y
    -o PATH, --out PATH                 Output file (full path) where cutout images are saved [default: ./cutout.fits]
    -w, --nowcscopy                     Don't get correct wcs of cutout saved in output header [default: False]
"""
import docopt
import astropy.io.fits as fits
from astropy import wcs
import numpy
import sys

def create_cutout_bypixel(image,saveloc,xmax,xmin,ymax,ymin,nowcscopy=False,verbose=False):
    """Create cutouts based on max, min x-coords and max, min y-coords. Cutout fits file will have the right wcs information if nowcopy=False."""

    printme = '\nimage={} \nsaveloc={} \nxmax={} \nxmin={} \nymax={} \nymin={}'.format(
                image,saveloc,xmax,xmin,ymax,ymin)
    if verbose:
        print(printme)

    xmax    = int(xmax)
    xmin    = int(xmin)
    ymax    = int(ymax)
    ymin    = int(ymin)   

    # read in images
    f = fits.open(image)
    w = wcs.WCS(f[0].header)

    # Grab cutout of larger image, including cutout wcs
    newf=fits.PrimaryHDU()
    newf.header=f[0].header
    newf.data=f[0].data[ymin:ymax,xmin:xmax]
    if nowcscopy:
        print('wcs not saved in output image')
    else:
        newf.header.update(w[ymin:ymax,xmin:xmax].to_header())

    # Save cutout image with correct wcs
    newf.writeto(saveloc, overwrite=True)
    print('File saved at: {}\n'.format(saveloc))

#################### Main part of the program ####################

if __name__=='__main__':

    # Read in arguments
    arguments   = docopt.docopt(__doc__)

    image   = arguments['<image>']
    saveloc = arguments['--out']
    
    xmax    = arguments['--Xmax']
    xmin    = arguments['--xmin']
    ymax    = arguments['--Ymax']
    ymin    = arguments['--ymin']   

    nowcscopy = arguments['--nowcscopy']
    verbose = arguments['--verbose']

    create_cutout_bypixel(image,saveloc,xmax,xmin,ymax,ymin,nowcscopy=nowcscopy,verbose=verbose)
