#!/usr/bin/env python

"""plot_crosssection.py -- read in fits file and 1D plot the values of pixels along a specified cross section. 

Usage: plot_crosssection [-h] [-v] [-o STRING] [--direction STRING] [--hline ARRAY] [--vline ARRAY] [--gridoff] <fitsfile> <slicenumber>

Options:
    -h, --help                          Show this screen
    -v, --verbose                       Show extra information [default: False]      
    -o STRING, --outfile STRING         Where plot is saved [default: ./crossection.png]
    --direction STRING                  horizontal or vertical [default: horizontal]
    --hline ARRAY                       E.g. 4
    --vline ARRAY                       E.g. 10,15
    --gridoff                           Major and minor grids off [default: False]

Examples:
    python plot_crosssection.py ref_psf.fits 12 --hline [0.25] --vline [9,14]  
"""

import docopt, sys
import matplotlib.pyplot as plt
import astropy.io.fits as fits
import numpy as np

def plot_crosssection(fits_name,slicenumber,
                      direction='horizontal',saveloc='./crosssection.png',
                      horizline=False, vertline=False,
                      gridoff=False,
                      title=False):
    d = fits.getdata(fits_name)
    if direction == 'horizontal':
        s = d[slicenumber,:]
    elif direction == 'vertical':
        s = d[:,slicenumber]
    else:
        sys.exit('Please specify horizontal or vertical for the direction.')
    l = len(s)
    plt.plot(np.linspace(0,l-1,l),s)
    if not gridoff:
        plt.minorticks_on()
        plt.grid(True,which='Major')
        plt.grid(True,which='Minor',linestyle='-', alpha=0.35)
    if horizline:
        for i in horizline:
            plt.axhline(i,alpha=0.4)
    if vertline:
        for i in vertline:
            plt.axvline(i,alpha=0.4)
    if title:
        plt.title(title)
    if saveloc:
        plt.savefig(saveloc)
    else:
        plt.show()
    return None

####################### BODY OF PROGRAM STARTS HERE ########################

if __name__ == "__main__":

    arguments = docopt.docopt(__doc__)

    # Mandatory argument
    fits_name     = arguments['<fitsfile>']
    slicenumber   = int(arguments['<slicenumber>'])

    # Non-mandatory options without arguments
    verbose     = arguments['--verbose']
    outfile     = arguments['--outfile']
    direction   = arguments['--direction']
    hlinearray  = [float(x) for x in arguments['--hline'].split(',')]
    vlinearray  = [float(x) for x in arguments['--vline'].split(',')]
    gridoff     = arguments['--gridoff']

    if verbose:
        print(arguments)

    plot_crosssection(fits_name,slicenumber,
                      direction=direction,saveloc=outfile,
                      horizline=hlinearray, vertline=vlinearray,
                      gridoff=gridoff)
