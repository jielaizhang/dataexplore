#!/usr/bin/env python

""" -- 
Usage: extract_extensions.py [-h] [-v] <multiextfits> <extensions>
        
Options:
    -h, --help                              Show this screen
    -v, --verbose                           Show extra information [default: False]   

Examples:
python extract_extensions.py f.fits 1,2,3,4
"""

import astropy.io.fits as fits
import docopt

def extract_extension(input_path, output_path, n):

    with fits.open(input_path, mode='readonly') as hdu:
        primary_header   = hdu[0].header
        extension_data   = hdu[n].data
        extension_header = hdu[n].header
        extension_header += primary_header

    fits.writeto(output_path, extension_data, extension_header,
                            output_verify='fix')

def extract_extensions(multiextfits,extensions,verbose=False):
    list_exts = [int(x) for x in extensions.split(',')]
    for n in list_exts:
        out = multiextfits.replace('.fits',f'_ext{n}.fits')
        extract_extension(multiextfits, out, n)
        if verbose:
            print(f'Saved: {out}')

if __name__ == "__main__":

    # Read in input arguments
    arguments       = docopt.docopt(__doc__)
    verbose         = arguments['--verbose']
    multiextfits    = arguments['<multiextfits>']
    extensions      = arguments['<extensions>']

    _ = extract_extensions(multiextfits,extensions,verbose=verbose)

