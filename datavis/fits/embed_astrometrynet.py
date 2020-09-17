#!/usr/bin/env python

"""embed_astrometrynet - store astrometric solution in all light images in a directory

Usage:
    embed_astrometrynet [-h] [-v] [-o directory] [-p ARCSEC] <fitsfiles>...

Options:
    -h, --help                              Show this screen
    -v, --verbose                           Show extra information [default: False]
    -o DIRECTORY, --outputdir DIRECTORY     Output directory name [default: PlateSolved]
    -p ARCSEC, --pixscale ARCSEC            Pixel scale in arcsec per pixel lower and upper limits 
"""

import os, sys
import subprocess
from glob import glob

import docopt

def mkdirp(dirpath):
    if not os.path.isdir(dirpath):
        os.makedirs(dirpath)
    return None

def plate_solve(fileloc,pixscale,outdir,verbose=False):
    if pixscale:
        Lpixscale        = pixscale.split(',')[0]
        Hpixscale        = pixscale.split(',')[1]

    print('solve-field for: '+fileloc)
    # Run solve-field command
    if pixscale: 
        command = 'solve-field -u arcsecperpix -L %s -H %s %s --overwrite' %(Lpixscale,Hpixscale,fileloc)
    else:
        command = 'solve-field %s --overwrite' %(fileloc)
    if verbose:
        print("Solving for WCS for %s" % fileloc)
    subprocess.call(command,shell=True)

    # Determine directory in which solve-field operated
    directory = fileloc.split('/')[0:-1]
    directory = '/'.join(directory)+'/'
    # Determine file name stub (with .fits)
    filestub  = fileloc.split('/')[-1].split('.')[0]

    # Look for .solved file, if found, move fits file with wcs to output directory.
    solved_file  = fileloc.split('.')[0]+'.solved'
    if os.path.exists(solved_file):
        command = 'mv {d}'.format(d=directory)  +filestub + '.new {d}/'.format(d=outdir) + filestub + '_wcs.fits'
        print(command)
        subprocess.call(command,shell=True)
        command = 'mv {d}'.format(d=directory) + filestub + '-objs.png {d}/'.format(d=outdir) + filestub + '-objs.png'
        print(command)
        subprocess.call(command,shell=True)
        command = 'rm {d}*rdls {d}*xyls {d}*axy {d}*wcs {d}*corr {d}*match {d}*solved {d}*png'.format(d=directory)
        print(command)
        subprocess.call(command,shell=True)
        return True
    else:
        command = 'rm {d}*axy'.format(d=directory)
        print(command)
        subprocess.call(command,shell=True)
        return False


####################### BODY OF PROGRAM STARTS HERE ########################

if __name__ == "__main__":

    arguments = docopt.docopt(__doc__)

    # Mandatory argument
    lightfilelist = arguments['<fitsfiles>']

    # Non-mandatory options without arguments
    verbose = arguments['--verbose']

    # Non-mandatory options with arguments
    output_directory = arguments['--outputdir']
    mkdirp(output_directory)
    pixscale         = arguments['--pixscale']

    if verbose:
        print( "")
        print( "**** PLATE SOLVING USING ASTROMETRY.NET ****")
        print( "")

    # Create the output directory if necessary
    if not os.path.exists(output_directory):
        os.makedirs(output_directory)

    # Create a list of all light frames elegantly exit if no images found in input directory.
    if len(lightfilelist) == 0:
	    print( 'There are no files in this directory. Likely no dark subtracted flats exist.')
	    sys.exit('ERROR: No files found.')

    # Loop over the lights and plate solve each one
    for fname in lightfilelist:
        if verbose:
            print('================')
            print(f'Working on {fname}')
            print('================')

        print(' ')

        # attempt to plate solve
        if plate_solve(fname,pixscale,output_directory,verbose=verbose):
            print( "Astrometric solution succeded for: " +fname)
        else:
            print( "Astrometric solution failed for: " +fname)
