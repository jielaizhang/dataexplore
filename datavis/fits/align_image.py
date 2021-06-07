#!/usr/bin/python

#!/usr/bin/env python

""" align_image.py -- Input two fits images, determine their maximum overlapping area, swarp resample to align.  

Usage: align_image.py [-h] [-v] [--debug] [-q] [--swarp LOC] [-o SAVELOC] <fitsfile1> <fitsfile2>

Arguments:
    fitsfile1 (string)
    fitsfile2 (string)

Options:
    -h, --help                              Show this screen
    -v, --verbose                           Show extra information [default: False]   
    -q, --quietmode                         Minimize print to screen. This is useful when this function is called in another function. [default: False]  
    --debug                                 Output more for debugging [default: False]
    -o SAVEDIR, --out SAVEDIR               Saved output as [default: ./]
    --swarp LOC                             Location of swarp [default: /opt/local/bin/swarp]

Examples:
    Bash: python align_image.py f1.fits f2.fits -o ./
    Python: from datavis.fits.align_image import align_image
            f_out1,f_out2 = align_image(f1,f2,savedir='./',swarploc='/opt/local/bin',verbose=False,debugmode=False,quietmode=False)
    
"""
import docopt
import sys, os, ntpath
import subprocess

# Jielai's modules
from datavis.fits.determine_imageOverlap import determine_imageOverlap


__author__      = "Jielai Zhang"
__license__     = "MIT"
__version__     = "1.0.1"
__date__        = "2020-06-26"
__maintainer__  = "Jielai Zhang"
__email__       = "zhang.jielai@gmail.com"

#########################################
# ======= House Keeping Functions =======
#########################################
'''These functions standardise verbose, debug printing'''
def print_verbose_string(printme,verbose=False,underscores=False):
    if verbose:
        if underscores:
            print("-" * len(f"VERBOSE: {printme}"),file=sys.stdout)
            print(f"VERBOSE: {printme}",file=sys.stdout)
            print("-" * len(f"VERBOSE: {printme}"),file=sys.stdout)
        else:
            print(f"VERBOSE: {printme}",file=sys.stdout)

def print_debug_string(printme,debugmode=False,underscores=False):
    if debugmode:
        if underscores:
            print("-" * len(f"DEBUG  : {printme}"),file=sys.stdout)
            print(f"DEBUG  : {printme}",file=sys.stdout)
            print("-" * len(f"DEBUG  : {printme}"),file=sys.stdout)
        else:
            print(f"DEBUG  : {printme}",file=sys.stdout)

'''These functions help organise or delete files'''
def clearit(fnames,debugmode=False):
    for fname in fnames:
        if os.path.isfile(fname):
            os.remove(fname)
            printme = f'Temporary file deleted: {fname}.'
            print_debug_string(printme,debugmode=debugmode)
    return None

##############################################################
####################### Main Function ########################
##############################################################

def align_image(f1,f2,savedir='./',swarploc='/opt/local/bin',verbose=False,debugmode=False,quietmode=False):

    # Write down what will be saved, even if temporarily 
    fname1          = ntpath.basename(f1)
    fname2          = ntpath.basename(f2)
    f_out1          = savedir+os.sep+fname1.replace('.fits','.resamp.fits')
    f_out2          = savedir+os.sep+fname2.replace('.fits','.resamp.fits')
    f_outweight1    = savedir+os.sep+fname1.replace('.fits','.resamp.weight.fits')
    f_outweight2    = savedir+os.sep+fname2.replace('.fits','.resamp.weight.fits')
    f_outxml        = '.'+os.sep+'swarp.xml'
    swarpdefault = './default.swarp'

    # Determine maximum overlap amount of two images in pixel coordinates
    XLEN,YLEN = determine_imageOverlap( f1,f2,plotsave=False,
                                        verbose=verbose,debugmode=debugmode,quietmode=quietmode)

    # SWarp align verbosity set
    if quietmode == True:
        verbosity = 'QUIET'
    elif debugmode == True:
        verbosity = 'FULL'
    else:
        verbosity = 'NORMAL'

    # Create default swarp config file, so that no error is printed
    command      = f'{swarploc} -d > {swarpdefault}'
    subprocess.call(command,shell=True)
    printme = f'Temporarily saved: {swarpdefault}.'
    print_debug_string(printme,debugmode=debugmode)

    # Run SWarp
    # swarpcommand = f'{swarploc} {f1} {f2} -SUBTRACT_BACK N -RESAMPLE Y -COMBINE N -CENTER_TYPE MOST -IMAGE_SIZE {XLEN},{YLEN} -RESAMPLE_DIR {savedir} -VERBOSE_TYPE {verbosity}'
    swarpcommand = f'{swarploc} {f1} {f2} -VMEM_MAX 10000 -MEM_MAX 5000 -COMBINE_BUFSIZE 5000 -SUBTRACT_BACK N -RESAMPLE Y -COMBINE N -CENTER_TYPE MOST -IMAGE_SIZE {XLEN},{YLEN} -RESAMPLE_DIR {savedir} -VERBOSE_TYPE {verbosity}'
    subprocess.call(swarpcommand,shell=True)

    # If debug, print temporary files that were deleted
    printme = f'Temporarily saved: {f_outweight1}.'
    print_debug_string(printme,debugmode=debugmode)
    printme = f'Temporarily saved: {f_outweight2}.'
    print_debug_string(printme,debugmode=debugmode)
    printme = f'Temporarily saved: {f_outxml}.'
    print_debug_string(printme,debugmode=debugmode)
    
    # If debug, print swarp command
    printme = f'SWarp command used:'
    print_debug_string(printme,debugmode=debugmode)    
    print_debug_string(swarpcommand,debugmode=debugmode)    

    # Clean weights and xml files and swarp default
    clearit([f_outweight1,f_outweight2,f_outxml,swarpdefault],debugmode=debugmode)

    # Print what was saved
    if not quietmode:
        print(f'SAVED  : {f_out1}')
        print(f'SAVED  : {f_out2}')
  
    return f_out1,f_out2

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
    fitsfile1       = arguments['<fitsfile1>']
    fitsfile2       = arguments['<fitsfile2>']
    swarploc        = arguments['--swarp']

    _ = align_image(fitsfile1,fitsfile2,savedir=savedir,swarploc=swarploc,verbose=verbose,debugmode=debugmode,quietmode=quietmode)
