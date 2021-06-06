#!/usr/bin/python

#!/usr/bin/env python

""" convert_wavelength_frequency_energy.py -- Input a list of energies in eV, a list of wavelengths in meter or a list of frequencies in Hz. Print out a table of eV, meter and Hz. If used in python environment, directly access the conversion functions. 

Usage: convert_wavelength_frequency_energy.py [-h] [-v] [--debug] (eV | meter | Hz) [-s INT] <values>...

Arguments:
    values (float)
        input a list of energies in eV, a list of wavelengths in meter or a list of frequencies in Hz.

Options:
    -h, --help                              Show this screen
    -v, --verbose                           Show extra information [default: False]   
    -q, --quietmode                         Minimize print to screen. This is useful when this function is called in another function. [default: False]  
    --debug                                 Output more for debugging [default: False]
    -s INT, --sigfigs INT                   Number of sig figs to print [default: 3]

Examples: 
Bash:
Python:
    from misc.convert_wavelength_frequency_energy import m2eVandHz, eV2mandHz, Hz2eVandm
    from misc.convert_wavelength_frequency_energy import eVs_to_Hzs
"""
import docopt
import scipy 
import scipy.constants


__author__      = "Jielai Zhang"
__license__     = "MIT"
__version__     = "1.0.1"
__date__        = "2021-01-25"
__maintainer__  = "Jielai Zhang"
__email__       = "zhang.jielai@gmail.com"

##############################################################
####################### Main Functions #######################
##############################################################

def eVs_to_Hzs(eVs):
    nu_list = []
    for eV in eVs:
        _,nu = eV2mandHz(eV)
        nu_list.append(nu)
    return nu_list

def m2eVandHz(m,printit=False):
    
    h = scipy.constants.Planck
    c = scipy.constants.speed_of_light
    eV_per_J = scipy.constants.physical_constants['joule-electron volt relationship'][0]
    
    eV = h*c/m*eV_per_J
    Hz = c/m
    
    if printit:
        print('eV :',eV)
        print('Hz :',Hz)
    
    return eV, Hz

def eV2mandHz(eV,printit=False):
    
    h = scipy.constants.Planck
    c = scipy.constants.speed_of_light
    eV_per_J = scipy.constants.physical_constants['joule-electron volt relationship'][0]
    
    m  = h*c/eV*eV_per_J
    Hz = c/m
    
    if printit:
        print('m :',m)
        print('Hz :',Hz)
    
    return m, Hz

def Hz2eVandm(Hz,printit=False):
    
    h = scipy.constants.Planck
    c = scipy.constants.speed_of_light
    eV_per_J = scipy.constants.physical_constants['joule-electron volt relationship'][0]
    
    m  = c/Hz
    eV = h*Hz*eV_per_J
    
    if printit:
        print('eV :',eV)
        print('m :',m)
    
    return eV, m

############################################################################
####################### BODY OF PROGRAM STARTS HERE ########################
############################################################################

if __name__ == "__main__":

    # Read in input arguments
    arguments       = docopt.docopt(__doc__)
    verbose         = arguments['--verbose']
    debugmode       = arguments['--debug']
    input_eV        = arguments['eV']
    input_m         = arguments['meter']
    input_Hz        = arguments['Hz']
    no_sig_figs     = int(arguments['--sigfigs'])
    
    if debugmode:
        print(arguments)   

    if input_eV:
        eVs = arguments['<values>']
        print('eV, m, Hz')
        for eV in eVs:
            m,Hz = eV2mandHz(float(eV))
            print(f"{eV}, {m:.4E}, {Hz:.4E}")

    if input_m:
        meters = arguments['<values>']

    if input_Hz:
        Hzs = arguments['<values>']

