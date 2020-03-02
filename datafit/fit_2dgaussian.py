#!/usr/bin/env python

"""fit_2dgaussian.py -- given  

Usage: fit_2dgaussian [-h] [-v] [-s STRING] [-o] <fitsfile> <init_params>

Arguments:
    fitsfile STRING
        Containing the Gaussian you want to fit
    init_params STRING
        amplitude,x0,y0,sigma_x,sigma_y,theta (from +y axis, clockwise),offset. 

Options:
    -h, --help                                  Show this screen
    -v, --verbose                               Show extra information (including fitted parameters) [default: False]   
    -s STRING, --saveloc STRING                 Save model as fits file, e.g. ./model.fits 
    -o, --overwrite                             Overwrite if saveloc exists.   

Examples:
    python fit_2dgaussian.py gaussian2fit.fits 3,100,100,20,40,0,10 
"""

import docopt
import astropy.io.fits as fits

import numpy as np
from scipy.optimize import curve_fit

def twoD_Gaussian(xytuple, amplitude, xo, yo, sigma_x, sigma_y, theta, offset):
    (x,y) = xytuple
    xo    = float(xo)
    yo    = float(yo)    
    a     = (np.cos(theta)**2)/(2*sigma_x**2) + (np.sin(theta)**2)/(2*sigma_y**2)
    b     = -(np.sin(2*theta))/(4*sigma_x**2) + (np.sin(2*theta))/(4*sigma_y**2)
    c     = (np.sin(theta)**2)/(2*sigma_x**2) + (np.cos(theta)**2)/(2*sigma_y**2)
    g     = offset + amplitude*np.exp( - (a*((x-xo)**2) + 2*b*(x-xo)*(y-yo) 
                            + c*((y-yo)**2)))
    return g.ravel()

def fit_2dgaussian(z,parameters_initial_guess,verbose=False):
    '''Fit 2D gaussian to z (lenx by leny array). 
    The parameters_initial_guess should include:
    amplitude, 
    x0, 
    y0, 
    sigma_x, 
    sigma_y, 
    theta (from +y axis, clockwise), 
    offset.'''
    
    # Create x and y indices
    lenx,leny = np.shape(z)
    x         = np.linspace(0, lenx-1, lenx)
    y         = np.linspace(0, leny-1, leny)
    x, y      = np.meshgrid(x, y)
    xytuple   = tuple((x,y))
    
    # Unwrap z
    z=z.ravel()
    
    # Do fit
    popt, pcov = curve_fit(twoD_Gaussian, xytuple, z, p0=parameters_initial_guess)

    # Print parameters
    if verbose:
        print('The model parameters are: ',popt)
        print('amplitude: ',popt[0])
        print('x0: ',popt[1])
        print('y0: ',popt[2])
        print('sigma_x: ',popt[3])
        print('sigma_y: ',popt[4])
        print('theta (clockwise from y-axis: ',popt[5])
        print('offset: ',popt[6])

    # Calculate model (i.e. create gaussian with fitted params)
    model            = twoD_Gaussian(xytuple, *popt)
    model            = model.reshape(lenx,leny)
    model_parameters = popt
    
    return model, model_parameters, x, y

def fit_2dgaussian_fitsinput(fitsfile,init_params,save=False,overwrite=False,verbose=False):
    # Read in fits file
    z                 = fits.getdata(fitsfile)
    # Unwrap initial parameters and make floats, then wrap it up
    [amplitude,x0,y0,sigma_x,sigma_y,theta,offset] = [float(x) for x in init_params.split(',')]
    initial_guess     = (amplitude,x0,y0,sigma_x,sigma_y,theta,offset)
    # Do fit
    data_fitted, popt = fit_2dgaussian(z,initial_guess,verbose=verbose)
    hduP=fits.PrimaryHDU(im)
    hdulist=fits.HDUList(hduP)
    keynames = ['amplitude','x0','y0','sigma_x','sigma_y','theta','offset']
    for name, value in zip(keynames,popt):
        hdulist[0].header.update(name,value)
    if overwrite:
        hdulist.writeto(saveloc,overwrite=True)
        print('File saved (overwrite=True) ',saveloc)
    else:
        hdulist.writeto(saveloc)
        print('File saved (overwrite=False) ',saveloc)
    return

####################### BODY OF PROGRAM STARTS HERE ########################

if __name__ == "__main__":

    arguments = docopt.docopt(__doc__)

    # Mandatory argument
    fitsfile    = arguments['<fitsfile>']
    init_params = argumetns['<init_params>']

    # Non-mandatory options without arguments
    verbose     = arguments['--verbose']
    save        = arguments['--saveloc']
    overwrite   = arguments['--overwrite']
    
    # Perform fit
    fit_2dgaussian_fitsinput(fitsfile,init_params,save=save,overwrite=overwrite,verbose=verbose)

