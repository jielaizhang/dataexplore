#!/usr/bin/env python

"""plot_histogram.py -- Plot a single column of an input text file as a histogram. 

Usage: predict_class [-h] [-v] [--debug] [-o FILE] [--delimiter STRING] [-b INTEGER] <textfile> <column>

Arguments:
    textfile (string)
        Text file containing a column you want to make a histogram of.
    column (integer)
        Python 0-index column number (the one you want to make a histogram of). 

Options:
    -h, --help                              Show this screen
    -v, --verbose                           Show extra information [default: False]     
    --debug                                 Output more for debugging [default: False]
    -o FILE, --out FILE               		Save output as [default: ./histogram.png]
    --delimiter string                      Text file delimiter [default: None]
    -b NUMBER, --bins number                The number of bins in histogram [default: None]

Examples:
from datavis.txt.plot_histogram import plot_histogram_array
"""

import docopt
import os
import numpy as np
import matplotlib.pyplot as plt
import pylab
from matplotlib.ticker import FormatStrFormatter

def plot_histogram_array(data_toplot, outfile=False, n_bins=None, delimiter=None, verbose=False):
    '''output png file containing 
    a histogram of given column of text file'''


    # Determine n_bins
    if n_bins == None:
        n_bins = int(len(data_toplot)/10.)
    else:
        n_bins = int(n_bins)

    # Plot histogram
    print('Creating histogram, please wait.')
    fig = plt.figure(1)
    ax  = fig.add_subplot(111)
    n,bins,patches = pylab.hist(data_toplot,n_bins)
    ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    print('Histograph created: ',outfile)
    
    if outfile:
        plt.savefig(outfile)
    else:
        plt.show()

    # Print out info
    if verbose:
        print('VERBOSE: File saved: ',outfile)

    return None

def plot_histogram(textfile, column, outfile='histogram.png', n_bins=None, delimiter=None, verbose=False):
    '''Input text file and column number, output png file containing 
    a histogram of given column of text file'''


    # Read in text file, and pick out th relevant column. 
    data        = np.transpose(np.genfromtxt(textfile,delimiter=delimiter,dtype='str'))
    data_toplot = [float(x) for x in data[int(column)]]

    # Determine n_bins
    if n_bins == None:
        n_bins = int(len(data_toplot)/10.)
    else:
        n_bins = int(n_bins)

    # Plot histogram
    print('Creating histogram, please wait.')
    fig = plt.figure(1)
    ax  = fig.add_subplot(111)
    n,bins,patches = pylab.hist(data_toplot,n_bins)
    ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    print('Histograph created.')
    plt.savefig(outfile)

    # Print out info
    if verbose:
        print('VERBOSE: File saved: ',outfile)

    return None

####################### BODY OF PROGRAM STARTS HERE ########################

if __name__ == "__main__":

    arguments = docopt.docopt(__doc__)

    # Mandatory argument
    textfile  = arguments['<textfile>']
    column    = int(arguments['<column>'])

    # Non-mandatory options
    verbose     = arguments['--verbose']
    debugmode   = arguments['--debug']
    outfile     = arguments['--out']
    delimiter   = arguments['--delimiter']
    n_bins      = arguments['--bins']

    if debugmode:
        print(arguments)

    # Convert None strings to NoneType
    delimiter = None if delimiter == 'None' else delimiter
    n_bins = None if n_bins == 'None' else n_bins

    plot_histogram(textfile, column, outfile=outfile, n_bins=n_bins, delimiter=delimiter, verbose=verbose)
