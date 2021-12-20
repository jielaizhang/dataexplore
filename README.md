# dataexplore
Scripts I keep using over and over again in different projects and disciplines for visualizing and exploring data. 

This repository is written in the vision of providing a set of standalone scripts that can be used in the terminal to do specific data visualisation tasks that are done over and over again. They are written to allow the user to do simple things without having to leave the bash environment. The scripts are written to be importable into a python environment if the user wants to. 

There are three modules. `datastats`, `datavis`, `datafit`. `datastats` contain functions that calculate statistical properties of input images. `datavis` contains functions that help visualise input data, and is split into folders depending on the type of input data format: ascii, fits, txt. `datafit` contains functions that help fit the input data. 

Example import command:
`from datavis.fits.convolve import convolve as dvfconv`

Example of how to get help on how to use scripts in command line:
`python convolve.py -h`

Anaconda install of python 3.6

Required modules:
docopt
astropy
matplotlib
numpy
os
sys
scipy
pylab

# Notes on Image Subtraction
Given two images with overlapping field of view, and comparible (but not identical) quality, you can perform image subtraction in the following way.

1. Install- downlaod the dataexplore repository, and add it to your PYTHONPATH (put `export PYTHONPATH=$PYTHONPATH:/path/to/dataexplore` in your .bashrc, or type into the command line upon launch). 
2. Align two images and crop out maximum overlapping area using dataexplore/datavis/fits/align_image.py (`python align_image.py -h` for instructions)
3. Perform image subtraction of two aligned images using dataexplore/datavis/subtract_image.py (`python subtract_image.py` -h for instructions)

Required software for above to work: SWarp (https://www.astromatic.net/software/swarp/), Source Extractor (https://www.astromatic.net/software/sextractor/), High Order Transform of Psf ANd Template Subtraction (https://github.com/acbecker/hotpants)


