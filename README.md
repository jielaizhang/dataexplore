# dataexplore
Scripts I keep using over and over again in different projects and disciplines for visualizing and exploring data. 

This repository is written in the vision of providing a set of standalone scripts that can be used in the terminal to do specific data visualisation tasks that are done over and over again. They are written to allow the user to do simple things without having to leave the bash environment. The scripts are written to be importable into a python environment if the user wants to. 

The different directories contain scripts that work on different input file formats. 

Example import command:
from datavis.fits.convolve import convolve as dvfconv

Example of how to get help on how to use scripts in command line:
python convolve.py -h

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
