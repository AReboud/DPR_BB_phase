# DPR_BB_phase
This code use the 2ADPR GPM products (available on the NASA platform https://storm.pps.eosdis.nasa.gov/storm/) to automatically retrieve the bright band elevations (height, bottom, and top).

## How to use
The code is developed for python 3.8+ and should run on any recent Windows system (and most likely also Linux, but not tested).

The following python packages are required:
  * numpy
  * matplotlib
  * proplot
  * h5py (to read the HDF files)
  * pandas
  * os
  * re

Other useful functions are available in the funcDPR.py file.

## What does this code do?
The code DPR_melting_layer.py iterates over all the HDF files placed in the 'Data' directory.
It calculates the elevation of the bright band (height, bottom and top) over the region specified by the user. It computes the median values, standard deviation (and other statistics) of the DPR variables inside a circular area defined by the user.
In the example provided, calculations are performed in a radius of 0.5° (about 50km) centered in Grenoble (France).
The lists of elevations retrieved are exported in a CSV file in the 'Outputs' directory.
Several visualizations of the spatialized variables are done (bright band elevation, Vertical profile, phase near the surface, etc.) and more can be done by adapting this code to other regions with other variables.
Plots are automatically exported in the 'Plots' directory.

## Data
A sample of 2ADPR product is provided in the 'Data/Alps/' folder to test the code. More data are opnely available after previous registration on the NASA platform (https://storm.pps.eosdis.nasa.gov/storm/).

## Questions
In case of any questions, please don't hesitate to contact Arnaud Reboud: arnaud [dot] reboud [at] univ [dash] grenoble [dash] alpes [dot] fr
