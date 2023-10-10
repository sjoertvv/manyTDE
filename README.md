# manyTDE
Collection of optically-selected TDEs, with black hole mass measurements from their late-time plateau luminosity 

## Paper
This repository containts the data that was used in [Mummery, van Velzen et al. 2023][https://ui.adsabs.harvard.edu/abs/2023arXiv230808255M/abstract]. Please cite this work if you use any of these products. 

## Catalog
The catalog of TDEs is contained in a Python pickle file. An script to read and reproduce some key plots from the paper is provided. 

>> python3 read_catalog_plot.py

This catalog contains the basic info (name, coordinates, redshift), host galaxy info (eg, mass,  velocity dispersion), the features of the TDE lightcurve (eg, peak luminosity, decay time, plateau luminosity), and our estimates of the black hole mass (eg, based on the plateau luminosity). The full list of columns is described below.  

## Lightcurves
For all TDEs we also provide the optical/UV light curves in the `/data/lightcurves/` folder. A script to plot these data is also provided. 

>> python3 plot_example_lightcurve.py

## Dynamical black hole mass measurements 
The dynamical black hole mass measurements from the review by [Greene, Strader, & Ho][https://ui.adsabs.harvard.edu/abs/2020ARA%26A..58..257G/abstract]
