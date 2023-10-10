# manyTDE
Collection of optically-selected TDEs, with black hole mass measurements based on their late-time plateau luminosity. 

## Paper
This repository containts the data that was used in [https://ui.adsabs.harvard.edu/abs/2023arXiv230808255M/abstract][Mummery, van Velzen et al. 2023]. Please cite this paper if you use these products. 

## Catalog
The catalog of TDEs is contained in a Python pickle file. An script to read and reproduce some key plots from the paper is provided. 

`>> python3 read_catalog_plot.py`

This catalog contains: 
-  the basic info (name, coordinates, redshift);
-  host galaxy info (eg, mass,  velocity dispersion);
-  the features of the TDE lightcurve (eg, peak luminosity, decay time, plateau luminosity);
-  our estimates of the black hole mass (eg, based on the plateau luminosity).

The full list of columns is described below.  

## Lightcurves
For all TDEs we also provide the optical/UV light curves in the `/data/lightcurves/` folder. A script to plot these data is also provided. 

`>> python3 plot_example_lightcurve.py`

## Dynamical black hole mass measurements 
The dynamical black hole mass measurements from the review by [https://ui.adsabs.harvard.edu/abs/2020ARA%26A..58..257G/abstract][Greene, Strader, & Ho (2020)] are stored in `data/dynamical/`

# Full catalog describtion

Mosf of the keys below are measured in log10; first entry is the median of the posterior, the second, third and the 68% CL interval (lower, upper).

```
gal_mass        : galaxy mass from photometry (Msun)
sigma           : velocity dispersion, from SDSS or literature (km/s)
plat_lum_exp    : plateau nu*Lnu at nu_kc when using exponential for early-part (erg/s)
plat_lum_exp_gr : plateau nuLnu at nu_gr when using exponential for early-part (erg/s)
nu_kc           : rest-frame frequency used for model curves (Hz), except for *_gr
nu_gr           : "optical" frequency used for plat_lum_exp_gr (6e14 Hz)
...
```
