# manyTDE
Collection of 63 optically-selected TDEs, with black hole mass measurements based on their late-time plateau luminosity. 

## Paper
This repository containts the data that was used in [Mummery, van Velzen et al. 2023](https://ui.adsabs.harvard.edu/abs/2023arXiv230808255M/abstract). Please cite this paper if you use these products. 

## Catalog
The main catalog is stored in a Python pickle file. A script to read and reproduce some key figures from the paper is provided. 

`>> python3 read_catalog_and_plot.py`

This catalog contains: 
-  the basic info (name, coordinates, redshift);
-  host galaxy info (eg, mass, velocity dispersion);
-  TDE spectral type (if available)
-  the features of the TDE lightcurve (eg, peak luminosity, decay time, plateau luminosity);
-  our estimates of the black hole mass (eg, based on the plateau luminosity).

The full list of catalog columns is described below.  

## Lightcurves
For all TDEs we also provide the optical/UV light curves in the `/data/lightcurves/` folder. A script to plot these data is also provided. 

`>> python3 plot_example_lightcurve.py`


## Dynamical black hole mass measurements 
The dynamical black hole mass measurements from the review by [Greene, Strader, & Ho (2020)](https://ui.adsabs.harvard.edu/abs/2020ARA%26A..58..257G/abstract) are stored in `data/dynamical/`

***

## Full catalog describtion

Mosf of the parameters below are measured in log10; first entry is the median of the posterior, the second, third and the 68% CL interval (lower, upper).

```
gal_mass        : galaxy stellar mass, estimated from the host galaxy photometry (Msun)
sigma           : velocity dispersion, if available (km/s)
plat_lum_exp    : plateau nu*Lnu at nu_kc when using exponential for early-part (erg/s)
plat_lum_exp_gr : plateau nu*Lnu at nu_gr when using exponential for early-part (erg/s)
nu_kc           : rest-frame frequency used for model curves (Hz), except for *_gr
nu_gr           : "optical" frequency used for plat_lum_exp_gr (6e14 Hz)
...
(to do, finish this documentation)
```
