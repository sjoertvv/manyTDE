# manyTDE
Collection of 63 optically-selected TDEs, with black hole mass measurements based on their light curve properties. 

## Paper
This repository containts the data that was used in [Mummery, van Velzen et al. 2023](https://ui.adsabs.harvard.edu/abs/2023arXiv230808255M/abstract). Please cite this paper if you use these products. 

## Catalog
The main catalog is stored in a Python pickle file. A script to read and reproduce some key figures from the paper is provided. 
We split the catalogs into those TDEs with and without measured plateaus for convinience. 

`>> python3 read_catalog_and_plot.py`

This catalog contains: 
-  the basic info (name, coordinates, redshift);
-  host galaxy info (eg, mass, velocity dispersion);
-  TDE spectral type (if available)
-  the features of the TDE lightcurve (eg, peak luminosity, decay time, plateau luminosity);
-  our estimates of the black hole mass (eg, based on the plateau luminosity).

The full list of catalog columns is described below.  

## Lightcurves
For all TDEs we also provide the optical/UV light curves in the `/data/lightcurves/` folder. A script to get or plot these data is also provided. 

`>> python3 tde_lightcurve.py`


## Dynamical black hole mass measurements 
The dynamical black hole mass measurements from the review by [Greene, Strader, & Ho (2020)](https://ui.adsabs.harvard.edu/abs/2020ARA%26A..58..257G/abstract) are stored in `data/dynamical/greene_sigmas.pkl` and `data/dynamical/greene_galaxy_masses.pkl`. 

Scripts for loading and plotting these data can be found in

`>> python3 dynamical_and_tdes.py`

***

## Full catalog describtion

All of the physical parameters below are measured in log10; first entry is the median of the posterior. The second and third columns are the 68% confidence level interval, ordered by (lower, upper).

```
Name                                 : TDE name
Redshift                             : Redshift (dimensionless)
RA                                   : RA of source 
DEC                                  : DEC of source
Spectral type                        : TDE spectral type, either H, H+He, He, Featureless or Unknown.
Plateau u-band                       : vL_v of the plateau at v = 10^15 Hz. (Rest frame)
Plateau g-band                       : vL_v of the plateau at v = 6 x 10^14 Hz. (Rest frame)
Peak u-band                          : vL_v at peak at v = 10^15 Hz. (Rest frame)
Peak g-band                          : vL_v at peak at v = 6 x 10^14 Hz. (Rest frame)
Peak bolometric                      : Peak bolometric luminosity (from black body fit).
Peak temperature                     : Temperature (K) from blackbody fit to early times.
Peak radius                          : Radius (cm) from blackbody fit to early times.
Energy radiated g-band               : Energy radiated in the early time g-band light curve.
Rise timescale                       : Rise timescale parameter (days) from early time model.
Decay timescale                      : Exponential decay timescale parameter (days) from early time model.
Galaxy mass                          : Mass of host galaxy from spectroscopic fitting (Solar masses). 
Velocity dispersion                  : Velocity dispersion of host galaxy (km/s).
Plateau black hole mass              : Black hole mass inferred from plateau luminosity (Solar masses). 
Peak black hole mass                 : Black hole mass inferred from peak luminosity (Solar masses). 
Energy radiated  black hole mass     : Black hole mass inferred from green band energy radiated (Solar masses).
Units                                : List of the dimensions of all quantities. 
```
