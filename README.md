# manyTDE
A collection optically-selected TDEs. Including light curves and host galaxy properties (mass, velocity dispersion, etc).

For a subset of these sources, we also include black hole mass measurements (see next sections). 

## Paper
This repository contains black hole mass estimates from [Mummery, van Velzen et al. 2023](https://ui.adsabs.harvard.edu/abs/2023arXiv230808255M/abstract). Please cite this paper if you use these products. 


## Catalogue
The main catalog of black hole mass catalogue is stored in a Python pickle file. A script to read and reproduce some key figures from the paper is provided. 

`>> python3 read_catalog_and_plot.py`

This requires the `corner` and `emcee` packages (install with pip). 

We split the catalogues into those TDEs with and without measured plateaus for convinience. The catalog contains: 
-  the basic info (name, coordinates, redshift);
-  host galaxy info (eg, mass, velocity dispersion);
-  TDE spectral type (if available)
-  the features of the TDE lightcurve (eg, peak luminosity, decay time, plateau luminosity);
-  our estimates of the black hole mass (eg, based on the plateau luminosity).

The full list of catalogue columns is described below.  

## Lightcurves
For all TDEs we also provide the optical/UV light curves in the `/data/sources/` folder. A script to get or plot these data is also provided. 

`>> python3 tde_lightcurve.py`


## Dynamical black hole mass measurements 
The dynamical black hole mass measurements from the review by [Greene, Strader, & Ho (2020)](https://ui.adsabs.harvard.edu/abs/2020ARA%26A..58..257G/abstract) are stored in `data/dynamical/greene_sigmas.pkl` and `data/dynamical/greene_galaxy_masses.pkl`. 

Scripts for loading and plotting these data can be found in

`>> python3 dynamical_and_tdes.py`

***

## Full catalogue describtion

All of the physical parameters below are measured in log10; first entry is the median of the posterior. The second and third columns are the 68% confidence level interval, ordered by (lower, upper).
The first 5 entries (Name, Redshift, RA, DEC and Spectral type) do not have uncertainties and are one dimensional. 

```
Name                                 : TDE name
Redshift                             : Redshift (dimensionless)
RA                                   : RA of source 
DEC                                  : DEC of source
Spectral type                        : TDE spectral type, either H, H+He, He, Featureless or Unknown.
Plateau u-band                       : vL_v of the plateau at v = 10^15 Hz. (Rest frame, erg/s)
Plateau g-band                       : vL_v of the plateau at v = 6 x 10^14 Hz. (Rest frame, erg/s)
Peak u-band                          : vL_v at peak at v = 10^15 Hz. (Rest frame, erg/s)
Peak g-band                          : vL_v at peak at v = 6 x 10^14 Hz. (Rest frame, erg/s)
Peak bolometric                      : Peak bolometric luminosity (from black body fit, erg/s).
Peak temperature                     : Temperature (K) from blackbody fit to early times.
Peak radius                          : Radius (cm) from blackbody fit to early times.
Energy radiated g-band               : Energy radiated in the early time g-band light curve (erg).
Rise timescale                       : Rise timescale parameter (days) from early time model.
Decay timescale                      : Exponential decay timescale parameter (days) from early time model.
Galaxy mass                          : Stellar mass of host galaxy from its photometry (Solar masses). 
Velocity dispersion                  : Velocity dispersion of host galaxy (km/s).
Plateau black hole mass              : Black hole mass inferred from plateau luminosity (Solar masses). 
Peak black hole mass                 : Black hole mass inferred from peak luminosity (Solar masses). 
Energy radiated  black hole mass     : Black hole mass inferred from g-band energy radiated (Solar masses).
Units                                : List of the dimensions of all quantities. 
```
