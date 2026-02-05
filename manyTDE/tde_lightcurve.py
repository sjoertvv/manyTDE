"""
Authors: Sjoert van Velzen and Andrew Mummery. 

Provides a python function "plot_tde" which gets the lightcurve data and plots. 
Second function "get_lightcurve_data" returns the data in various filters.  

The data is returned/plotted at full time resolution (no binning). 

This example also shows how to correct the TDE light for extinction. 
"""

import json
import matplotlib.pyplot as plt 
import numpy as np 

import plot_utils

from pathlib import Path
we_are_here = Path(__file__).resolve().parent # note, this only work as path to data if we install with pip -e

def main():
	
	#plot_utils.format_plots_like_paper() # toggle this to match the formatting of the paper
	
	"""
	An example using the TDE AT2019dsg.  
	"""

	# Make a Figure
	_ = plot_tde('AT2019dsg')


	### Or you might want to get the data yourself ...
	lc_dict, filters, frequency_Hz = get_lightcurve_data('AT2019dsg')

	#  ... and make your own plots
	# plt.figure()
	# plt.errorbar(lc_dict[filters[-1]][0], lc_dict[filters[-1]][1], yerr=lc_dict[filters[-1]][2], fmt='o')


	plt.pause(0.2)    
	key = input('done (press any key to exit)')

	# plt.show()



def plot_tde(tde_name='ASASSN-14li'):
	"""
	Plot the lightcurve of the tde with name tde_name. Returns the lightcurves as a figure.  
	"""
	fname = we_are_here / "data" / "sources" / "{0}.json".format(tde_name)
	tde_data = json.load(open(fname,'r'))# Load data. 

	# These conversion are needed because json doesn't store tuples.
	dt = [tuple(x) for x in tde_data['lightcurve']['dtype']]
	lc_obj = [tuple(x) for x in tde_data['lightcurve']['data']] 

	# Make a recarray. 
	lc_rec = np.array(lc_obj, dtype=dt)
	mjd0 = tde_data['peak_mjd']

	fig = plt.figure()
	ax = fig.add_subplot()

	for flt in tde_data['lightcurve']['filters']:
		idx = lc_rec['filter']==flt

		flux = lc_rec[idx]['flux_Jy']*1e6
		flux_corr = flux / tde_data['extinction']['linear_extinction'][flt]# Correct for extinction. 

		ax.errorbar(lc_rec[idx]['mjd']-mjd0, 
				flux_corr, 
				lc_rec[idx]['e_flux_Jy']*1e6,
				fmt=plot_utils.marker_dict[flt], 
				alpha=0.9,
				color=plot_utils.lc_color_dict[flt],
				label=flt)

	ax.set_title(tde_name)
	ax.legend()
	ax.set_xlabel('MJD-{0:0.1f}'.format(mjd0))
	ax.set_ylabel(r'Flux ($\mu$Jy)')
	
	return fig

def get_lightcurve_data(tde_name = 'ASASSN-14li'):
	"""
	Input: 
		The TDEs name

	Returns:
		1. A dictionary with all of the light curve data, labelled by observing band. 
		2. A list of lightcurve filters with available data. 
	"""

	fname = we_are_here / "data" / "sources" / "{0}.json".format(tde_name)
	tde_data = json.load(open(fname,'r'))# Load data. 

	# These conversion are needed because json doesn't store tuples.
	dt = [tuple(x) for x in tde_data['lightcurve']['dtype']]
	lc_obj = [tuple(x) for x in tde_data['lightcurve']['data']] 

	# Make a recarray. 
	lc_rec = np.array(lc_obj, dtype=dt)
	mjd0 = tde_data['peak_mjd']

	lc_dict = {}
	filters = tde_data['lightcurve']['filters']
	frequency_Hz = tde_data['lightcurve']['frequency_Hz']

	for flt in filters:
		idx = lc_rec['filter']==flt

		flux = lc_rec[idx]['flux_Jy']*1e6
		flux_corr = flux / tde_data['extinction']['linear_extinction'][flt]# Correct for extinction. 

		lc_dict[flt] = [lc_rec[idx]['mjd']-mjd0, flux_corr, lc_rec[idx]['e_flux_Jy']*1e6]
	
	return lc_dict, filters, frequency_Hz


if __name__ == "__main__":
	main()

