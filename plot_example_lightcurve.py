'''
Read the lightcurve data and plot

This shows the data at the full time resolution (no binning). 

This example also shows how to correct the TDE light for extinction
'''
import json
import matplotlib.pyplot as plt 
import numpy as np 

from plot_utils import marker_dict, lc_color_dict

def main():
	# tde_name = 'AT2019dsg' # aka Bran
	_ = plot_tde()
	plt.show()

def plot_tde(tde_name='ASASSN-14li'):
	fname = './data/lightcurves/{0}.json'.format(tde_name)
	tde_data = json.load(open(fname,'r'))

	# these conversion are needed because json doesn't store tuples...
	dt = [tuple(x) for x in tde_data['lightcurve']['dtype']]
	lc_obj = [tuple(x) for x in tde_data['lightcurve']['data']] 

	# make a nice recarray
	lc_rec = np.array(lc_obj, dtype=dt)
	mjd0 = tde_data['peak_mjd']

	fig = plt.figure()
	ax = fig.add_subplot()
	for flt in tde_data['lightcurve']['filters']:
		idx = lc_rec['filter']==flt

		flux = lc_rec[idx]['flux_Jy']*1e6
		flux_corr = flux / tde_data['extinction']['linear_extinction'][flt] 

		ax.errorbar(lc_rec[idx]['mjd']-mjd0, 
				flux_corr, 
				lc_rec[idx]['e_flux_Jy']*1e6,
				fmt=marker_dict[flt], 
				alpha=0.9,
				color=lc_color_dict[flt],
				label=flt)

	ax.set_title(tde_name)
	ax.legend()
	ax.set_xlabel('MJD-{0:0.1f}'.format(mjd0))
	ax.set_ylabel(r'Flux ($\mu$Jy)')
	
	return fig


if __name__ == "__main__":
	main()