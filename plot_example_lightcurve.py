'''
Read the lightcurve data and plot

This shows the data at the full time resolution (no binning). 

This example also shows how to correct the TDE light for extinction
'''
import json

from plot_utils import marker_dict, lc_color_dict

tde_name = 'AT2018lna' # aka Cersei
#tde_name = 'AT2019dsg' # aka Bran

fname = './data/lightcurves/{0}.json'.format(tde_name)
tde_data = json.load(open(fname,'r'))

# these conversion are needed because json doesn't store tuples...
dt = [tuple(x) for x in tde_data['lightcurve']['dtype']]
lc_obj = [tuple(x) for x in tde_data['lightcurve']['data']] 

# make a nice recarray
lc_rec = np.array(lc_obj, dtype=dt)
mjd0 = tde_data['peak_mjd']

plt.clf()
for flt in tde_data['lightcurve']['filters']:
	idx = lc_rec['filter']==flt

	flux = lc_rec[idx]['flux_Jy']*1e6
	flux_corr = flux / tde_data['extinction']['linear_extinction'][flt] 

	plt.errorbar(lc_rec[idx]['mjd']-mjd0, 
			flux_corr, 
			lc_rec[idx]['e_flux_Jy']*1e6,
			fmt=marker_dict[flt], 
			alpha=0.9,
			color=lc_color_dict[flt],
			label=flt)

plt.title(tde_name)
plt.legend()
plt.xlabel('MJD-{0:0.1f}'.format(mjd0))
plt.ylabel('Flux (uJy)')
plt.show()
plt.pause(0.1)
key = input('done.')
