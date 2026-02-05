"""
Authors:  Andrew Mummery and Sjoert van Velzen. 

Two functions which reproduce analysis from the paper.  

simple_fit() fits a straight line to the data, and reports information about the Kendall tau statistic and scatter. 

mcmc_fit() runs the MCMC analysis from the paper.  

Both functions are described in more detail below. 
"""


import numpy as np
import emcee 
from scipy.stats import linregress, kendalltau
from scipy.optimize import minimize 
import warnings## to deal with scipy nonsense. 

def simple_fit(kx, ky, dp, dn,
               quiet=True):
    """
    Inputs: 

        kx = label of x-variable. 
        ky = label of y-variable. 
        dp = dictionary of TDEs with plateaus. 
        dn = dictionary of TDEs without plateaus. 
        quiet = whether to print values (False) or not (True). 

    Note: both dp and dn are returned from the function
    dp, dn =  read_catalog_and_plot.load_dicts() 

    Returns: 

        fit_pars = a list of parameters from simple fitting between the x and y variables.  

        fit_pars is of the form [m, c, p, ktv, ktp, scatter], where 
        m, c -> parameters from a straight line fit of the form y = mx + c. 
        ktv, ktp -> Kendall tau parameter (ktv) and p-value (ktp). 
        scatter -> average scatter (in dex) between y-values and best fit straight line. 

    """
    
    dx_p = dp[kx]
    dy_p = dp[ky]

    j = 0
    for i in range(len(dp['Name'])):
        if (dx_p[i][0] > 0) and (dy_p[i][0] > 0):
            if j == 0:
                data_x = [dx_p[i][0]]
                data_y = [dy_p[i][0]]
                j+=1
            else:
                data_x += [dx_p[i][0]]
                data_y += [dy_p[i][0]]


    if ('Plateau' not in kx) and ('Plateau' not in ky): 
        dx_n = dn[kx]
        dy_n = dn[ky]
        for i in range(len(dn['Name'])):
            if (dx_n[i][0] > 0) and (dy_n[i][0] > 0):
                data_x += [dx_n[i][0]]
                data_y += [dy_n[i][0]]
    
    data_x = np.asarray(data_x)
    data_y = np.asarray(data_y)

    fit = linregress(data_x, data_y)
    if not quiet:
        print('\n Variable x = %s, variable y = %s'%(kx, ky))
        print('\n', 'Linear regression slope = ', fit.slope, '±',  fit.stderr, '\n', 'Linear regression p-value = ', fit.pvalue)

    m, c, p = fit.slope, fit.intercept, fit.pvalue

    kt = kendalltau(data_x, data_y)
    if not quiet:
        print('\n Kendall tau statistic = ',kt[0],'\n','Kendall tau p-value = ', kt[1])
    ktv, ktp = kt[0], kt[1]

    
    scatter= np.std(data_y-(m*data_x+c))
    if not quiet:
        print ('\n rsm y-model(x):', scatter, '\n')


    return [m, c, p, ktv, ktp, scatter]



def mcmc_fit(kx, ky, dp, dn,
             quiet = True, 
             nstep = 10000,
             mcmc_scale=0,
             n_walkers=32,
             f_discard=0.5):
    """
    Inputs: 

        kx = label of x-variable. 
        ky = label of y-variable. 
        dp = dictionary of TDEs with plateaus. 
        dn = dictionary of TDEs without plateaus. 

        nstep = number of steps in MCMC chain. 
        mcmc_scale = log_10 of a scale for the y-value, i.e., MCMC uses Y = data_y / 10**mcmc_scale as input. 
        nwalkers = number of walkers in MCMC chain. 
        f_discard = fraction of steps to discard for burn in. 

        quiet = whether to print results from analysis (False) or not (True). 

    Note: both dp and dn are returned from the function
    dp, dn =  read_catalog_and_plot.load_dicts() 

    Returns: 

        samples = the full MCMC chain output. 
        flat_samples = the flattened and burn in discarded chain, for plotting in a corner plot.  
    
    This MCMC fit uses the likelihood function described in the paper.  
    """
    

    dx_p = dp[kx]
    dy_p = dp[ky]

    j = 0
    for i in range(len(dp['Name'])):
        if (dx_p[i][0] > 0) and (dy_p[i][0] > 0):
            if j == 0:
                data_x = [dx_p[i][0]]
                data_y = [dy_p[i][0]]
                err_y = [np.mean([dy_p[i][1], dy_p[i][2]])]
                j+=1
            else:
                data_x += [dx_p[i][0]]
                data_y += [dy_p[i][0]]
                err_y += [np.mean([dy_p[i][1], dy_p[i][2]])]


    if ('Plateau' not in kx) and ('Plateau' not in ky): 
        dx_n = dn[kx]
        dy_n = dn[ky]
        for i in range(len(dn['Name'])):
            if (dx_n[i][0] > 0) and (dy_n[i][0] > 0):
                data_x += [dx_n[i][0]]
                data_y += [dy_n[i][0]]
                err_y += [np.mean([dy_p[i][1], dy_p[i][2]])]
    
    data_x = np.asarray(data_x)
    data_y = np.asarray(data_y)
    err_y = np.asarray(err_y)

    def log_lklihood(pars, data):## The log likelihood from the paper. 
        a, b, eps = pars[0], pars[1], pars[2]
        x_dat, y_dat, y_err = data[0], data[1], data[2]
        mod = a * x_dat + b - a * mcmc_scale
        sigma2 = y_err**2 + eps**2
        lkl = - np.sum((mod - y_dat)**2/sigma2 + np.log(sigma2))
        return lkl

    def log_prior(pars):## Trivial prior, do not want negative intrinsic scatter. 
        a, b, eps = pars[0], pars[1], pars[2]
        if eps < 0:
            return -np.inf
        return 0

    def log_prob(pars, data):## Probability = Prior * Likelihood. 
        lp = log_prior(pars)
        if not np.isfinite(lp):
            return -np.inf
        return lp + log_lklihood(pars, data)

    def n_log_prob(pars, data):
        return -log_prob(pars, data)

    pars = simple_fit(kx, ky, dp, dn, quiet=True)## Simple fit for initial guess. 
    
    warnings.filterwarnings("ignore")## scipy.minimize() throws meaningless warnings all the time. 
    soln = minimize(n_log_prob, [pars[0], pars[1]+mcmc_scale*pars[0], 1], args=[data_x, data_y, err_y])
    
    pos = soln.x[:3] + 1e-4 * np.random.randn(n_walkers, 3) 
    nwalkers, ndim = pos.shape


    sampler = emcee.EnsembleSampler(
        nwalkers, ndim, log_prob, args=[[data_x, data_y, err_y]])## Set up MCMC sampler. 
    sampler.run_mcmc(pos, nstep, progress=True)## Run chain.  
    
    warnings.filterwarnings("default")## Switch warnings back on now we have done scipy nonsense. 

    n_discard = int(f_discard*nstep)

    samples = sampler.get_chain()
    flat_samples = samples[n_discard:].reshape((-1, ndim))

    if not quiet:
            print('\n Results from MCMC analysis \n', 
            'alpha = ', np.median(flat_samples[:, 1]), ' ± ', np.std(flat_samples[:, 1]), '\n', 
            'beta = ', np.median(flat_samples[:, 0]), ' ± ', np.std(flat_samples[:, 0]), '\n',
            'epsilon = ', np.median(flat_samples[:, 2]), ' ± ', np.std(flat_samples[:, 2]), '\n',
            )


    return samples, flat_samples


