"""
Authors:  Andrew Mummery and Sjoert van Velzen. 

Functions for loading and plotting derived quantities from the paper. 

Loading: 
    load_dicts() loads the raw dictionaries of derived quantities from the paper. 

    get_pars(k) returns numpy arrays of a specific (labelled by k) parameter from the paper. 

Plotting:
    plot_me(kx, ky) plots the variables with labels kx and ky against each other, 
    and includes the options of running analysis scripts from the paper.  

main() shows some examples from the paper. 

"""


import numpy as np 
import matplotlib.pyplot as plt
import pickle 
import corner 


import plot_utils## plot styling 
import analysis_functions## analysis functions from the paper 

def main():
    
    # plot_utils.format_plots_like_paper() # uncomment for fancier plotting (but also slower due to useTex)
    
    ## Plot galaxy mass versus black hole mass measured form the TDE plateau luminosity. 
    kx = 'Galaxy mass'## X-variable
    ky = 'Plateau black hole mass'## Y-variable. 

    f = plot_me(kx, ky, 
                       do_fit=True, ### Also fit a straight line, and report kendall tau statistics. 
                       do_mcmc=True, ## Run an mcmc chain 
                       quiet=False, 
                       mcmc_scale=6, ## Normalise plateau black hole mass by 10^6 M_sun. 
                       nstep_mcmc=10000)


    ## Or you might want to learn about all the TDEs in the sample. 
    n = get_pars('Name')
    ra = get_pars('RA')
    dec = get_pars('DEC')
    z = get_pars('Redshift')
    sp = get_pars('Spectral type')

    print('\n In this sample we have the following TDEs, with parameters', 
          '\n', 
           'Name    RA     DEC     Redshift    Spectral type',
           '\n')
    
    for i  in range(len(n)):
        str_ = "%s     %s     %s     %s     %s \n"%(n[i], ra[i], dec[i], z[i], sp[i])
        print(str_)


    ## Or you might want to do your own analsis.  
    x, x_err_low, x_err_up = get_pars(kx)
    y, y_err_low, y_err_up = get_pars(ky)

    fig = plt.figure()
    ax = fig.add_subplot()

    i_plot = (x>0) * (y>0)# only plot TDEs with relevant data. 

    ax.plot(x[i_plot], y[i_plot], 'o')

    plt.show()

    plt.pause(0.2)    
    key = input('done (press any key to exit)')




def get_pars(k):
    """
    Inputs: 
        k = label of the variable you want from the Table in the paper. 
    
    Returns: 
        vals = the values of the variable. 
        err_low = the lower error on the variable. 
        err_up = the upper error on the variable. 
    
    Notes:
        Quantities returned in log_10, and errors represent the 68% confidence level interval. 

        If k is one of 'Name', 'RA', 'DEC', 'Redshift', 'Spectral type', then no errors returned. 

    """

    dp, dn = load_dicts()

    if k in ['Name', 'RA', 'DEC', 'Redshift', 'Spectral type']:
        vals = []
        for n in dp[k]:
            vals += [n]
        for n in dn[k]:
            vals += [n]
        return vals


    val = np.asarray([dp[k][i][0] for i in range(len(dp[k]))])
    err_down = np.asarray([dp[k][i][1] for i in range(len(dp[k]))])
    err_up = np.asarray([dp[k][i][2] for i in range(len(dp[k]))])
    
    tmp = np.asarray([0 for _ in range(len(dn[k]))])
    tmp_err_down = np.asarray([0 for _ in range(len(dn[k]))])
    tmp_err_up = np.asarray([0 for _ in range(len(dn[k]))])

    if ('Plateau' not in k): 
        tmp = np.asarray([dn[k][i][0] for i in range(len(dn[k]))])
        tmp_err_down = np.asarray([dn[k][i][1] for i in range(len(dn[k]))])
        tmp_err_up = np.asarray([dn[k][i][2] for i in range(len(dn[k]))])

    val = np.append(val, tmp)
    err_down = np.append(err_down, tmp_err_down)
    err_up = np.append(err_up, tmp_err_up)

    return val, err_down, err_up


def load_dicts():
    with open('data/inferred_params/plateau_tdes.pkl', 'rb') as f:
        dp = pickle.load(f)    
    with open('data/inferred_params/no_plateau_tdes.pkl', 'rb') as f:
        dn = pickle.load(f)    
    return dp, dn

def plot_me(kx, ky, 
            quiet=True, 

            lx=None, 
            ly=None, 
            cp='darkblue', 
            cn='red',

            do_fit=False, 
            plot_fit=None, 

            do_mcmc=False,
            plot_mcmc_corner=None,
            plot_mcmc_samples=None,
            plot_mcmc_walkers=None,
            nstep_mcmc=10000,
            nwalkers_mcmc=32,
            mcmc_scale=0):
    
    """
    Inputs: 
        kx = label of the x-variable you want from the Table in the paper. 
        ky = label of the y-variable you want from the Table in the paper. 

        lx = x-axis label for plots, if None then uses kx.
        ly = y-axis label for plots, if None then uses ky.

        cp = color of markers for TDEs with plateaus. 
        cn = color of markers for TDEs without plateaus. 

        do_fit = run analys_functions.simple_fit() on the data. 
        plot_fit = plot best fitting straight line to the data.  If None then set to do_fit. 

        do_mcmc = run analys_functions.mcmc_fit() on the data. 
        plot_mcmc_corner = corner plot of MCMC results. If None then set to do_mcmc. 
        plot_mcmc_samples = plot some samples of MCMC results on the data. If None then set to do_mcmc. 
        plot_mcmc_walkers = show walker evolution from MCMC results. If None then set to do_mcmc. 

        nstep_mcmc, nwalkers_mcmc and mcmc_scale are parameters from analys_functions.mcmc_fit(). 

        quiet = passed to analys_functions.simple_fit(). If False prints analysis results.  

    
    Returns: 
        Up to 3 figures (depending on settings), of various results for parameters x vs y. 


    """

    if plot_fit is None:
        plot_fit = do_fit
    
    if plot_mcmc_corner is None:
        plot_mcmc_corner = do_mcmc
    if plot_mcmc_samples is None:
        plot_mcmc_samples = do_mcmc
    if plot_mcmc_walkers is None:
        plot_mcmc_walkers = do_mcmc
    

    dp, dn = load_dicts()

    if (kx not in list(dp.keys())) or (ky not in list(dp.keys())):
        if (kx not in list(dp.keys())) and (ky not in list(dp.keys())):
            raise ValueError('Both of your chosen variables "%s", "%s", are not allowed. Allowed values are: %s'%(kx, ky, [*dp.keys()]))
        elif (kx not in list(dp.keys())):
            raise ValueError('Your chosen variable "%s" is not allowed. Allowed values are: %s'%(kx, [*dp.keys()]))
        elif (ky not in list(dp.keys())):
            raise ValueError('Your chosen variable "%s" is not allowed. Allowed values are: %s'%(ky, [*dp.keys()]))
        return 

    if lx is None:
        lx = kx 
        tmp_a=list(dp.keys())
        tmp_b=[kx for _ in range(len(dp['Units']))]
        tmp_i=0
        while tmp_a[tmp_i] != tmp_b[tmp_i]:
            tmp_i+=1 
        lx += ' (' + dp['Units'][tmp_i] + ')'
    if ly is None:
        ly = ky
        tmp_a=list(dp.keys())
        tmp_b=[ky for _ in range(len(dp['Units']))]
        tmp_i=0
        while tmp_a[tmp_i] != tmp_b[tmp_i]:
            tmp_i+=1 
        ly += ' (' + dp['Units'][tmp_i] + ')'

    fig = plt.figure()
    ax = fig.add_subplot()

    dx_p = dp[kx]
    dy_p = dp[ky]

    _min = 1e10
    _max = 0

    for i in range(len(dp['Name'])):
        if (dx_p[i][0] > 0) and (dy_p[i][0] > 0):
            ax.errorbar(dx_p[i][0], 
                        dy_p[i][0], 
                        xerr=np.array([[dx_p[i][1], dx_p[i][2]]]).T, 
                        yerr=np.array([[dy_p[i][1], dy_p[i][2]]]).T, 
                        c=cp, 
                        fmt='o'                    
            )
            if dx_p[i][0]<_min:
                _min=dx_p[i][0]
            if dx_p[i][0]>_max:
                _max = dx_p[i][0]


    if ('Plateau' not in kx) and ('Plateau' not in ky): 
        dx_n = dn[kx]
        dy_n = dn[ky]
        for i in range(len(dn['Name'])):
            if (dx_n[i][0] > 0) and (dy_n[i][0] > 0):
                ax.errorbar(dx_n[i][0], 
                            dy_n[i][0], 
                            xerr=np.array([[dx_n[i][1], dx_n[i][2]]]).T, 
                            yerr=np.array([[dy_n[i][1], dy_n[i][2]]]).T, 
                            c=cn, 
                            fmt='s'                    
                )
                if dx_n[i][0]<_min:
                    _min=dx_n[i][0]
                if dx_n[i][0]>_max:
                    _max = dx_n[i][0]


    ax.set_xlabel(lx)
    ax.set_ylabel(ly)
    ax.grid(True)

    if do_fit:
        pars=analysis_functions.simple_fit(kx, ky, dp, dn, quiet)
        if plot_fit:
            xx = np.linspace(_min-0.1*(_max-_min), _max+0.1*(_max-_min))
            ax.plot(xx, xx*pars[0]+pars[1], c='k', ls='--')

    if do_mcmc:
        samples, flat_samples = analysis_functions.mcmc_fit(kx, ky, dp, dn, quiet=quiet, nstep=nstep_mcmc, mcmc_scale=mcmc_scale, n_walkers=nwalkers_mcmc)
        if plot_mcmc_samples:
            xx = np.linspace(_min-0.1*(_max-_min), _max+0.1*(_max-_min))
            for j in range(100):
                a, b, _ = flat_samples[-j-1]
                ax.plot(xx, xx*a+b-a*mcmc_scale, c='lightgreen', alpha=0.1, zorder=-100)
        if plot_mcmc_corner:
            fig2 = plt.figure()
            fig2 = corner.corner(flat_samples, fig=fig2, plot_contours=True, color='b')
            axes = fig2.get_axes()

            axes[3].set_ylabel(r'$\alpha$')
            axes[6].set_ylabel(r'$\epsilon$')
            axes[6].set_xlabel(r'$\beta$')
            axes[7].set_xlabel(r'$\alpha$')
            axes[8].set_xlabel(r'$\epsilon$')

            axes[6].set_ylim(max([0, axes[6].get_ylim()[0]]))
            axes[7].set_ylim(max([0, axes[7].get_ylim()[0]]))
            axes[8].set_xlim(max([0, axes[8].get_xlim()[0]]))
            if plot_mcmc_walkers:
                fig3, axes3 = plt.subplots(3, sharex=True)
                labels = [r"$\alpha$", r"$\beta$", r"$\epsilon$"]
                for i in range(3):
                    ax_ = axes3[i]
                    ax_.plot(samples[:, :, i], "k", alpha=0.3)
                    ax_.set_xlim(0, nstep_mcmc)
                    ax_.set_ylabel(labels[i])
                    ax_.yaxis.set_label_coords(-0.1, 0.5)

                axes3[-1].set_xlabel("step number")


                return fig, fig2, fig3
            return fig, fig2
    return fig


if __name__ == "__main__":
    main()
