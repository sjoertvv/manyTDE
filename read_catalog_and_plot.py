import numpy as np 
import matplotlib.pyplot as plt
import pickle 
import emcee 
import corner 
from scipy.stats import linregress, kendalltau
from scipy.optimize import minimize 
import plot_utils

def main():

    k1 = 'Galaxy mass'
    k2 = 'Plateau black hole mass'

    f = plot_me(k1, k2, 
                       do_fit=True, 
                       plot_fit=True, 
                       do_mcmc=True, 
                       plot_mcmc_samples=True, 
                       plot_mcmc_corner=True, 
                       plot_mcmc_walkers=True,
                       quiet=False, 
                       mcmc_scale=11,
                       nstep_mcmc=10000)


    x, _, _ = get_pars(k1)
    y, _, _ = get_pars(k2)
    n = get_pars('Name')
    ra = get_pars('RA')
    dec = get_pars('DEC')
    z = get_pars('Redshift')
    sp = get_pars('Spectral type')

    for i  in range(len(n)):
        print(n[i], ra[i], dec[i], z[i], sp[i])


    fig = plt.figure()
    ax = fig.add_subplot()

    ax.plot(x, y, 'o')

    plt.show()



def get_pars(k1):
    dp, dn = load_dicts()

    if k1 in ['Name', 'RA', 'DEC', 'Redshift', 'Spectral type']:
        vals = []
        for n in dp[k1]:
            vals += [n]
        for n in dn[k1]:
            vals += [n]
        return vals


    val = np.asarray([dp[k1][i][0] for i in range(len(dp[k1]))])
    err_down = np.asarray([dp[k1][i][1] for i in range(len(dp[k1]))])
    err_up = np.asarray([dp[k1][i][2] for i in range(len(dp[k1]))])
    
    tmp = np.asarray([0 for _ in range(len(dn[k1]))])
    tmp_err_down = np.asarray([0 for _ in range(len(dn[k1]))])
    tmp_err_up = np.asarray([0 for _ in range(len(dn[k1]))])

    if ('Plateau' not in k1): 
        tmp = np.asarray([dn[k1][i][0] for i in range(len(dn[k1]))])
        tmp_err_down = np.asarray([dn[k1][i][1] for i in range(len(dn[k1]))])
        tmp_err_up = np.asarray([dn[k1][i][2] for i in range(len(dn[k1]))])

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
            lx=None, 
            ly=None, 
            cp='darkblue', 
            cn='red',
            do_fit=False, 
            plot_fit=False, 
            quiet=True, 
            do_mcmc=False,
            plot_mcmc_corner=True,
            plot_mcmc_samples=False,
            plot_mcmc_walkers=True,
            nstep_mcmc=10000,
            nwalkers_mcmc=32,
            mcmc_scale=0):
    
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
        pars=simple_fit(kx, ky, quiet)
        if plot_fit:
            xx = np.linspace(_min-0.1*(_max-_min), _max+0.1*(_max-_min))
            ax.plot(xx, xx*pars[0]+pars[1], c='k', ls='--')

    if do_mcmc:
        samples, flat_samples = mcmc_fit(kx, ky, nstep=nstep_mcmc, mcmc_scale=mcmc_scale, n_walkers=nwalkers_mcmc)
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


def simple_fit(kx, ky, 
               quiet=True):
    
    dp, dn = load_dicts()

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



def mcmc_fit(kx, ky, 
             nstep = 10000,
             mcmc_scale=0,
             n_walkers=32,
             f_discard=0.5):
    
    dp, dn = load_dicts()

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

    def log_lklihood(pars, data):
        a, b, eps = pars[0], pars[1], pars[2]
        x_dat, y_dat, y_err = data[0], data[1], data[2]
        mod = a * x_dat + b - a * mcmc_scale
        sigma2 = y_err**2 + eps**2
        lkl = - np.sum((mod - y_dat)**2/sigma2 + np.log(sigma2))
        return lkl

    def log_prior(pars):
        a, b, eps = pars[0], pars[1], pars[2]
        if eps < 0:
            return -np.inf
        return 0

    def log_prob(pars, data):
        lp = log_prior(pars)
        if not np.isfinite(lp):
            return -np.inf
        return lp + log_lklihood(pars, data)

    def n_log_prob(pars, data):
        return -log_prob(pars, data)

    pars = simple_fit(kx, ky, quiet=True)
    
    import warnings
    warnings.filterwarnings("ignore")
    soln = minimize(n_log_prob, [pars[0], pars[1]+mcmc_scale*pars[0], 1], args=[data_x, data_y, err_y])
    
    pos = soln.x[:3] + 1e-4 * np.random.randn(n_walkers, 3) 
    nwalkers, ndim = pos.shape


    sampler = emcee.EnsembleSampler(
        nwalkers, ndim, log_prob, args=[[data_x, data_y, err_y]])
    sampler.run_mcmc(pos, nstep, progress=True)
    warnings.filterwarnings("default")

    n_discard = int(f_discard*nstep)

    samples = sampler.get_chain()
    flat_samples = samples[n_discard:].reshape((-1, ndim))

    return samples, flat_samples




if __name__ == "__main__":
    main()