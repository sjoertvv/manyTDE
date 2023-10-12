import numpy as np 
import matplotlib.pyplot as plt
import pickle 
import emcee 
import corner 
from scipy.stats import linregress, kendalltau
from scipy.optimize import minimize 
import plot_utils
import read_catalog_and_plot


def main():
    s, esl, esu, m, msl, msu = get_greene_m_sigmas()
    g, egl, egu, mg, mgsl, mgsu = get_greene_m_galaxy()

    stde, stdeel, stdeeu = read_catalog_and_plot.get_pars('Velocity dispersion')
    gtde, gtdeel, gtdeeu = read_catalog_and_plot.get_pars('Galaxy mass')
    mtde, mtdeel, mtdeeu = read_catalog_and_plot.get_pars('Plateau black hole mass')

    plt.scatter(s, m)
    plt.scatter(stde[(stde>0)*(mtde>0)], mtde[(stde>0)*(mtde>0)])

    plt.figure()
    plt.scatter(g, mg)
    plt.scatter(gtde[(gtde>0)*(mtde>0)], mtde[(gtde>0)*(mtde>0)])

    plt.show()



def load_greene_dicts():
    with open('data/dynamical/greene_sigmas.pkl', 'rb') as f:
        dgs = pickle.load(f)    
    with open('data/dynamical/greene_galaxy_masses.pkl', 'rb') as f:
        dgg = pickle.load(f)    
    return dgs, dgg

def get_greene_m_sigmas():
    dgs, _ = load_greene_dicts()

    s_val = np.asarray([dgs['Velocity dispersion'][i][0] for i in range(len(dgs['Velocity dispersion']))])
    s_err_down = np.asarray([dgs['Velocity dispersion'][i][1] for i in range(len(dgs['Velocity dispersion']))])
    s_err_up = np.asarray([dgs['Velocity dispersion'][i][2] for i in range(len(dgs['Velocity dispersion']))])

    m_val = np.asarray([dgs['Black hole mass'][i][0] for i in range(len(dgs['Black hole mass']))])
    m_err_down = np.asarray([dgs['Black hole mass'][i][1] for i in range(len(dgs['Black hole mass']))])
    m_err_up = np.asarray([dgs['Black hole mass'][i][2] for i in range(len(dgs['Black hole mass']))])
    
    return s_val, s_err_down, s_err_up, m_val, m_err_down, m_err_up

def get_greene_m_galaxy():
    _, dgg = load_greene_dicts()

    g_val = np.asarray([dgg['Galaxy mass'][i][0] for i in range(len(dgg['Galaxy mass']))])
    g_err_down = np.asarray([dgg['Galaxy mass'][i][1] for i in range(len(dgg['Galaxy mass']))])
    g_err_up = np.asarray([dgg['Galaxy mass'][i][2] for i in range(len(dgg['Galaxy mass']))])

    m_val = np.asarray([dgg['Black hole mass'][i][0] for i in range(len(dgg['Black hole mass']))])
    m_err_down = np.asarray([dgg['Black hole mass'][i][1] for i in range(len(dgg['Black hole mass']))])
    m_err_up = np.asarray([dgg['Black hole mass'][i][2] for i in range(len(dgg['Black hole mass']))])
    
    return g_val, g_err_down, g_err_up, m_val, m_err_down, m_err_up

if __name__ == "__main__":
    main()