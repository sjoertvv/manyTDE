"""
Authors:  Andrew Mummery. 

This script provides functions for loading in the Greene et al. 2020 TDE sample.  

The function main() shows an example of how to use the various functions in this script.  
main() also loads in the new TDE black hole masses and plots both samples together.  

"""


import numpy as np 
import matplotlib.pyplot as plt
import pickle 
import plot_utils
import read_catalog_and_plot


def main():
    ## Load Greene et al. 2020 black hole masses and velocity dispersions.  
    s, err_s_low, err_s_up, m, err_m_low, err_m_up = get_greene_m_sigmas()
    
    ## Note that Greene et al. provide no errors on Galaxy mass measurements. 
    ## Our scipt simply returns 0 for uniformity of syntax. 
    g, err_g_low, err_g_up, mg, err_mg_low, err_mg_up = get_greene_m_galaxy()

    ## Load the TDE measurements from the paper. 
    s_tde, s_tde_err_low, s_tde_err_up = read_catalog_and_plot.get_pars('Velocity dispersion')
    g_tde, g_tde_err_low, g_tde_err_up = read_catalog_and_plot.get_pars('Galaxy mass')
    m_tde, m_tde_err_low, m_tde_err_up = read_catalog_and_plot.get_pars('Plateau black hole mass')


    i_plot_s = (s_tde>0)*(m_tde>0)## Only plot TDEs with relevant data. 
    plt.scatter(s, m)
    plt.scatter(s_tde[i_plot_s], m_tde[i_plot_s])

    plt.figure()
    i_plot_g = (g_tde>0)*(m_tde>0)## Only plot TDEs with relevant data. 
    plt.scatter(g, mg)
    plt.scatter(g_tde[i_plot_g], m_tde[i_plot_g])

    plt.show()



def load_greene_dicts():
    """ Loads in dictionaries from stored pickle files. 
    
        Dictionaries also contain galaxy names for each blck hole mass, 
        these names are not returned in the functions 'get_greene_m_X'. 
    """
    with open('data/dynamical/greene_sigmas.pkl', 'rb') as f:
        dgs = pickle.load(f)    
    with open('data/dynamical/greene_galaxy_masses.pkl', 'rb') as f:
        dgg = pickle.load(f)    
    return dgs, dgg

def get_greene_m_sigmas():
    """  Returns 6 numpy arrays. In returned order these are 
    
        1. velocity dispersion, 
        2. lower error on velocity dispersion, 
        3. upper error on velocity disperion, 
        4. black hole mass, 
        5. lower error on black hole mass, 
        6. upper error on black hole mass. 

        Each returned quantity is measured in log_10.  

        These values are taken from Greene et al. 2020, and represent the values in their tables 
        which appear in their figures.   
    """
    dgs, _ = load_greene_dicts()

    s_val = np.asarray([dgs['Velocity dispersion'][i][0] for i in range(len(dgs['Velocity dispersion']))])
    s_err_down = np.asarray([dgs['Velocity dispersion'][i][1] for i in range(len(dgs['Velocity dispersion']))])
    s_err_up = np.asarray([dgs['Velocity dispersion'][i][2] for i in range(len(dgs['Velocity dispersion']))])

    m_val = np.asarray([dgs['Black hole mass'][i][0] for i in range(len(dgs['Black hole mass']))])
    m_err_down = np.asarray([dgs['Black hole mass'][i][1] for i in range(len(dgs['Black hole mass']))])
    m_err_up = np.asarray([dgs['Black hole mass'][i][2] for i in range(len(dgs['Black hole mass']))])
    
    return s_val, s_err_down, s_err_up, m_val, m_err_down, m_err_up

def get_greene_m_galaxy():
    """  Returns 6 numpy arrays. In returned order these are 
    
        1. Galaxy mass, 
        2. lower error on galaxy mass, 
        3. upper error on galaxy mass, 
        4. black hole mass, 
        5. lower error on black hole mass, 
        6. upper error on black hole mass. 

        Each returned quantity is measured in log_10.  

        These values are taken from Greene et al. 2020, and represent the values in their tables 
        which appear in their figures.   
    """

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