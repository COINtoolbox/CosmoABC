import numpy as np
from numpy.random import normal
from scipy.stats import uniform

def my_simulation(v):
    """
    Toy model simulator.
    Samples a normally distributed random variable 
    v['n'] times, having  mean =  v['mean'] and 
    variance = v['std']. 

    input: v -> dictionary of input parameters

    output: scalar 
    """

    l1 = normal(loc=v['mean'],
                scale=v['std'],
                size=v['n'])

    return np.atleast_2d(l1).T


def my_prior(par, func=False):
    """
    Flat prior.
    If func=False returns a random number beteween 
    par[0] and par[1]. 

    Otherwise, returns the corresponding uniform 
    distribution function. 

    input: par -> dictionary of input parameters
    
    output: scalar (if func=False)
            uniform function (if func=True)
    """

    gap = par['max'] - par['min']
    pdf = uniform(loc=par['min'], scale=gap)

    if func == False :
        draw = pdf.rvs()
        return draw
    else:
        return pdf

def my_distance(d2, p):
    """
    Distance between observed and simulated catalogues. 

    input: d2 -> array of simulated catalogue
           p -> dictonary of input parameters

    output: list of 1 scalar (distance)
    """

    mean_obs = np.mean(p['dataset1'])
    std_obs = np.std(p['dataset1'])

    gmean = abs((mean_obs - np.mean(d2))/mean_obs)
    gstd = abs((std_obs - np.std(d2))/std_obs)

    rho = gmean + gstd


    return np.atleast_1d(rho)
