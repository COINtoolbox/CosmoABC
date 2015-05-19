"""
Functions for initial prior distributions. 

***
All functions must have take as input a dictionary of parameters 
and have one optional variable called 'func'
If 'func' is False -> return one sampling of the pdf
If 'func' is True -> return the pdf itself
the dictionary of input parameters must contain 2 keywords:
min and max, which are reserved for delimiting the parameter domain. 
Other variables should have different names.
***

"""

import numpy as np

from scipy.stats import beta
from scipy.stats import uniform
from scipy.stats import norm

def gaussian_prior(par, func=False):
    """
    Gaussian prior.
  
    input: par -> dictionary of parameter values
                  keywords: mean, standard_devitation, 
                            min and max

                  values: all scalars 
           func -> boolean (optional)
                   if True returns the pdf random variable. 
                   Default is False.

    output: scalar (if func=False)
            gaussian probability distribution function (if func=True)

    """

    np.random.seed()    
    dist = norm(loc=par['pmean'], scale=par['pstd'])
    flag = False  
    while flag == False:   
        draw = dist.rvs() 
        if par['min'] < draw and draw < par['max']:
            flag = True
     
    if func == False:
        return draw
    else:
        return dist


def flat_prior(par, func=False):
    """
    Flat prior.
  
    input: par ->  dictionary of parameter values
                   keywords: min and max (scalars)
           func -> boolean (optional)
                   if True returns the pdf random variable. 
                   Default is False.

    output: scalar (if func=False)
            uniform probability distribution function (if func=True)

    """

    np.random.seed()
    dist = uniform(loc=par['pmin'], scale=par['pmax']-par['pmin'])

    draw = dist.rvs()
              
    if func == False:
        return draw
    else:
        return dist

def beta_prior(par, func=False):
    """
    Beta prior.
  
    input: par ->  dictionary of parameter values
                   keywords: alpha, belta
           func -> boolean (optional)
                   if True returns the pdf random variable. 
                   Default is False.

    output: scalar (if func=False)
            beta probability distribution function (if func=True)

    """
    np.random.seed()
    rv = beta(par['alpha'], par['beta'])
    draw = beta.rvs(par['alpha'], par['beta']) 

    if func == False:
        return draw
    else:
        return rv
 

def main():
  print(__doc__)

if __name__=='__main__':
  main()
