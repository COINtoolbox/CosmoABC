"""
Functions for initial prior distributions. 

***
To insert a new family of prior distribution, insert your function in this file.
All functions must depend on 2 variables:
One defining the parameters characterizing the distribution (par)
and another determining the acceptable extreme bounds ( par_lim ). 
There must also be an optional variable regulating the output (rather a draw or a pdf)
***

"""

import numpy as np

from scipy.stats import norm
from scipy.stats import uniform
from scipy.stats import beta

def gaussian_prior( par, par_lim, func=False ):
    """
    Draw a parameter value from a Gaussian prior.
  
    input: 	par          -> vector of parameter required by the corresponding distribution family.
                                format: [ mean, standard_devitaion]

                par_lim      -> physical reasonable limits for the cosmological parameters
                                2-dimensional vector with [min_value, max_value] for each parameter

		func (optional)	->   return the pdf random variable (boolean). Default is False.
 
                

    output: scalar	-> one random draw
            or
            pdf		-> probability distribution function

    """
    np.random.seed()
   
    #check dimension of feature vector  defining prior distribution 
    if len( par ) == 2:

        flag = False

        #draw a parameter value physically meaningfull
        while flag == False:
            draw = np.random.normal( loc=par[0], scale=par[1] ) 
                
            if par_lim[0] <= draw and draw < par_lim[1]:
                flag = True

    else:
        raise ValueError("Gaussian distribution requires 2-dimensional parameter vector: [mean, standard_deviation].")


    if func == False:
        return draw
    else:
        return norm( loc=par[0],  scale=par[1])


def flat_prior( par, par_lim, func=False ):
    """
    Draw a parameter value from a flat prior.
  
    input: 	par          -> vector of parameter required by the corresponding distribution family.
                                format: [ lower_bound, upper_bound]

                par_lim      -> physical reasonable limits for the cosmological parameters
                                2-dimensional vector with [min_value, max_value] for each parameter  

                func (optional)	->   return the pdf random variable (boolean). Default is False.

    output:     scalar	     -> draw number or 1
    """
    np.random.seed()

    #check dimensional of feature vector defining distribution
    #if distribution is flat there is no need to check the bounderies
    if len( par ) == 2 and par[0] < par[1]:
        draw = np.random.uniform( low=par[0], high=par[1] )
    else:
        raise ValueError("Flat distribution requires 2-dimensional parameter vector: [lower_bound, upper_bound].")
               
    if func == False:
        return draw
    else:
        return uniform( loc=par_lim[0], scale=par_lim[1]-par_lim[0])

def beta_prior( par, par_lim, func=False):
    """
    Draw a parameter value from a beta prior.
  
    input: 	par          -> vector of parameter required by the corresponding distribution family.
                                format: [ lower_bound, upper_bound]

                par_lim      -> physical reasonable limits for the cosmological parameters
                                2-dimensional vector with [min_value, max_value] for each parameter  

                func (optional)	->   return the pdf random variable (boolean). Default is False.

    output:     scalar	     -> draw number or 1
    """
    np.random.seed()

    #check dimension of feature vector  defining prior distribution 
    if len( par ) == 2:

        rv = beta(par[0], par[1])
        draw = beta.rvs( par[0], par[1]) 

    else:
        raise ValueError("Beta distribution requires 2-dimensional parameter vector: [mean, standard_deviation].")


    if func == False:
        return par[2] * draw
    else:
        return (1.0/abs(par[2])) * rv
  



def main():
  print(__doc__)

if __name__=='__main__':
  main()
