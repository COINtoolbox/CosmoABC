#!/usr/bin/env python

"""
This file contains a simple example of simulation and distance function. 
Use this as an example to implement your own functions. 

.. warning::
    For consistency, the 
    -  **simulation** function **must be called *'simulation'* and take a dictionary as input

    - the **distance function **must be called *'distance'* and take 2 data sets 
       and an extra input parameter as input**

    - the **prior** function **must take 2 arrays and 1 bool as input**:
          first input: list of prior parameters (one vector for each parameter being fitted)
          second input: list of extreme limits, defining the interval where the prior is 
                        expected to be valid (this was included to avoid unecessary and 
                                              too long calculations)
          third input: a boolean variable called *func*.
                       it regulates wether the user wants as output a sampling of the PDF or the 
                       PDF function itself. 
                       *This is important for the calculation of weights!*

"""



import numpy
from CosmoABC.priors import gaussian_prior
from scipy.stats import norm


def simulation( v ):
    """
    Simple simulation function.
    Generates one realization of a Gaussian probability distribution function with mean v['mean'] and standard deviation v['sigma'] containing v['n'] data points.

    input:   v -> dictionary of parameter values
                  keys: mean  -> mean of Gaussian PDF
                        std   -> standard deviation of Gaussian PDF
                        n     -> number of points in the data catalog    

    output:  array of simulated catalog
    """

    l1 = numpy.random.normal( loc=v['mean'], scale=v['std'], size=v['n'] )
    
    return numpy.atleast_2d( l1 ).T 


def distance( dataset1, dataset2, s ):
    """
    Simple distance function from 2 independent data catalogs generated with function my_sim.simulatio.

    input:    dataset1, dataset2  -> catalogs to be compared 
              s                   -> extra parameter.
                                     *included for consistency with other functions but is not used in this specific example.*  

    output:   distance between dataset1 and dataset2 (float)
              
    """
    
    #difference between estimated mean 
    t1 = abs( numpy.mean( dataset1.flatten() ) - numpy.mean( dataset2.flatten() )) 

    #difference between estimated standard deviation
    t2 = abs( numpy.std(  dataset1.flatten() ) - numpy.std(  dataset2.flatten() ))

    #difference between number of data points
    t3 = abs( len( dataset1 ) - len( dataset2 ) )

    return t1 + t2 + t3


def prior( prior_param, prior_lim, func=False):
    """
    Example of a simple prior function.
    Draw a parameter value from a flat prior.
  
    input: 	par          -> vector of parameter required by the corresponding distribution family.
                                format: [ lower_bound, upper_bound]

                par_lim      -> physical reasonable limits for the cosmological parameters
                                2-dimensional vector with [min_value, max_value] for each parameter  

    output:     scalar	     -> draw number or 1
    """

    draw = numpy.random.uniform( low=prior_param[0], high=prior_param[1] )
             
    if func == False:
        return draw
    else:
        return norm( loc=prior_lim[0], scale=prior_lim[1]-prior_lim[0])


def main():
  print(__doc__)



if __name__=='__main__':
  main() 
