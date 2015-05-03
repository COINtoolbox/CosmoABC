#!/usr/bin/env python

"""
Approximate Bayesian Computation package. 


Usage:
 - 	python run_ABC.py  -i <user_input_file>

Files description:
 - 	priors.py:			Functions for initial prior PDF.
 - 	distances.py:			Functions for determining distances between catalogs.
 -  	ABC_sampler.py:			The ABC class. 
 - 	run_ABC.py:			Warp of the functionalities of the ABC class.
 
Tests:
 -	toy_model.py:			Simple Gaussian example.


Libraries used:

 -	scipy.stats:			For probability distribution functions.
 -	numpy:				For basic math.
 -	time:				For tracking.
 -	statsmodels.stats.weightstats:	For weighted covariances matrixes.
 -	sys:				For maintainence.
 
Optional external dependences:

 -	NumCosmo:			For cosmological simulations.	
"""

import argparse
import imp

from CosmoABC.distances import distance_quantiles, distance_GRBF
from CosmoABC.priors import flat_prior, gaussian_prior, beta_prior
from CosmoABC.ABC_sampler import ABC
from CosmoABC.plots import plot_1p, plot_2p, plot_3p, plot_4p  
from CosmoABC.ABC_functions import read_input, get_cores, DrawAllParams, SetDistanceFromSimulation, SelectParamInnerLoop

def main( args ):

    user_input = read_input(args.input)
    m1 = imp.load_source(args.functions[:-3], args.functions)

    user_input['simulation_func'] = m1.simulation


    if  isinstance(user_input['distance_func'], list):
        user_input['distance_func'] = m1.distance

    for l1 in range(user_input['npar']):
        if isinstance(user_input['prior_func'][ l1 ], str):            
            user_input['prior_func'][ l1 ] = getattr(m1, user_input['prior_func'][ l1 ])
            
    #check if observed data exist, simulate in case negative
    if user_input['path_to_obs'] == 'None':
        user_input['dataset1'] = user_input['simulation_func'](user_input['simulation_input'])

     
    if 'dist_dim' not in user_input.keys():    
        user_input['dist_dim'] = len(user_input['distance_func'](user_input['dataset1'], user_input))

    #initiate ABC construct
    sampler_ABC = ABC(params=user_input) 

    #build first particle system
    sys1 = sampler_ABC.BuildFirstPSystem()

    #update particle system until convergence
    sampler_ABC.fullABC()

    #plot results
    if len(user_input['param_to_fit'] ) == 1 :
        plot_1p(sampler_ABC.T, 'results.pdf', user_input)

    elif len(user_input['param_to_fit'] ) == 2 :
        plot_2p(sampler_ABC.T, 'results.pdf', user_input)      

    elif len(user_input['param_to_fit'] ) == 3 :
        plot_3p(sampler_ABC.T, 'results.pdf', user_input) 

    elif len(user_input['param_to_fit']) == 4:
        plot_4p(sampler_ABC.T, 'results.pdf', user_input)

    else:
        raise ValueError('Only 1, 2, 3, and 4 dimensional plots are implemented so far!')

if __name__=='__main__':
  
    #get user input file name
    parser = argparse.ArgumentParser(description='Approximate Bayesian Computation code.')
    parser.add_argument('-i','--input', dest='input', help='User input file name.',required=True)
    parser.add_argument('-f','--functions',  dest='functions', help='File name for user defined functions.', required=True)
    args = parser.parse_args()
   
    main( args )


