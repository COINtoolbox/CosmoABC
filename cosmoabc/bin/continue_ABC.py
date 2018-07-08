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

from cosmoabc.distances import distance_quantiles, distance_GRBF
from cosmoabc.priors import flat_prior, gaussian_prior, beta_prior
from cosmoabc.ABC_sampler import ABC
from cosmoabc.ABC_functions import SelectParamInnerLoop, DrawAllParams, SetDistanceFromSimulation, read_input 

def main(args):

    user_input = read_input(args.input)

    if user_input['path_to_obs'] == 'None':
        raise IOError('It is not possible to continue a process without determining a static data set.')    

    m1 = imp.load_source(args.functions[:-3], args.functions)

    user_input['simulation_func'] = getattr(m1, user_input['simulation_func'][0])

    if isinstance(user_input['distance_func'][0], str):
        user_input['distance_func'] = getattr(m1, user_input['distance_func'][0])

    #check dimension of distance output
    dtemp = user_input['distance_func'](user_input['dataset1'], user_input)
    user_input['dist_dim'] = len(dtemp)
    
    for l1 in range(user_input['npar']):
        par = user_input['param_to_fit'][l1]
        if isinstance(user_input['prior'][par]['func'], str):            
            user_input['prior'][par]['func'] = getattr(m1, user_input['prior_func'][l1])
            
    #initiate ABC construct
    sampler_ABC = ABC(params=user_input) 

    #define finished particle system index
    sampler_ABC.T = int(args.PS)
    
    #continue from previous run
    sampler_ABC.ContinueStoppedRun(sampler_ABC.T, nruns=int(user_input['nruns'][0]))



if __name__=='__main__':

    parser = argparse.ArgumentParser(description='Approximate Bayesian Computation code.')
    parser.add_argument('-i','--input', dest='input', help='User input file name.',required=True)
    parser.add_argument('-f','--functions',  dest='functions', help='File name for user defined functions.', required=True)
    parser.add_argument('-p','--PS', dest='PS', help='Last particle system index completed.', required=True)
    args = parser.parse_args()
   
    main( args )


