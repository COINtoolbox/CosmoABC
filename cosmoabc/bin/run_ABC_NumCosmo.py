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
import numpy
import imp
import os

from cosmoabc.distances import distance_quantiles, distance_GRBF 
from cosmoabc.priors import flat_prior, gaussian_prior, beta_prior
from cosmoabc.ABC_sampler import ABC
from cosmoabc.ABC_functions import SelectParamInnerLoop, DrawAllParams, SetDistanceFromSimulation, read_input 

try: 
    from gi.repository import NumCosmo as Nc
    from cosmoabc.sim_NumCosmo_cluster  import NCountSimul, ChooseParamsInput, numcosmo_sim_cluster
except ImportError:
    raise ImportError('You must have NumCosmo running to use the sim_NumCosmo simulation!' + 
                      '\n Please check your NumCosmo instalation.')

def main( args ):

    user_input = read_input(args.input)

    if args.functions != None:
        m1 = imp.load_source( args.functions[:-3], args.functions )

        if isinstance(user_input['distance_func'][0], str):
            user_input['distance_func'] = getattr(m1, user_input['distance_func'][0])
    
    for l1 in range(user_input['npar']):
        par = user_input['param_to_fit'][l1]
        if isinstance(user_input['prior'][par]['func'], str):            
            user_input['prior'][par]['func'] = getattr(m1, user_input['prior_func'][l1])

    #initiate ABC construct
    sampler_ABC = ABC(params=user_input) 

    #build first particle system
    sys1 = sampler_ABC.BuildFirstPSystem()

    #update particle system until convergence
    sampler_ABC.fullABC(nruns=int(user_input['nruns'][0]))
     
         

if __name__=='__main__':
  
    #get user input file name
    parser = argparse.ArgumentParser(description='Approximate Bayesian Computation code.')
    parser.add_argument('-i','--input', help='Input file name',required=True)
    parser.add_argument('-f', '--functions', help='User defined functions', required=False, default=None)
    args = parser.parse_args()
   
    main( args )


