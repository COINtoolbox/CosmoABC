#!/usr/bin/env python

"""
Approximate Bayesian Computation package. 


Usage:
 - 	python run_ABC.py  -i <user_input_file> -f <user_function_files>

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
from cosmoabc.plots import plot_1p, plot_2p, plot_3p 
from cosmoabc.ABC_functions import read_input, get_cores, DrawAllParams, SetDistanceFromSimulation, SelectParamInnerLoop

def main( args ):

    user_input = read_input(args.input)
    m1 = imp.load_source(args.functions[:-3], args.functions)

    user_input['simulation_func'] = getattr(m1, user_input['simulation_func'][0])


    if  isinstance(user_input['distance_func'], list):
        user_input['distance_func'] = getattr(m1, user_input['distance_func'][0])

    for element in user_input['param_to_fit']:
        if isinstance(user_input['prior'][element]['func'], str):            
            user_input['prior'][element]['func'] = getattr(m1, user_input['prior'][element]['func'])
            for pvar in user_input[element + '_prior_par_name'][:user_input[element + '_prior_par_name'].index('#')]:
                indx = user_input[element + '_prior_par_name'].index(pvar)
                user_input['prior'][element][pvar] = float(user_input[element + '_prior_par_val'][indx])
            
    #check if observed data exist, simulate in case negative
    if user_input['path_to_obs'] == 'None':
        user_input['dataset1'] = user_input['simulation_func'](user_input['simulation_input'])

        #save data to file 
        op1 = open('synthetic_data.dat', 'w')
        for line in user_input['dataset1']:
            for item in line:
                op1.write(str(item) + '    ')
            op1.write('\n')
        op1.close()

    dist_try = user_input['distance_func'](user_input['dataset1'], user_input)
    if 'dist_dim' not in user_input.keys():    
        user_input['dist_dim'] = len(dist_try)


    #initiate ABC construct
    sampler_ABC = ABC(params=user_input) 

    #build first particle system
    sys1 = sampler_ABC.BuildFirstPSystem()

    #update particle system until convergence
    
    sampler_ABC.fullABC(nruns=int(user_input['nruns'][0]))

    #plot results
    if len(user_input['param_to_fit'] ) == 1 :
        plot_1p(sampler_ABC.T, 'results.pdf', user_input)

    elif len(user_input['param_to_fit'] ) == 2 :
        plot_2p(sampler_ABC.T, 'results.pdf', user_input)      

    elif len(user_input['param_to_fit'] ) == 3 :
        plot_3p(sampler_ABC.T, 'results.pdf', user_input)

    else:
        raise ValueError('Only 1, 2, and 3 dimensional plots are implemented so far!')

if __name__=='__main__':
  
    #get user input file name
    parser = argparse.ArgumentParser(description='Approximate Bayesian Computation code.')
    parser.add_argument('-i','--input', dest='input', help='User input file name.',required=True)
    parser.add_argument('-f','--functions',  dest='functions', help='File name for user defined functions.', required=True)
    args = parser.parse_args()
   
    main( args )


