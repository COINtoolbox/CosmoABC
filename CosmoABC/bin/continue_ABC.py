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


__author__ = "E. E. O. Ishida, S. D. P. Vitenti, M. Penna-Lima,  R. S. de Souza, J. Cisewski, E. Cameron, V. C. Busti"
__maintainer__ = "E. E. O. Ishida"
__copyright__ = "Copyright 2015"
__version__ = "0.1"
__email__ = "emilleishida@gmail.com"
__status__ = "Prototype"
__license__ = "GPL"

import argparse
import imp

from CosmoABC.distances import distance_quantiles, summ_quantiles, distance_grbf, SumGRBF 
from CosmoABC.priors import flat_prior, gaussian_prior, beta_prior
from CosmoABC.ABC_sampler import ABC
from CosmoABC.ABC_functions import SelectParamInnerLoop, DrawAllParams, SetDistanceFromSimulation, read_input 

def main( args ):

    user_input = read_input( args.input )

    m1 = imp.load_source( args.functions[:-3], args.functions )

    user_input['simulation_func'] = m1.simulation

    if 'distance_func' not in user_input.keys():
        user_input['distance_func'] = m1.distance
    
    for l1 in range( user_input['npar'] ):
        if isinstance( user_input['prior_func'][ l1 ], str):            
            user_input['prior_func'][ l1 ] = getattr( m1, user_input['prior_func'][ l1 ] )
            
    #initiate ABC construct
    sampler_ABC = ABC( params=user_input ) 

    #define finished particle system index
    sampler_ABC.T = int( args.PS )
    
    #continue from previous run
    sampler_ABC.ContinueStoppedRun( sampler_ABC.T , user_input['file_root'] )



if __name__=='__main__':

    parser = argparse.ArgumentParser(description='Approximate Bayesian Computation code.')
    parser.add_argument('-i','--input', dest='input', help='User input file name.',required=True)
    parser.add_argument('-f','--functions',  dest='functions', help='File name for user defined functions.', required=True)
    parser.add_argument('-p','--PS', dest='PS', help='Last particle system index completed.', required=True)
    args = parser.parse_args()
   
    main( args )


