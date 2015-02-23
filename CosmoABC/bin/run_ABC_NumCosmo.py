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


__author__ = "E. E. O. Ishida"
__maintainer__ = "E. E. O. Ishida"
__copyright__ = "Copyright 2015"
__email__ = "emilleishida@gmail.com"
__status__ = "Prototype"
__license__ = "GPL"



import argparse
import numpy
import imp
import os

from CosmoABC.distances import distance_quantiles, summ_quantiles, distance_grbf, SumGRBF 
from CosmoABC.priors import flat_prior, gaussian_prior, beta_prior
from CosmoABC.ABC_sampler import ABC
from CosmoABC.plots import plot_1D, plot_2D, plot_3D, plot_4D
from CosmoABC.ABC_functions import SelectParamInnerLoop, DrawAllParams, SetDistanceFromSimulation, read_input 


try: 
    from gi.repository import NumCosmo as Nc
    from CosmoABC.sim_NumCosmo  import NCountSimul, ChooseParamsInput, numcosmo_simulation
except ImportError:
    raise ImportError('You must have NumCosmo running to use the sim_NumCosmo simulation!' + 
                      '\n Please check your NumCosmo instalation.')

def main( args ):

    user_input = read_input( args.input )

    #initiate NumCosmo object necessary for simulation
    Cosmo=ChooseParamsInput()

 
    Cosmo.params = user_input['simulation_input']
    Cosmo.params["OL"]  = 1.- Cosmo.params['Om']-Cosmo.params['Ob']

    #assign input for simulation
    user_input['simulation_params'] = Cosmo.params

    if args.functions != None:
        m1 = imp.load_source( args.functions[:-3], args.functions )

        if 'distance_func' not in user_input.keys():
            user_input['distance_func'] = m1.distance
    
    for l1 in range( user_input['npar'] ):
        if isinstance( user_input['prior_func'][ l1 ], str):            
            user_input['prior_func'][ l1 ] = getattr( m1, user_input['prior_func'][ l1 ] )


    #initiate ABC construct
    sampler_ABC = ABC(params=user_input) 

    
    #build first particle system
    sys1 = sampler_ABC.BuildFirstPSystem()

    #update particle system until convergence
    sampler_ABC.fullABC()
     
         

if __name__=='__main__':
  
    #get user input file name
    parser = argparse.ArgumentParser(description='Approximate Bayesian Computation code.')
    parser.add_argument('-i','--input', help='Input file name',required=True)
    parser.add_argument('-f', '--functions', help='User defined functions', required=False, default=None)
    args = parser.parse_args()
   
    main( args )


