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


__author__ = "E. E. O. Ishida, S. D. P. Vitenti, M. Penna-Lima"
__maintainer__ = "E. E. O. Ishida"
__copyright__ = "Copyright 2015"
__version__ = "0.1"
__email__ = "emilleishida@gmail.com"
__status__ = "Prototype"
__license__ = "GPL"



import argparse
import numpy
from CosmoABC.distances import * 
from CosmoABC.priors import *
from CosmoABC.ABC_sampler import *
from CosmoABC.plots import *
import imp



def read_input( filename ):
    """
    Read user input from file and construct initial dictionary parameter. 

    input:    filename (string) -> user input file parameter 

    output:   dictionary with formated user choices
    """

    
 
    #read user input data
    op1 = open( filename, 'r')
    lin1 = op1.readlines()
    op1.close()

    data1 = [ elem.split() for elem in lin1 ]     

    #store options in params dictionary
    params_ini = dict( [ ( line[0], line[2:] )  for line in data1 if len( line ) > 1 ] )
    

    #read observed data
    op2 = open( params_ini['path_to_obs'][0], 'r' )
    lin2 = op2.readlines()
    op2.close()

    data2 = [ elem.split() for elem in lin2[1:] ]


    params = {}
    params['path_to_obs'] = params_ini['path_to_obs'][0] 
    params['dataset1'] = numpy.array([ [ float( item ) for item in line ] for line in data2 ])  
    params['param_to_fit'] = [ params_ini['param_to_fit'][ i ] for i in xrange( params_ini['param_to_fit'].index('#') ) ]
 
    params['npar'] = len( params['param_to_fit'] )
    params['prior_par'] = [ [ float( params_ini[ params[ 'param_to_fit' ][ i ] + '_prior_par' ][ j ] ) for j in xrange(2) ] for i in xrange( params['npar'] ) ]
    params['param_lim'] = [ [ float( params_ini[ params[ 'param_to_fit' ][ i ] + '_lim' ][ j ] ) for j in xrange(2) ] for i in xrange( params['npar'] ) ]
    params['M'] = int( params_ini['M'][0] )
    params['epsilon1'] = float( params_ini[ 'epsilon1' ][0] )
    params['qthreshold'] = float( params_ini['qthreshold'][0])
    params['delta'] = float( params_ini['delta'][0] )
    params['s'] =  float( params_ini['s'][0] )
    params['file_root'] = params_ini['file_root'][0]  


    #fiducial extra parameters
    sim_par = {}  
    for item in params_ini.keys():
        if item not in params.keys():
            try:
                float( params_ini[ item ][0] )
                sim_par[ item ] = float( params_ini[ item ][0] )
            except ValueError:
                sim_par[ item ] = params_ini[ item ][0] 

    params['simulation_params'] = sim_par

    return params
    

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
    sampler_ABC = ABC( dataset1=user_input['dataset1'], params=user_input, simulation_func=user_input['simulation_func'], prior_func=user_input['prior_func'], distance_func=user_input['distance_func']) 

    #define finished particle system index
    sampler_ABC.T = int( args.PS )
    
    #continue from previous run
    sampler_ABC.ContinueStoppedRun( sampler_ABC.T , user_input['file_root'] )

    #plot results
    user_input['param_lim']=[ [ item/2.0 for item in line] for line in user_input['param_lim']]
    if len( user_input['param_to_fit'] ) == 1 :
        plot_1D( sampler_ABC.T, 'results.pdf', user_input)

    elif len( user_input['param_to_fit'] ) == 2 :
        plot_2D( sampler_ABC.T, 'results.pdf', user_input )   

    elif len( user_input['param_to_fit'] ) == 3 :
        plot_3D( sampler_ABC.T, 'results.pdf', user_input )       

    else:
        raise ValueError('Only 1, 2 and 3 dimensional plots are implemented so far!')

if __name__=='__main__':
  
    #get user input file name
    parser = argparse.ArgumentParser(description='Approximate Bayesian Computation code.')
    parser.add_argument('-i','--input', dest='input', help='User input file name.',required=True)
    parser.add_argument('-f','--functions',  dest='functions', help='File name for user defined functions.', required=True)
    parser.add_argument('-p','--PS', dest='PS', help='Last particle system index completed.', required=True)
    args = parser.parse_args()
   
    main( args )


