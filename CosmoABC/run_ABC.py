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



from distances import *
from priors import *
import argparse

def read_input( filename ):
    """
    Read user input from file and construct initial dictionary parameter. 

    input:    filename (string) -> user input file parameter 

    output:   dictionary with formated user choices
    """

    #read user input data
    op1 = open( filename, 'r')
    lin1 = readlines()
    op1.close()

    data1 = [ elem.split() for elem in lin1 ]     

    #store options in params dictionary
    params_ini = dict( [ [ data1[0], data1[2:] ]  for line in data1 ] )
    

    #read observed data
    op2 = open( params_ini['path_to_obs'], 'r' )
    lin2 = op2.readlines()
    op2.close()

    data2 = [ elem.split() for elem in lin2[1:] ]


    params = {}
    params['dataset1'] = numpy.array([ [ float( item ) for item in line ] for line in data2 ])  
    params['param_to_fit'] = [ params_ini['param_to_fit'][ i ] for i in xrange( params_ini['param_to_fit'].index('#') ) ]
 
    params['npar'] = len( params['param_to_fit'] )
    params['prior_par'] = [ [ float( params_ini[ params[ 'param_to_fit' ] + '_prior_par' ][ j ] ) for j in xrange(2) ] for i in xrange( npar ) ]
    params['param_lim'] = [ [ float( params_ini[ params[ 'param_to_fit' ] + '_lim' ][ j ] ) for j in xrange(2) ] for i in xrange( npar ) ]
    params['M'] = int( params_ini['M'][0] )
    params['epsilon1'] = float( params_ini[ 'epsilon1' ][0] )
    params['qthreshold'] = float( params_ini['qthreshold'][0])
    params['delta'] = float( params_ini['delta'][0] )
    params['s'] =  float( params_ini['s'][0] )
    params['file_root'] = params_ini['file_root'][0]  

    #functions
    ###### Update this if you include any new functions!!!!!  ##############
    dispatcher = {'NumCosmo_simulation': NumCosmo_simulation, 'flat_prior': flat_prior, 'gaussian_prior': gaussian_prior, 'distance_GRBF':distance_GRBF}
    
    params['simulation_func'] = dispatcher[ params_ini['simulation_func'][0] ]
    params['prior_func'] = dispatcher[ params_ini[ 'prior_func' ][0] ]
    params['distance_func'] = dispatcher[ params_ini[ 'distance_func' ][0] ]

    #fiducial extra parameters
    sim_par = {}  
    for item in params_ini.keys():
        if item not in params.keys():
            sim_par[ item ] = float( params_ini[ item ][0] ) 

    params['simulation_params'] = sim_par

    return params
    

def main( args ):

    user_input = read_input( args.input )

    #initiate ABC construct
    sampler_ABC = ABC( dataset1=user_input['dataset1'], params=user_input, simulation_func=user_input['simulation_func'], prior_func=user_input['prior_func'], distance_func=user_input['distance_func']) 

    #build first particle system
    sys1 = sampler_ABC.BuildFirstPSystem( filename=user_input['file_root'] + '0.dat' )

    #update particle system until convergence
    sampler_ABC.fullABC(  user_input['file_root'] )
        
         

if __name__=='__main__':
  
    #get user input file name
    parser = argparse.ArgumentParser(description='Approximate Bayesian Computation code.')
    parser.add_argument('-i','--input', help='Input file name',required=True)
    args = parser.parse_args()
   
    main( args )


