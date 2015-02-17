#!/usr/bin/env python

"""
Auxiliary script to test distance before using it in the ABC sampler. 

Usage:  test_ABC_distance.py -i <input_parameters_file> -f <user_function_file>

        ***WARNING***
        at this point the distance function must be named 'distance'
                      the simulation function must be named 'simulation'
                      and if a personalized prior function is defined 
                      it must be named 'prior'.                   

"""

import numpy as np
import matplotlib.pyplot as plt
import argparse
import imp
import os

from CosmoABC.distances import distance_quantiles, summ_quantiles, distance_grbf, SumGRBF 
from CosmoABC.priors import flat_prior, gaussian_prior, beta_prior
from CosmoABC.ABC_sampler import ABC
from CosmoABC.ABC_functions import SelectParamInnerLoop, DrawAllParams, SetDistanceFromSimulation, read_input 
from CosmoABC.sim_NumCosmo import NCountSimul, ChooseParamsInput, numcosmo_simulation

        

def main(args):

    params = read_input( args.input )
   
    if args.functions != None:
        m1 = imp.load_source( args.functions[:-3], args.functions )

        if 'distance_func' not in params:
            params['distance_func'] = m1.distance

        if 'simulation_func' not in params:
            params['simulation_func'] = m1.simulation

        for i1 in xrange(params['npar']):
            if isinstance(params['prior_func'][i1], str):
                params['prior_func'][i1] = m1.prior

    if not args.output:
        output_file = raw_input('Enter output file name:  ')
    else:
        output_file = args.output

    if 'dataset1' not in params:
        params['dataset1'] = params['simulation_func'](params['simulation_input'])

        if params['distance_func'] == distance_grbf:
            params['extra'] =  SumGRBF(params['dataset1'], params['dataset1'], params)

        elif params['distance_func'] == distance_quantiles:
            params['dist_dim'] = len(params['dataset1'][0]) + 1
            params = summ_quantiles(params['dataset1'], params)

    #test distance between identical cataloges
    distance_equal = np.atleast_1d(params['distance_func'](params['dataset1'], params))

    if distance_equal.any() != 0.0:
        raise ValueError('Distance for 2 identical cataloges = ' + str(distance_equal))
    else:
        print 'Distance between identical cataloges = ' + str(distance_equal)

    #test a single distance calculation
    new_sim_input = DrawAllParams(params)

    for j in xrange(params['npar']):
        params['simulation_input'][params['param_to_fit'][j]] = new_sim_input[j]
 
    data_simul = params['simulation_func'](params['simulation_input'])   

    try:
        distance_single = params['distance_func'](data_simul, params)
        print 'New parameter value = ' + str(new_sim_input)
        print 'Distance between observed and simulated data = ' + str(distance_single) 
    except ValueError:
        print 'Error in calculating single distance with parameters:'
        for item in params['param_to_fit']:
            print item + '=' + str(params['simulation_input'][item])
    
    #generate grid for distance behaviour inspection
    ngrid = int(raw_input('Enter number of draws in parameter grid: '))    
    param_grid = [DrawAllParams(params) for j1 in xrange(ngrid)]

    grid = []
    for pars in param_grid:
        
        print 'Particle index: ' + str(len(grid)+1)
 
        grid_element = list(pars)
        for j1 in xrange(params['npar']):
            params['simulation_input'][params['param_to_fit'][j1]] = pars[j1]
         
        data_simul_grid = params['simulation_func'](params['simulation_input']) 
        distance_grid = np.atleast_1d(params['distance_func'](data_simul_grid, params))
             
        if sum(distance_grid) < 9**9:
            for elem in distance_grid:
                grid_element.append(elem)

            grid.append(grid_element)                   

    grid = np.array(grid)

    #plot distance behaviour
    plt.figure()
    plt.suptitle(params['distance_func'].__name__)
    plt.subplots_adjust(top = 0.95, right=0.985, left=0.075, bottom=0.075, wspace=0.25, hspace=0.25)
    
    n=0
    for par_indx in xrange(params['npar']):
        for dist_indx in xrange(len(distance_single)):

            n = n + 1
            ylim = np.std( grid[:,params['npar'] + dist_indx])
            plt.subplot(params['npar'], len(distance_single), n)
            plt.scatter(grid[:,par_indx], grid[:,params['npar'] + dist_indx])
            plt.xlabel(params['param_to_fit'][par_indx])
            plt.ylabel('distance' + str(dist_indx + 1), fontsize = 10)
            plt.xticks(np.arange(params['param_lim'][par_indx][0], params['param_lim'][par_indx][1], 
                       (params['param_lim'][par_indx][1]-params['param_lim'][par_indx][0])/5), 
                       fontsize=8)
            plt.yticks(fontsize = 8)
            plt.ylim(-0.25*ylim, ylim)
   
    plt.savefig(output_file)
    plt.close()
  
    print '\n Figure containing distance results is stored in ' + output_file

if __name__=='__main__':
  
    #get user input file name
    parser = argparse.ArgumentParser(description='Approximate Bayesian Computation code.')
    parser.add_argument('-i', '--input', dest='input', help='User input file name.',required=True)
    parser.add_argument('-f', '--functions',  dest='functions', help='File name for user defined functions.', required=False)
    parser.add_argument('-o', '--output', dest='output', help='Output plot file name.', required=False)
    args = parser.parse_args()
   
    main( args )

