#!/usr/bin/env python3

"""
Auxiliary script to test distance before using it in the ABC sampler. 

Usage:  test_ABC_distance.py -i <input_parameters_file> -f <user_function_file>
"""

import numpy as np
import matplotlib.pyplot as plt
import argparse
import imp
import os

from cosmoabc.distances import distance_quantiles, distance_GRBF 
from cosmoabc.priors import flat_prior, gaussian_prior, beta_prior
from cosmoabc.ABC_sampler import ABC
from cosmoabc.ABC_functions import SelectParamInnerLoop, DrawAllParams, SetDistanceFromSimulation, read_input 
        
import sys

def main(args):

    params = read_input( args.input )
   
    if args.functions != None:
        m1 = imp.load_source( args.functions[:-3], args.functions )

        if isinstance(params['distance_func'], list):
            params['distance_func'] = getattr(m1, params['distance_func'][0])

        if isinstance(params['simulation_func'], list):
            params['simulation_func'] = getattr(m1, params['simulation_func'][0])

        for par in params['param_to_fit']:
            if isinstance(params['prior'][par]['func'], str):
                params['prior'][par]['func'] = getattr(m1, params['prior_func'][params['param_to_fit'].index(par)])

    if not args.output:
        output_file = raw_input('Enter root for output files (no extension):  ')
    else:
        output_file = args.output

    if 'dataset1' not in params:
        params['dataset1'] = params['simulation_func'](params['simulation_input'])

    if 'cov' in params['simulation_input']:
        try:
            fname_cov = args.covariance
            params['simulation_input']['cov'] = np.loadtxt(fname_cov)
            params['cov'] = np.loadtxt(fname_cov)
        except IOError:
            print('Provide name of file containing covariance matrix!')

    #test distance between identical cataloges
    distance_equal = np.atleast_1d(params['distance_func'](params['dataset1'], params))

    if distance_equal.any() != 0.0:
        raise ValueError('Distance for 2 identical cataloges = ' + str(distance_equal))
    else:
        print('Distance between identical cataloges = ' + str(distance_equal))

    #test a single distance calculation
    new_sim_input = DrawAllParams(params['prior'])

    for j in range(params['npar']):
        params['simulation_input'][params['param_to_fit'][j]] = new_sim_input[j]
 
    data_simul = params['simulation_func'](params['simulation_input'])   

    try:
        distance_single = params['distance_func'](data_simul, params)
        print('New parameter value = ' + str(new_sim_input))
        print('Distance between observed and simulated data = ' + str(distance_single))
    except ValueError:
        print('Error in calculating single distance with parameters:')
        for item in params['param_to_fit']:
            print(item + '=' + str(params['simulation_input'][item]))

    if str(distance_single) is 'nan':
        print('NaN found!')
    
    #generate grid for distance behaviour inspection
    ngrid = int(raw_input('Enter number of draws in parameter grid: '))    
    param_grid = [DrawAllParams(params['prior']) for j1 in range(ngrid)]

    #open output data file
    op = open(output_file + '.dat', 'w')
    for par in params['param_to_fit']:
        op.write(par + '    ')
    for i1 in range(len(distance_equal)):
        op.write('dist' + str(i1 + 1) + '    ')
    op.write('\n')

    grid = []
    for pars in param_grid:

        if params['screen']:         
            print('Particle index: ' + str(len(grid)+1))
 
        grid_element = list(pars)
 
        for j1 in range(params['npar']):
            params['simulation_input'][params['param_to_fit'][j1]] = pars[j1]
         
        data_simul_grid = params['simulation_func'](params['simulation_input']) 
        distance_grid = np.atleast_1d(params['distance_func'](data_simul_grid, params))
             
        if sum(distance_grid) < 9**9:
            for item in grid_element:
                op.write(str(item) + '   ')
 
            for elem in distance_grid:
                grid_element.append(elem)
                op.write(str(elem) + '    ')
            op.write('\n')
  
            grid.append(grid_element)                   

    op.close()

    grid = np.array(grid)

    #plot distance behaviour
    plt.figure()
    plt.subplots_adjust(top = 0.95, right=0.985, left=0.075, bottom=0.075, 
                        wspace=0.25, hspace=0.25)
    
    n=0
    for par_indx in range(params['npar']):
        p1 = params['param_to_fit'][par_indx]
        for dist_indx in range(len(distance_single)):

            n = n + 1
            ylim = np.std( grid[:,params['npar'] + dist_indx])
            plt.subplot(params['npar'], len(distance_single), n)
            plt.scatter(grid[:,par_indx], grid[:,params['npar'] + dist_indx])
            plt.xlabel(params['param_to_fit'][par_indx], fontsize=14)
            plt.ylabel('distance' + str(dist_indx + 1), fontsize=14)
            #plt.xticks(range(int(params['prior'][p1]['min']) - 1, 
            #           int(params['prior'][p1]['max']) + 1, int(params['prior'][p1]['max'] - params['prior'][p1]['min'] +2)/5), fontsize=8)
            plt.yticks(fontsize = 8)
            plt.ylim(-0.25*ylim, ylim)
   
    plt.savefig(output_file + '.pdf')
    plt.close()
  
    print('\n Figure containing distance results is stored in ' + output_file + '.pdf')

if __name__=='__main__':
  
    #get user input file name
    parser = argparse.ArgumentParser(description='Approximate Bayesian Computation code.')
    parser.add_argument('-i', '--input', dest='input', help='User input file name.',required=True)
    parser.add_argument('-f', '--functions',  dest='functions', help='File name for user defined functions.', required=False)
    parser.add_argument('-o', '--output', dest='output', help='Output plot file name.', required=False)
    parser.add_argument('-c', '--covariance', dest='covariance', help='File name for covariance matrix.', required=False)
    args = parser.parse_args()
   
    main( args )

