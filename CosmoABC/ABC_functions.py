#!/usr/bin/env python

"""
Functions to optimize parallelization.
"""
import numpy as np
import time
import sys
	
from inspect import isfunction

from scipy.stats import multivariate_normal
from scipy.interpolate import Rbf
from scipy.stats.mstats import mquantiles
from scipy.stats import norm

from statsmodels.stats.weightstats import DescrStatsW

from CosmoABC.distances import distance_quantiles, summ_quantiles, distance_grbf, SumGRBF 
from CosmoABC.priors import flat_prior, gaussian_prior, beta_prior


################################################

def get_cores():
    """
    Ask the user to input number of cores.
    """
    cores = raw_input('Please enter the number of cores: ')

    return int(cores) 

def read_input(filename):
    """
    Read user input from file and construct initial dictionary parameter. 

    input:    filename (string) -> user input file parameter 

    output:   dictionary with formated user choices
    """

    #read user input data
    op1 = open(filename, 'r')
    lin1 = op1.readlines()
    op1.close()

    data1 = [elem.split() for elem in lin1]     

    #store options in params dictionary
    params_ini = dict([(line[0], line[2:])  for line in data1 if len(line) > 1])
    
    params = {}
    params['path_to_obs'] = params_ini['path_to_obs'][0] 

    #check if ``observer'' data already exists
    if params['path_to_obs'] != 'None':
        #read observed data
        op2 = open(params_ini['path_to_obs'][0], 'r')
        lin2 = op2.readlines()
        op2.close()

        data2 = [elem.split() for elem in lin2[1:]]
        params['dataset1'] = np.array([[float(item) for item in line ] for line in data2])  


    params['param_to_fit'] = [params_ini['param_to_fit'][i] 
                             for i in xrange(params_ini['param_to_fit'].index('#'))]
    params['param_to_sim'] = [params_ini['param_to_sim'][i] 
                             for i in xrange(params_ini['param_to_sim'].index('#'))]
    params['npar'] = len(params['param_to_fit'])
    params['prior_par'] = [[float(params_ini[params['param_to_fit'][i] + '_prior_par'][j]) 
                          for j in xrange(2) ] for i in xrange(params['npar'])]
    params['param_lim'] = [[float(params_ini[params['param_to_fit'][i] + '_lim'][j]) 
                          for j in xrange(2) ] for i in xrange(params['npar'])]
    params['M'] = int(params_ini['M'][0])
    params['Mini'] = int(params_ini['Mini'][0])
    params['qthreshold'] = float(params_ini['qthreshold'][0])
    params['delta'] = float(params_ini['delta'][0])
    params['file_root'] = params_ini['file_root'][0]  
    params['screen'] = bool(params_ini['screen'][0])
    params['ncores'] = int(params_ini['ncores'][0])

    #functions
    ###### Update this if you include any new functions!!!!!  ##############
    dispatcher = {'flat_prior': flat_prior, 
                  'gaussian_prior': gaussian_prior, 'beta_prior':beta_prior, 
                  'distance_grbf':distance_grbf, 'distance_quantiles': distance_quantiles}
    
    if params_ini['simulation_func'][0] in dispatcher:
        params['simulation_func'] = dispatcher[params_ini['simulation_func'][0]]

    params['prior_func'] = [dispatcher[params_ini['prior_func'][k]] 
                           if params_ini['prior_func'][k] in dispatcher.keys() 
                           else params_ini['prior_func'][k] for k in xrange(params['npar'])]

    #fiducial extra parameters
    sim_par = {}  
    for item in params['param_to_sim']:
        try:
            par = float(params_ini[ item ][0])
            sim_par[ item ] = par
        except ValueError:
            sim_par[ item ] = params_ini[ item ][0] 

    if params_ini['simulation_func'] == 'numcosmo_simulation':

        try: 
            from gi.repository import NumCosmo as Nc
            from CosmoABC.sim_NumCosmo import NCountSimul, ChooseParamsInput, numcosmo_simulation
        except ImportError:
            raise ImportError( 'You must have NumCosmo running to use the sim_NumCosmo simulation!' +
                                '\n Please check your NumCosmo instalation.' )

        sim_par["OL"] = 1. - sim_par["Om"] - sim_par["Ob"]
        Cosmo=ChooseParamsInput()
        Cosmo.params = sim_par 
        params['simulation_input'] = Cosmo.params

    else:
        params['simulation_input'] = sim_par

    if params_ini['distance_func'][0] in dispatcher.keys():
        params['distance_func'] = dispatcher[params_ini['distance_func'][0]]

        if params['distance_func'] == distance_grbf:
            params['s'] = float(params_ini['s'][0])
            params['kernel_func'] = str(params_ini['kernel_func'][0])

            if 'dataset1' in params:
                params['extra'] =  SumGRBF(params['dataset1'], params['dataset1'], params)

        elif 'dataset1' in params and params['distance_func'] == distance_quantiles:
            params['dist_dim'] = len(params['dataset1'][0]) + 1
            params = summ_quantiles(params['dataset1'], params)

    return params

def SelectParamInnerLoop(var1):
    """
    Draw model parameters based on previous particle system and return those satisfying distance threshold.

    :param	previous_particle_system:  model parameters surviving previous distance threshold
                collumns -> [ param_to_fit[0], param_to_fit[1], .., param_to_fit[n], distance_from_dataset1 ] 	

    :param	W:  vector of weights 

    :param	previous_cov_matrix: covariance matrix based on previous particle system results
    :param	cosmological_parameters: dictionary of necessary cosmological parameters
		keys must include: [H0, Omegab, Omegam, Tgamma0, ns, sigma8, w]	

    :param	epsilon: list of distance threshold to be satisfied	

    :returns: vector -> [ surviving_model_parameters, distance, number_necessary_draws, computational_time (s), distance_threshold ]
    """

    dist = [10**10 for i in xrange(var1['params']['dist_dim'])]
    K = 0

    #time marker
    time_start = time.time()

    #determine distance threshold
    epsilon = [mquantiles(var1['previous_particle_system'][:, kk ], prob=var1['params']['qthreshold'])[0] 
                  for kk in xrange(len(var1['params']['param_to_fit']), 
                                   len(var1['params']['param_to_fit']) + var1['params']['dist_dim'])]

    flag = [False for ii in xrange(var1['params']['dist_dim'])]

    while False in flag:               
 
        #update counter
        K = K + 1 

        #draw model parameters to serve as mean 
        index_theta0 = np.random.choice(xrange(len(var1['W'])), p=var1['W'])
        theta0 = np.atleast_2d(var1['previous_particle_system'][index_theta0][:len( var1['params']['param_to_fit'])])

        #initialize boolean parmeter vector
        theta_t = [False for i in xrange(len(var1['params']['param_to_fit']))]
           
        while False in theta_t:
          
            #draw model parameter values for simulation
            mvn = multivariate_normal.rvs(theta0.flatten(), var1['previous_cov_matrix'])

            try:
                len(mvn)
                theta_t_try = list(mvn)

            except TypeError:
                theta_t_try = [mvn]

            theta_t = []
            for k1 in xrange(len(var1['params']['param_to_fit'])):

                if  theta_t_try[k1] >= var1['params']['param_lim'][k1][0] and theta_t_try[k1] < var1['params']['param_lim'][k1][1]:
                    theta_t.append(True)

                else:
                    theta_t.append(False)         


        #update parameter values in dictionary
        for i1 in range(len(var1['params']['param_to_fit'])):
            var1['params']['simulation_input'][var1['params']['param_to_fit'][i1]] = theta_t_try[i1]

        #generate simulation
        DataSimul = var1['params']['simulation_func'](var1['params']['simulation_input'])
      
        #calculate distance
        dist = var1['params']['distance_func'](DataSimul, var1['params'])
              
        #check if it satisfies distance thresholds 
        flag = [] 
        for ll in xrange(len( dist )):
            if dist[ll] > epsilon[ll]:
                flag.append(False)
            else:
                flag.append(True)

    #store results
    for d2 in dist:
        theta_t_try.append(d2)
    theta_t_try.append(K)    
    theta_t_try.append(time.time() - time_start)

    for d3 in epsilon:
        theta_t_try.append(d3)

    if var1['params']['screen']:
        print 'Number of draws = ' + str(K)

    return theta_t_try

def DrawAllParams(params):
    """
    Draw complete set of  parameters from prior.

    :returns: array of parameters sampled from the prior
    """

    pars = []
    for j in range(len(params['param_to_fit'])):
        p1 = params['prior_func'][j](params['prior_par'][j], params['param_lim'][j])   
        pars.append(p1)

    return np.array(pars)

def SetDistanceFromSimulation(var):
    """
    Draw cosmological parameter values from prior, generate simulation and calculate distance from a given comparison  catalog. 
    In the context of Ishida et al., 2015 that translates into calculating the distance from the real to the simulated cataloges.

    :returns: scalar (distance between dataset1 and dataset2)
    
    """

    time1 = time.time()

    #draw parameters from prior
    ParTry = DrawAllParams(var)
  
    for i1 in range(len(var['param_to_fit'])):
        var['simulation_input'][var['param_to_fit'][i1]] = ParTry[i1]

    #generate simulation
    DataSimul = var['simulation_func'](var['simulation_input'])
   
    #calculate distance
    dist = var['distance_func'](DataSimul, var)

    if dist > 0:

        total_time = time.time() - time1
        return dist, total_time, var
    else:
        print 'dist = ' + str(dist)
        print 'DataSimul = ' + str(DataSimul)
        print 'ParTry = ' + str(ParTry)

        op1=open('simulation.dat', 'w')
        for line in DataSimul:
            for item in line:
                op1.write( str(item) + '    ')
            op1.write('\n')
        op1.close() 
        raise ValueError('Found negative distance value!')


def main():
  print(__doc__)

if __name__=='__main__':
  main()    

    




