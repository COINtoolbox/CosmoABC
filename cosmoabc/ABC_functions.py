#!/usr/bin/env python3

"""
Functions to optimize parallelization.
"""
import numpy as np
import time
import sys
import os

from inspect import isfunction

from scipy.stats import multivariate_normal
from scipy.interpolate import Rbf
from scipy.stats.mstats import mquantiles
from scipy.stats import norm

from statsmodels.stats.weightstats import DescrStatsW

from .distances import distance_quantiles, distance_GRBF 
from .priors import flat_prior, gaussian_prior, beta_prior


################################################

def get_cores():
    """
    Ask the user to input number of cores.
    """
    cores = input('Please enter the number of cores: ')

    return int(cores) 

def read_input(filename):
    """
    Read user input from file and construct initial dictionary parameter. 

    input: filename (string) -> user input file parameter 

    output: dictionary with formated user choices
    """

    #read user input data
    op1 = open(filename, 'r')
    lin1 = op1.readlines()
    op1.close()

    data1 = [elem.split() for elem in lin1]     

    #store options in params dictionary
    params_ini = dict([(line[0], line[2:]) 
                        for line in data1 if len(line) > 1])
    
    params = {}
    params['path_to_obs'] = params_ini['path_to_obs'][0] 
    params['param_to_fit'] = [params_ini['param_to_fit'][i] 
                             for i in range(
                                 params_ini['param_to_fit'].index('#'))]
   
    params['param_to_sim'] = [params_ini['param_to_sim'][i] 
                             for i in range(
                                 params_ini['param_to_sim'].index('#'))]
    params['npar'] = len(params['param_to_fit'])
    params['M'] = int(params_ini['M'][0])
    params['Mini'] = int(params_ini['Mini'][0])
    params['qthreshold'] = float(params_ini['qthreshold'][0])
    params['delta'] = float(params_ini['delta'][0])
    params['file_root'] = params_ini['file_root'][0]  
    params['screen'] = bool(int(params_ini['screen'][0]))
    params['ncores'] = int(params_ini['ncores'][0])
    params['prior_func'] = [params_ini['prior_func'][i] 
                             for i in range(
                                 params_ini['prior_func'].index('#'))]
    params['dist_dim'] = int(params_ini['dist_dim'][0])
    
    #functions
    from distances import distance_GRBF
    dispatcher = {'flat_prior': flat_prior, 
                  'gaussian_prior': gaussian_prior, 'beta_prior': beta_prior, 
                  'distance_quantiles': distance_quantiles, 
                  'distance_GRBF': distance_GRBF}
    
    if params_ini['simulation_func'][0] in dispatcher:
        params['simulation_func'] = dispatcher[params_ini['simulation_func'][0]]

    if len(params['prior_func']) == params['npar']:

        prior_dic = {}
        prior_dic['sequence'] = params['param_to_fit']
        for i1 in range(params['npar']):
            el = params['param_to_fit'][i1]
            prior_dic[el] = {}
            if el + '_lim' in params_ini.keys():
                prior_dic[el]['min'] = float(params_ini[el + '_lim'][0])
                prior_dic[el]['max'] = float(params_ini[el + '_lim'][1])

            indx = params_ini[el + '_prior_par_name'].index('#')
            for name in params_ini[el + '_prior_par_name'][:indx]:
                indx = params_ini[el + '_prior_par_name'].index(name)
                prior_dic[el][name] = float(params_ini[el + 
                                            '_prior_par_val'][indx])

            if params['prior_func'][i1] == 'gaussian_prior':    
                prior_dic[el]['func'] = dispatcher['gaussian_prior']

            elif params['prior_func'][i1] == 'flat_prior':
                prior_dic[el]['func'] = dispatcher['flat_prior']
                
            elif params['prior_func'][i1] == 'beta_prior':
                prior_dic[el]['func'] = dispatcher['beta_prior']
               
            else:
                prior_dic[el]['func'] = params['prior_func'][i1]

        params['prior'] = prior_dic
                    
    else:
        raise ValueError('Number of prior functions does not ' +
                          'match number of parameters!') 


    #fiducial extra parameters
    sim_par = {}  
    for item in params['param_to_sim']:
        try:
            par = float(params_ini[item][0])
            sim_par[item] = par
        except ValueError:
            sim_par[item] = params_ini[ item ][0] 

    if params_ini['simulation_func'][0] == 'numcosmo_sim_cluster':

        try: 
            from gi.repository import NumCosmo as Nc
            from sim_NumCosmo_cluster import NCountSimul
            from sim_NumCosmo_cluster import ChooseParamsInput 
            from sim_NumCosmo_cluster import numcosmo_sim_cluster
        except ImportError:
            raise ImportError('You must have NumCosmo running to use the ' + 
                              'sim_NumCosmo simulation!' +
                              '\n Please check your NumCosmo instalation.')

        sim_par["OL"] = 1. - sim_par["Om"] - sim_par["Ob"]
        Cosmo=ChooseParamsInput()
        Cosmo.params = sim_par 
        params['simulation_input'] = Cosmo.params
        params['simulation_func'] = numcosmo_sim_cluster
        
    else:
        params['simulation_input'] = sim_par


    #check if ``observer'' data already exists, simulate in case negative
    if params['path_to_obs'] != 'None':
        #read observed data
        op2 = open(params_ini['path_to_obs'][0], 'r')
        lin2 = op2.readlines()
        op2.close()

        data2 = [elem.split() for elem in lin2[1:]]
        params['dataset1'] = np.array([[float(item) for item in line] 
                                        for line in data2])  

    elif 'simulation_func' in params.keys():
        params['dataset1'] = params['simulation_func'](
                                    params['simulation_input']
                                    )
    
    if params_ini['distance_func'][0] in dispatcher.keys():
        params['distance_func'] = dispatcher[params_ini['distance_func'][0]]

        if params['distance_func'] == distance_GRBF:
            from distances import GRBF, logf, norm_GRBF, prep_GRBF
            from distances import distance_GRBF
            params['s'] = float(params_ini['s'][0])
            params = prep_GRBF(params)

        elif (('dataset1' in params) and 
              (params['distance_func'] == distance_quantiles)):
            from distances import summ_quantiles
            params['dist_dim'] = len(params['dataset1'][0]) + 1
            params['quantile_nodes'] = int(params_ini['quantile_nodes'][0])
            params = summ_quantiles(params['dataset1'], params)

    for key in params_ini.keys():
        if key not in params.keys():
            params[key] = params_ini[key]

    return params

def SelectParamInnerLoop(var1):
    """
    Draw model parameters based on previous particle system and return those 
    satisfying distance threshold.

    input: var1 -> dictionary of input parameters	

    output: vector of [surviving_model_parameters, distance, 
                       number_necessary_draws, computational_time (s), 
                       distance_threshold]
    """

    try:
        dist = [10**10 for i in range(var1['params']['dist_dim'])]
        K = 0
        np.random.seed()

        #time marker
        time_start = time.time()

        #determine distance threshold
        epsilon = [
                   mquantiles(
                             var1['previous_particle_system'][:, kk], 
                             prob=var1['params']['qthreshold']
                             )[0] 
                      for kk in range(
                                      len(var1['params']['param_to_fit']), 
                                      len(var1['params']['param_to_fit']) + 
                                      var1['params']['dist_dim']
                                      )
                  ]

        flag = [False for ii in range(var1['params']['dist_dim'])]

        while False in flag:               
 
            #update counter
            K = K + 1 

            #draw model parameters to serve as mean 
            index_theta0 = np.random.choice(range(len(var1['W'])), 
                                            p=var1['W'])
            indx = len(var1['params']['param_to_fit'])
            theta0 = np.atleast_2d(var1['previous_particle_system'][index_theta0][:indx])

            #initialize boolean parmeter vector
            theta_t = [False for i in range(len(var1['params']['param_to_fit']))]
           
            while False in theta_t:
          
                #draw model parameter values for simulation
                mvn = multivariate_normal.rvs(theta0.flatten(), var1['previous_cov_matrix'])

                try:
                    len(mvn)
                    theta_t_try = list(mvn)

                except TypeError:
                    theta_t_try = [mvn]

                theta_t = []
                for par in var1['params']['param_to_fit']:
                   
                    k1 = var1['params']['param_to_fit'].index(par)
                    if (theta_t_try[k1] >= var1['params']['prior'][par]['min'] 
                    and theta_t_try[k1] <= var1['params']['prior'][par]['max']):
                        theta_t.append(True)

                    else:
                        theta_t.append(False)         
   

            #update parameter values in dictionary
            for i1 in range(len(var1['params']['param_to_fit'])):
                p = var1['params']['param_to_fit'][i1]
                var1['params']['simulation_input'][p] = theta_t_try[i1]

            #generate simulation
            DataSimul = var1['params']['simulation_func'](
                                          var1['params']['simulation_input']
                                                          )
      
            #calculate distance
            dist = var1['params']['distance_func'](DataSimul, var1['params'])
              
            #check if it satisfies distance thresholds 
            flag = [] 
            for ll in range(len( dist )):
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
            print('Number of draws = ' + str(K))

        return theta_t_try
    
    except KeyboardInterrupt:
        pass


def DrawAllParams(prior_dic):
    """
    Draw complete set of  parameters from prior.
  
    input: dictionary of input parameters
    output: array of parameters sampled from the prior
    """

    pars = []

    for element in prior_dic['sequence']:
        p1 = prior_dic[element]['func'](prior_dic[element])   
        pars.append(p1)

    return np.array(pars)
   
def SetDistanceFromSimulation(var):
    """
    Draw cosmological parameter values from prior, generate simulation 
    and calculate distance from a given comparison  catalog. 

    input: dictionary of input parameters
    output: distance between dataset1 and dataset2
    """

    try:
        time1 = time.time()

        #draw parameters from prior
        ParTry = DrawAllParams(var['prior'])
  
        for i1 in range(len(var['param_to_fit'])):
            var['simulation_input'][var['param_to_fit'][i1]] = ParTry[i1]

        #generate simulation
        DataSimul = var['simulation_func'](var['simulation_input'])
   
        #calculate distance
        dist = var['distance_func'](DataSimul, var)
        
        if isinstance(dist, float):
            check_dist = dist > 0
        else:
            dist = np.array(dist, dtype=object)
            check_dist = np.linalg.norm(dist) > 0

        if check_dist:
            total_time = time.time() - time1
            if var['screen']:
                print('Calculated distance: \n   ' + str(dist))

            result = (dist, total_time, ParTry)
 
            return result
        else:
            print('dist = ' + str(dist))
            print('DataSimul = ' + str(DataSimul))
            print('ParTry = ' + str(ParTry))

            op1=open('simulation.dat', 'w')
            for line in DataSimul:
                for item in line:
                    op1.write( str(item) + '    ')
                op1.write('\n')
            op1.close() 
            raise ValueError('Found negative distance value!')

    except KeyboardInterrupt:
        pass        

def main():
  print(__doc__)

if __name__=='__main__':
  main()    

    




