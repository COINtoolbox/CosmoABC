#!/usr/bin/env python3

import unittest 
import os
import numpy as np
import sys

from statsmodels.stats.weightstats import DescrStatsW

from cosmoabc.distances import distance_quantiles, summ_quantiles
from cosmoabc.priors import flat_prior
from cosmoabc.ABC_sampler import ABC
from cosmoabc.ABC_functions import SelectParamInnerLoop, SetDistanceFromSimulation, DrawAllParams, get_cores
from cosmoabc.plots import plot_1p, plot_2p, plot_3p


def ysim(v):

    l1 = np.random.normal(loc=v['mu'], scale=v['sigma'], size=v['n'])
    
    return np.atleast_2d(l1).T    

class TestABC(unittest.TestCase):

    def setUp(self):

        self.mu = 2.5
        self.sigma = 1.0
        self.n = 1000
        
        self.params = {}
        self.params['simulation_func'] = ysim
        self.params['simulation_input'] = {'mu': self.mu, 'sigma':self.sigma, 'n':self.n} 
        self.params['dataset1'] = self.params['simulation_func']( self.params['simulation_input'] )
        self.params['param_to_fit']=['mu', 'sigma']
        self.params['screen'] = 0
        self.params['Mini'] = 200 
        self.params['M'] = 100
        self.params['delta'] =0.1
        self.params['qthreshold'] = 0.75
        self.params['file_root'] = os.getcwd() + '/test_PS'
        self.params['distance_func'] =  distance_quantiles  
        self.params['quantile_nodes'] = 20
        self.params['split_output'] = [1]
        self.params['prior'] = {}
        self.params['prior']['sequence'] = self.params['param_to_fit']
        self.params['prior']['mu'] = {}
        self.params['prior']['mu']['func'] = flat_prior
        self.params['prior']['mu']['pmin'] = 1.0
        self.params['prior']['mu']['pmax'] = 4.0
        self.params['prior']['mu']['min'] = 1.0
        self.params['prior']['mu']['max'] = 4.0
        self.params['prior']['sigma'] = {}
        self.params['prior']['sigma']['func'] = flat_prior
        self.params['prior']['sigma']['pmin'] = 0.001
        self.params['prior']['sigma']['pmax'] = 3.0
        self.params['prior']['sigma']['min'] = 0.001
        self.params['prior']['sigma']['max'] = 3.0  

        #initiate ABC sampler
        self.sampler_ABC = ABC( self.params ) 
        self.W = [1.0/self.params['M'] for i in range( self.params['M'] )]
        self.params = summ_quantiles(self.params['dataset1'], self.params) 

    def test_DrawAllParams( self ):
         
        #draw parameters
        r1 = DrawAllParams(self.params['prior'])

        res = []
        for i1 in range(len(r1)):
            par = self.params['param_to_fit'][i1]
            if r1[i1] >= self.params['prior'][par]['min'] and r1[i1] <= self.params['prior'][par]['max']:
                res.append(True)

            else:
                res.append(False)

        self.assertEqual([True for item in r1], res)

    def test_SetDistanceFromSimulation(self):

        #set distance
        r2 = SetDistanceFromSimulation(self.params)
        
        if isinstance(r2, float):
            check_dist = r2 > 0
        else:
            check_dist = np.linalg.norm(r2[0]) > 0
        
        self.assertTrue(check_dist)

    def test_plot(self):

     
        self.params['ncores'] = get_cores()
        self.sampler_ABC.fullABC(build_first_system=True)

        if len(self.params['param_to_fit']) == 1:
            plot_1p(self.sampler_ABC.T, 'results.pdf', self.params)

        elif len(self.params['param_to_fit']) == 2:
            plot_2p(self.sampler_ABC.T, 'results.pdf', self.params) 

        elif len(self.params['param_to_fit']) == 3:
            plot_3p(self.sampler_ABC.T, 'results.pdf', self.params) 

        elif len(self.params['param_to_fit']) == 4:
            plot_4p(self.sampler_ABC.T, 'results.pdf', self.params) 
   

if __name__ == '__main__':

    unittest.main()
