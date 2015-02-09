import unittest 
import os
import numpy

from statsmodels.stats.weightstats import DescrStatsW

from CosmoABC.distances import distance_quantiles, summ_quantiles
from CosmoABC.priors import flat_prior
from CosmoABC.ABC_sampler import ABC
from CosmoABC.plots import plot_1D, plot_2D, plot_3D, plot_4D


def ysim(v):

    l1 = numpy.random.normal(loc=v['mu'], scale=v['sigma'], size=v['n'])
    
    return numpy.atleast_2d(l1).T 


class TestABC(unittest.TestCase):

    def setUp(self):

        self.mu = 2.5
        self.sigma = 0.1
        self.n = 100
        
        self.params = {}
        self.params['simulation_func'] = ysim
        self.params['simulation_params'] = {'mu': self.mu, 'sigma':self.sigma, 'n':self.n} 
        self.params['dataset1'] = self.params['simulation_func']( self.params['simulation_params'] )
        self.params['param_to_fit']=['mu', 'sigma', 'n']							
        self.params['param_lim']=[[2.0, 3.0],[0.001, 0.5],[50, 150]]	
        self.params['prior_par'] = [[2.0, 3.0],[0.001, 0.4],[60, 140]]
        self.params['screen'] = 0
        self.params['Mini'] = 200 
        self.params['s']=0.15					
        self.params['epsilon1'] = [ 0.5, 0.5]			
        self.params['M'] = 100				
        self.params['delta'] =0.2				
        self.params['qthreshold'] = 0.75
        self.params['file_root'] = os.getcwd() + '/test_PS'	
        self.params['distance_func'] =  distance_quantiles 
 	self.params['prior_func'] = [ flat_prior, flat_prior, flat_prior ]	

        #initiate ABC sampler
        self.sampler_ABC = ABC( self.params ) 

        self.W = [1.0/self.params['M'] for i in xrange( self.params['M'] )]
        self.params = summ_quantiles(self.params['dataset1'], self.params)
      
    def test_DrawAllParams( self ):
         
        #draw parameters
        r1 = self.sampler_ABC.DrawAllParams()

        res = []
        for i1 in xrange(len(r1)):
            if r1[i1] >= self.params['param_lim'][i1][0] and r1[i1] < self.params['param_lim'][i1][1]:
                res.append(True)

            else:
                res.append(False)

        self.assertEqual([True for item in r1], res)

    def test_SetDistanceFromSimulation(self):

        #set distance
        r2 = self.sampler_ABC.SetDistanceFromSimulation()

        self.assertTrue(r2 >= 0)

    def test_SelectParamInnerLoop(self):

        PS0 = self.sampler_ABC.BuildFirstPSystem(output=True)

        vv = numpy.array([[float(item) 
                           for item in line[:len(self.params['param_to_fit']) + len(self.params['epsilon1'])]] 
                                       for line in  PS0[0]])

        #determine previous covariance matrix
        ds = DescrStatsW(vv[:,:len(self.params['param_to_fit'])], weights=self.W)

        cov1 = ds.cov

        #Select Parameters for subsequent loop
        r4 = self.sampler_ABC.SelectParamInnerLoop(vv, self.W, cov1)

        res = []
        for i1 in xrange(len(self.params['param_to_fit'])):
            if r4[i1] >= self.params['param_lim'][i1][0] and r4[i1] < self.params['param_lim'][i1][1]:
                res.append(True)

            else:
                res.append(False)

        self.assertEqual([True for item in r4[:len(self.params['param_to_fit'])]], res)

    def test_BuildPSystem(self): 

        #build second particle system
        PS0 = self.sampler_ABC.BuildFirstPSystem(output=True)
        vv = numpy.array([[float(item) 
                        for item in line[:len(self.params['param_to_fit']) + len(self.params['epsilon1'])]] 
                                    for line in  PS0[0]])
  
        PS1 = self.sampler_ABC.BuildPSystem(vv, self.W, 1)

        self.assertTrue(len(PS1) == self.params['M'])

    def test_UpdateWeights(self):

        PS0 = self.sampler_ABC.BuildFirstPSystem(output=True)
        self.sampler_ABC.T = 1

        vv = numpy.array([[float(item) 
                        for item in line[:len(self.params['param_to_fit']) + len(self.params['epsilon1'])]] 
                                    for line in  PS0[0]])
        PS1 = self.sampler_ABC.BuildPSystem(vv, self.W, self.sampler_ABC.T)

        #update weights
        r6 = self.sampler_ABC.UpdateWeights(self.W, vv, PS1, output=True)

        self.assertTrue(len(r6) == self.params['M'])

    def test_plot(self):

        self.sampler_ABC.fullABC(build_first_system=False)

        if len(self.params['param_to_fit']) == 1:
            plot_1D(self.sampler_ABC.T, 'results.pdf', self.params)

        elif len(self.params['param_to_fit']) == 2:
            plot_2D(self.sampler_ABC.T, 'results.pdf', self.params) 

        elif len(self.params['param_to_fit']) == 3:
            plot_3D(self.sampler_ABC.T, 'results.pdf', self.params) 

        elif len(self.params['param_to_fit']) == 4:
            plot_4D(self.sampler_ABC.T, 'results.pdf', self.params) 
    
    def test_continueStoppedRun(self):

        self.sampler_ABC.ContinueStoppedRun(2) 

if __name__ == '__main__':
    unittest.main()
