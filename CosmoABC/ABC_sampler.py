#!/usr/bin/env python

"""
Approximate Bayesian Computation sampler.
"""

__author__ = "E. E. O. Ishida, S. D. P. Vitenti, M. Penna-Lima, R. S. de Souza, J. Cisewski, A. M. M. Trindade, V. C. Busti, E. Cameron"
__maintainer__ = "E. E. O. Ishida"
__copyright__ = "Copyright 2015"
__version__ = "0.1.9"
__email__ = "emilleishida@gmail.com"
__status__ = "Prototype"
__license__ = "GPL"


import numpy
import time
import sys

from scipy.stats import multivariate_normal
from scipy.interpolate import Rbf
from scipy.stats.mstats import mquantiles
from scipy.stats import norm

from statsmodels.stats.weightstats import DescrStatsW

from multiprocessing import Pool 

from ABC_functions import SelectParamInnerLoop, SetDistanceFromSimulation, DrawAllParams

################################################

class ABC(object):

    def __init__(self, params=False):

        """
        Constructor for the ABC class.
        This contains content related to:

        1. User input
        2. Acceptance/rejection ABC algorithm.
        3. Update weights.
        4. Iteration between particle systems. 
        5. Convergence


        Instance variables:
        -------------------------

        params		-	Complete set of input parameters
        data		-	"Real" data catalog.
	simulation	-	Simulation function.
	distance	-	Distance function.
        prior		-	Prior distribution function.
	delta		-	Parameter for convergence criteria.
	M		-	Number of particles in each particle system.
	qthreshold	-	Quantile for choosing subsequent distance thresholds.
        T		-	Total number of particle system at convergence
	

	Method attributes
	------------------------

	BuildFirstPSystem	-	Build first particle system
	BuildPSystem		-	Build subsequent particle system
	UpdateWeights		-	Update weights 
	fullABC			-	Run full ABC algorithm
        ContinueStoppedRun	-	Continue algorithm using previously sttoped run
        """ 

        self.params = params                                    #all parameters 
        self.data = params['dataset1']                          #set of 2-dimensional arrays of "real" catalog.         
        self.simulation = params['simulation_func']             #function which performs the simulation
        self.distance = params['distance_func']     	        #distance function
        self.prior = params['prior_func']                       #list of prior distribution function     
        self.delta = params['delta']                            #convergence criteria 
        self.M	= params['M']                                   #number of elements in each particle system
        self.qthreshold	= params['qthreshold']                  #quantile to define the distance threshold for subsequent particle system

        self.params['dist_dim'] = len(params['dataset1'][0]) + 1

        #check length of observed data
        if isinstance(self.data, bool):
            raise IOError('No real data catalog provided.')

        if len(self.data) < 2 :
            raise IOError('Real data catalog too short. Table must contain more than 1 element.')

        try:
            len(self.data[0])
        except TypeError:
            raise TypeError('Real data catalog must be at least 2 dimensional.')
                
        #check minimum keywords in params
        self.min_keys = ['param_to_sim', 'param_to_fit', 'prior_par', 'param_lim', 'M', 'qthreshold', 'delta','file_root']  
        for item in self.min_keys:
            if item not in self.params.keys():
                raise IOError('Keyword ' + str( item ) + '  is missing from inputed dictionary of parameters (params)!') 

        #check simulation function
        if not self.simulation:
            raise IOError('Please, provide a valid simulation function. \n See file "simulation.py".')

        #check prior function
        if not self.prior:
            raise IOError('Please provide a valid prior distribution function for each variable parameter. \n "See file priors.py"')

        #check distance function
        if not self.distance:
            raise IOError('Please provide a valid distance function. \n See file "distances.py"')

    def BuildFirstPSystem(self, output=True):
        """
        Build the first particle system, storing the parameter values satisfying distance threshold. 

        :param 	output: optional, boolean (choose to write output data file, default is True)
                 

        :returns:	[0] array containing first particle system
	                collumns -> [ param_to_fit[0], param_to_fit[1], .., param_to_fit[n], distance_from_dataset1 ] 	

		        [1] integer -> total number of draws necessary to build the first particle system
        """

        print 'Building first particle system:'

        pool = Pool(self.params['ncores'])

        time_ini = time.time()
        args = [self.params for item in xrange(self.params['Mini'])]
        dist = pool.map(SetDistanceFromSimulation, args) 
        #for i in args:
        #    print 'Calculated ' + str(args.index(i)) + 'from ' + str(len(args))

        pool.close()
        time_end = time.time() - time_ini

        total_time = time_end/self.M

        #initiate variables to store total number of draws and surviving parameters 
        theta = []    

        for line in dist:
            theta_t = [line[2]['simulation_params'][item] for item in line[2]['param_to_fit']]
                       
            for elem in line[0]:
                theta_t.append( elem )
            theta_t.append(str(self.params['Mini']/self.M))
            theta_t.append(total_time)
            theta.append(theta_t)
       
               
        #choose smaller distance 
        d1 = numpy.array([ numpy.sqrt(sum(line[j]**2 
                         for j in xrange(len(self.params['param_to_fit']),len(self.params['param_to_fit']) 
                                  + self.params['dist_dim']))) for line in theta])

        d1B = list(d1) 
        d1.sort()

        indx = [d1B.index(item) for item in d1]

        theta_new = [theta[elem] for elem in indx[:self.M]]    

        #write particle system to file
        if output:
            op = open(self.params['file_root'] + '0.dat', 'w')
            for item in self.params['param_to_fit']:
                op.write(item  + '    ')

            for i2 in xrange(self.params['dist_dim']):
                op.write('distance' + str(i2 + 1) + '    ')    
 
            op.write('NDraws    time       ')

            for i2 in xrange(self.params['dist_dim']):
                op.write('dist_threshold' + str(i2 + 1) + '    ')
            op.write('\n')

            for line in theta_new:
                for elem in line:
                    op.write(str( elem )  + '    ')
                for i3 in xrange(self.params['dist_dim']):
                    op.write(str( max(numpy.array(theta_new)[:,-self.params['dist_dim'] + i3])) + '    ')
                op.write('\n')
            op.close()

        #determine initial weights
        W1 = [1.0/self.M for i2 in xrange(self.M)]

        op2 = open(self.params['file_root'] + '0weights.dat', 'w')
        for item in W1:
            op2.write(str(item) + '\n')
        op2.close()
 
        return numpy.array(theta_new), self.params['Mini']

    def BuildPSystem(self, previous_particle_system, W, t):
        """
        Build particle system. 

        input: 	previous_particle_system -> model parameters surviving previous distance threshold
	       	collumns -> [ param_to_fit[0], param_to_fit[1], .., param_to_fit[n], distance_from_dataset1 ] 	

                W -> vector of weights

                t -> particle system index

        output: updated particle system
                collumns -> -> [ surviving_model_parameters, distance, number_necessary_draws, computational_time (s), distance_threshold ]

        """
        
        #calculate weighted covariance matrix from previous particle system
        ds = DescrStatsW(previous_particle_system[:,:len(self.params['param_to_fit'])], weights=W)
        cov1 = ds.cov
    
        var = {}
        var['params'] = self.params
        var['W'] = W
        var['previous_particle_system'] = previous_particle_system
        var['previous_cov_matrix'] = cov1

        args = [var for j in xrange(self.params['M'])]

        pool = Pool(self.params['ncores'])

        surv_param = pool.map(SelectParamInnerLoop, args)

        pool.close()

        #begin writing output file
        op = open(self.params['file_root'] + str(t) + '.dat' , 'w')
        for item in self.params['param_to_fit']:
            op.write(item  + '    ' )
        for jj in xrange(self.params['dist_dim']):
            op.write('distance' + str(jj) + '    ')     
        op.write('NDraws    time    ')
        for kk in xrange(self.params['dist_dim']): 
            op.write('dist_threshold' + str(kk + 1) + '    ')
        op.write('\n')
        for line in surv_param:
            for elem in line:
                op.write(str( elem ) + '    ')
            op.write('\n')
        op.close()

        return numpy.array(surv_param)


    def UpdateWeights(self, W, previous_particle_system, current_particle_system, output=True):   
        """
        Update weights given new particle system.

        input: 	W ->  vector of current weights

                previous_particle_system -> model parameters surviving previous distance threshold
	       	collumns -> [ param_to_fit[0], param_to_fit[1], .., param_to_fit[n], distance_from_dataset1 ] 	

                current_particle_system -> model parameters surviving current distance threshold
	       	collumns -> [ param_to_fit[0], param_to_fit[1], .., param_to_fit[n], distance_from_dataset1 ] 	

        output:   vector of updated weights      
                 
        """

        print 'update weights'

        #calculate weighted covariance matrix from previous particle system
        ds = DescrStatsW(previous_particle_system[:,:len(self.params['param_to_fit'])], weights=W)
        cov1 = ds.cov

        new_weights = []

        #determine prior distributions
        distributions = [self.prior[i1](self.params['prior_par'][i1], self.params['param_lim'][i1], func=True) 
                                        for i1 in xrange(len(self.params['param_to_fit']))]

        for i4 in range(len(current_particle_system)):
            nominator = numpy.prod([distributions[i2].pdf(current_particle_system[i4][i2])  
                                   for i2 in xrange(len(self.params['param_to_fit']))])

            denominator = sum(W[i3]*multivariate_normal.pdf(current_particle_system[i4][:len(self.params['param_to_fit'])], 
                              previous_particle_system[i3][:len(self.params['param_to_fit'])], cov=cov1)  
                                                               for i3 in xrange(len(W)))

            new_weights.append(nominator/denominator)

        final_weights = [item/sum(new_weights) for item in new_weights]

        if output == True:
            op = open(self.params['file_root'] + str(self.T) + 'weights.dat', 'w')
            for item in final_weights:
                op.write(str(item) + '\n')
            op.close()

        return final_weights
              


    def fullABC(self, build_first_system=False):
        """
        Run complete ABC sampler algorithm. 

        input:	root_file_name -> root of file name to be used in all runs (string)

		build_first_system (optional) -> boolean  (read or generate first particle system). Default is False.
   	
        output:	particle systems and corresponding weights written in data files.
        """

        #determine initial weights
        W1 = [1.0/self.M for i2 in xrange(self.M)]

        if build_first_system == True:
            #build first particle system
            sys0 = self.BuildFirstPSystem()

        #read first particle system from file
        op = open(self.params['file_root'] + '0.dat', 'r')
        lin = op.readlines()
        op.close()

        t1 = [elem.split() for elem in lin[1:]]

        sys1 = numpy.array([numpy.array([float(line[i1]) 
                           for i1 in xrange(len(self.params['param_to_fit']) + self.params['dist_dim'])]) 
                                            for line in t1 ])

        #determine number of draws in previous particle system generation
        K =  sum(int(line[len(self.params['param_to_fit']) + self.params['dist_dim']]) for line in t1)
        print 'number of draws PS0 = ' + str(K)

        #initiate iteration counter
        t = 0
        K = self.M

        while float(self.M)/K > self.delta:      

            t = t + 1

            self.T = t 

            sys_new = self.BuildPSystem(sys1, W1, t)
        
            W2 = self.UpdateWeights(W1, sys1, sys_new)

 
            K = sum(sys_new[:, len(self.params['param_to_fit']) + self.params['dist_dim']])

            del sys1, W1

            sys1 = sys_new
            W1 = W2

            del sys_new, W2 

            print ' finished PS ' + str(t) + ',    convergence = ' + str(float(self.M)/K)
           
        self.T = t

        
    def  ContinueStoppedRun(self, t):
        """
        Continue ABC sampler algorithm from a specific time-step (run). 

        input: 	t -> index of last completed particle system (int)

		root_file_name -> root of file name to be used in all subsequent runs (string)
	
	output:	subsequent particle systems and corresponding weights written in data files.
        """
 

        op = open(self.params['file_root'] + str(t) + '.dat', 'r' )
        lin = op.readlines()
        op.close()

        t1 = [elem.split() for elem in lin]
        
        sys1 = numpy.array([numpy.array([float(line[i1]) 
                          for i1 in xrange(len(self.params['param_to_fit']) + self.params['dist_dim'])]) 
                                          for line in t1[1:]])
        
        #determine number of draws in previous particle system generation
        K =  sum(int(line[list(t1[0]).index('NDraws')]) for line in t1[1:])
        print 'number of draws PS' + str(t) + ' = ' + str(K)
        
        if t > 0:        
            W1 = numpy.loadtxt(self.params['file_root'] + str(t) + 'weights.dat')
        elif t == 0:
            W1 = [1.0/self.M for i2 in xrange(self.M)]
    

        while float(self.M)/K > self.delta:

            t = t + 1

            self.T = t

            sys_new = self.BuildPSystem(sys1, W1, t)
        
            W2 = self.UpdateWeights(W1, sys1, sys_new)

 
            K = sum(sys_new[:, len(self.params['param_to_fit' ]) + self.params['dist_dim']])

            del sys1, W1

            sys1 = sys_new
            W1 = W2

            del sys_new, W2 


            print ' finished PS' + str(t) + ',    convergence = ' + str(float(self.M)/K)
        
                  
        self.T = t



def main():
  print(__doc__)



if __name__=='__main__':
  main()    

    




