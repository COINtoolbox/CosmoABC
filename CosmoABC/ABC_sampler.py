#!/usr/bin/env python

"""
Approximate Bayesian Computation library.
"""

__author__ = "E. E. O. Ishida, S. D. P. Vitenti, M. Penna-Lima,  R. S. de Souza, J. Cisewski, E. Cameron, V. C. Busti"
__maintainer__ = "E. E. O. Ishida"
__copyright__ = "Copyright 2015"
__version__ = "0.1.4"
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



################################################


class ABC( object ):

    def __init__( self, params=False ):

        """
        Constructor for the ABC class.
        This contains content related to:

        1. User input
        2. Comparison of distances between catalogs. 
        3. Acceptance/rejection ABC algorithm.
        4. Update weights.
        5. Iteration between particle systems. 
        6. Convergence


        Instance variables:
        -------------------------

        data		-	"Real" data catalog.
	simulation	-	Simulation function.
	params		-	Complete set of input parameters
	prior		-	Prior distribution function.
	distance	-	Distance function.
	delta		-	Parameter for convergence criteria.
	s		-	Smooth parameter
	M		-	Number of particles in each particle system.
        epsilon1	-	Distance threshold for first particle system.
	qthreshold	-	Quantile for choosing subsequent distance thresholds.
        T		-	Total number of particle system at convergence
	

	Method attributes
	------------------------

	DataSimul	-	Simulated catalog
	dist		-	Distance between simulated and observed catalogs
	K		-	Number of draws necessary to accept one particle
	firstPS		-	Initial particle system 
        epsilon		-	Distance threshold
        theta_t_try	-	Parameter set surviving acceptance/rejection algorithm
        particle_system	-	Complete particle system for iterations > 0
        cov1		-	Weighted covariance matrix
        updated_weights	-	Updated weights 
         
        """ 

        self.params = params                                    #all parameters 
        self.data = params['dataset1']                          #set of 2-dimensional arrays of "real" catalog.         
        self.simulation = params['simulation_func']             #function which performs the simulation
        self.distance = params['distance_func']     	        #distance function
        self.prior = params['prior_func']                       #list of prior distribution functions

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
        self.min_keys = ['simulation_params', 'param_to_fit', 'prior_par', 'param_lim', 'M', 'epsilon1', 'qthreshold', 'delta','s', 'file_root']  
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
           
        #define parameters separately so inputed values are not changed       
        self.delta = params['delta']             #convergence criteria 
        self.s = params['s']                     #smooth parameter (scalar)
        self.M	= params['M']                    #number of elements in each particle system
        self.epsilon1 = params['epsilon1']       #list distance threshold for first particle system
        self.qthreshold	= params['qthreshold']   #quantile to define the distance threshold for subsequent particle system
    


    def DrawAllParams(self):
        """
        Draw complete set of  parameters from prior.

        :returns: array of parameters sampled from the prior
        """

        pars = []
        for j in range(len(self.params['param_to_fit'])):
            p1 = self.prior[j]( self.params['prior_par'][j], self.params['param_lim'][j])   
            pars.append(p1)

        return numpy.array(pars)

   
    def SetDistanceFromSimulation(self):
        """
        Draw cosmological parameter values from prior, generate simulation and calculate distance from a given comparison  catalog. 
        In the context of Ishida et al., 2015 that translates into calculating the distance from the real to the simulated cataloges.

        :returns: scalar (distance between dataset1 and dataset2)
    
        """

        #draw parameters from prior
        ParTry = self.DrawAllParams()
  
        for i1 in range(len(self.params['param_to_fit'])):
            self.params['simulation_params'][self.params['param_to_fit'][i1]] = ParTry[i1]

        #generate simulation
        DataSimul = self.simulation(self.params['simulation_params'])
   
        #calculate distance
        dist = self.distance(DataSimul, self.params)

        if dist > 0:
            return dist
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
    

    def BuildFirstPSystem(self, output=True):
        """
        Build the first particle system, storing the parameter values satisfying distance threshold. 

        :param 	output: optional, boolean (choose to write output data file, default is True)
                 

        :returns:	[0] array containing first particle system
	                collumns -> [ param_to_fit[0], param_to_fit[1], .., param_to_fit[n], distance_from_dataset1 ] 	

		        [1] integer -> total number of draws necessary to build the first particle system
        """

        print 'Building first particle system:'

        #initiate variables to store total number of draws and surviving parameters 
        theta = []    
   
        while len( theta ) < self.params['Mini']:
 
            dist = [10**10]        
            time1 = time.time()
            dist = self.SetDistanceFromSimulation()
  
            theta_t = [self.params['simulation_params'][item] for item in self.params['param_to_fit']] 

            if dist[0] > 0:

                total_time = time.time() - time1

                for item in dist:  
                    theta_t.append(item)
                theta_t.append(str(self.params['Mini']/self.M))
                theta_t.append(total_time)
                theta.append(theta_t)
        
                if self.params['screen'] == True:
                    print '        particle index = ' + str(len(theta)) 
                    for ii in xrange(len(dist)) :
                        print '          distance' + str(ii) + '=' + str(dist[ii])
               
        #choose smaller distance 
        d1 = numpy.array(theta)[:,len(self.params['param_to_fit'])]   
        d1B = list(d1) 
        d1.sort()

        indx = [d1B.index(item) for item in d1]

        theta_new = [theta[elem] for elem in indx[:self.M]]    

        #write particle system to file
        if output:
            op = open(self.params['file_root'] + '0.dat', 'w')
            for item in self.params['param_to_fit']:
                op.write(item  + '    ')

            for i2 in xrange(len(self.params['epsilon1'])):
                op.write('distance' + str(i2 + 1) + '    ')    
 
            op.write('NDraws    time       ')

            for i2 in xrange(len(self.params['epsilon1'])):
                op.write('dist_threshold' + str(i2 + 1) + '    ')
            op.write('\n')

            for line in theta_new:
                for elem in line:
                    op.write(str( elem )  + '    ')
                for i3 in xrange(len( self.epsilon1 )):
                    op.write(str( max(numpy.array(theta_new)[:,-len(self.epsilon1) + i3])) + '    ')
                op.write('\n')
            op.close()
 
        return numpy.array(theta_new), self.params['Mini']


    def SelectParamInnerLoop( self, previous_particle_system, W, previous_cov_matrix ):
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

        dist = [10**10 for i in xrange( len( self.params['epsilon1'] ) )]
        K = 0

        #time marker
        time_start = time.time()

        #determine distance threshold
        epsilon = [mquantiles(previous_particle_system[:, kk ], prob=self.qthreshold)[0] 
                      for kk in xrange(len(self.params['param_to_fit']), 
                                       len(self.params['param_to_fit']) + len(self.params['epsilon1']))]

        flag = [False for ii in xrange(len(self.params['epsilon1']))]

        while False in flag:               
 
            #update counter
            K = K + 1 
  
            #draw model parameters to serve as mean 
            index_theta0 = numpy.random.choice(xrange(len(W)), p=W)
            theta0 = numpy.atleast_2d(previous_particle_system[index_theta0][:len( self.params['param_to_fit'])])

            #initialize boolean parmeter vector
            theta_t = [False for i in xrange(len(self.params['param_to_fit']))]
           
            while False in theta_t:
          
                #draw model parameter values for simulation
                mvn = multivariate_normal.rvs(theta0.flatten(), previous_cov_matrix)

                try:
                    len(mvn)
                    theta_t_try = list(mvn)

                except TypeError:
                    theta_t_try = [mvn]

                theta_t = []
                for k1 in xrange(len(self.params['param_to_fit'])):

                    if  theta_t_try[k1] >= self.params['param_lim'][k1][0] and theta_t_try[k1] < self.params['param_lim'][k1][1]:
                        theta_t.append(True)

                    else:
                        theta_t.append(False)         


            #update parameter values in dictionary
            for i1 in range(len(self.params['param_to_fit'])):
                self.params['simulation_params'][self.params['param_to_fit'][i1]] = theta_t_try[i1]

            #generate simulation
            DataSimul = self.simulation(self.params['simulation_params'])
      
            #calculate distance
            dist = self.distance(DataSimul, self.params)
              
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

        return theta_t_try


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

        #begin writing output file
        op = open(self.params['file_root'] + str(t) + '.dat' , 'w')
        for item in self.params['param_to_fit']:
            op.write(item  + '    ' )
        for jj in xrange(len(self.params['epsilon1'])):
            op.write('distance' + str(jj) + '    ')     
        op.write('NDraws    time    ')
        for kk in xrange(len(self.params['epsilon1'])): 
            op.write('dist_threshold' + str(kk + 1) + '    ')
        op.write('\n')

        particle_system = []

        for j in xrange(self.M):

            surv_param = self.SelectParamInnerLoop(previous_particle_system,  W, cov1)

            particle_system.append(surv_param)
  
            if self.params['screen'] == True:
                print 'particle index = ' + str(j) + ',    number of draws = ' + \
                       str(surv_param[len(self.params['param_to_fit']) + len(qself.params['epsilon1'])])
 
                for thr in xrange(len(self.epsilon1)):
                    ',    distance' + str(thr) + '=' + str(surv_param[len(self.params['param_to_fit']) + thr])

            for elem in surv_param:
                op.write(str( elem ) + '    ')
            op.write('\n')
        op.close()

        return numpy.array(particle_system)


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

        if build_first_system == True:
            #build first particle system
            sys0 = self.BuildFirstPSystem()
            
        #read first particle system from file
        op = open(self.params['file_root'] + '0.dat', 'r')
        lin = op.readlines()
        op.close()

        t1 = [elem.split() for elem in lin[1:]]

        sys1 = numpy.array([numpy.array([float(line[i1]) 
                           for i1 in xrange(len(self.params['param_to_fit']) + len(self.params['epsilon1']))]) 
                                            for line in t1 ])

        #determine number of draws in previous particle system generation
        K =  sum(int(line[len(self.params['param_to_fit']) + len(self.params['epsilon1'])]) for line in t1)
        print 'number of draws PS0 = ' + str(K)

        #determine initial weights
        W1 = [1.0/self.M for i2 in xrange(self.M)]

        #initiate iteration counter
        t = 0
        K = self.M

        while float(self.M)/K > self.delta:

            self.T = t

            t = t + 1

            sys_new = self.BuildPSystem(sys1, W1, t)
        
            W2 = self.UpdateWeights(W1, sys1, sys_new)

 
            K = sum(sys_new[:, len(self.params['param_to_fit']) + len( self.params['epsilon1'])])

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
                          for i1 in xrange(len(self.params['param_to_fit']) + len(self.params['epsilon1']))]) 
                                          for line in t1[1:]])
        
        #determine number of draws in previous particle system generation
        K =  sum(int(line[list(t1[0]).index('NDraws')]) for line in t1[1:])
        print 'number of draws PS' + str(t) + ' = ' + str(K)
        
        if t > 0:        
            W1 = numpy.loadtxt(self.params['file_root'] + str(t) + 'weights.dat')
        elif t == 0:
            W1 = [1.0/self.M for i2 in xrange(self.M)]
    

        while float(self.M)/K > self.delta:

            self.T = t

            t = t + 1

            sys_new = self.BuildPSystem(sys1, W1, t)
        
            W2 = self.UpdateWeights(W1, sys1, sys_new)

 
            K = sum(sys_new[:, len(self.params['param_to_fit' ]) + len(self.params['epsilon1'])])

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

    




