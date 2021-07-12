#!/usr/bin/env python3

"""
Approximate Bayesian Computation sampler.
"""

import numpy as np
import time
import sys
import os 

from scipy.stats import multivariate_normal
from scipy.interpolate import Rbf
from scipy.stats.mstats import mquantiles
from scipy.stats import norm

from statsmodels.stats.weightstats import DescrStatsW

from multiprocessing import Pool 

from .ABC_functions import SelectParamInnerLoop, SetDistanceFromSimulation 
from .ABC_functions import DrawAllParams

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
     	M		    -	Number of particles in each particle system.
     	qthreshold	-	Quantile for choosing subsequent distance thresholds.
        T		    -	Total number of particle system at convergence
	

	Method attributes
	------------------------

	BuildFirstPSystem	-	Build first particle system
	BuildPSystem		-	Build subsequent particle system
	UpdateWeights		-	Update weights 
	fullABC			-	Run full ABC algorithm
    ContinueStoppedRun	-	Continue algorithm from previous  one
        """ 

        self.params = params                             
        self.data = params['dataset1']                       
        self.simulation = params['simulation_func']      
        self.distance = params['distance_func']     	        
        self.prior = params['prior']                             
        self.delta = params['delta']                           
        self.M	= params['M']                                   
        self.qthreshold	= params['qthreshold']                  

        #check length of observed data
        if isinstance(self.data, bool):
            raise IOError('No real data catalog provided.')

        if len(self.data) < 2 :
            raise IOError('Real data catalog too short. '  +
                          'Table must contain more than 1 element.')

        try:
            len(self.data[0])
        except TypeError:
            raise TypeError('Real data catalog must be at ' + 
                            'least 2 dimensional.')
                
        #check minimum keywords in params
        self.min_keys = ['simulation_input', 'param_to_fit', 'prior', 
                         'M', 'qthreshold', 'delta','file_root']  
        for item in self.min_keys:
            if item not in self.params.keys():
                raise IOError('Keyword ' + str(item) + 
                              '  is missing from inputed dictionary of ' + 
                              'parameters (params)!') 

        #check simulation function
        if not self.simulation:
            raise IOError('Please, provide a valid simulation function.' + 
                          '\n See file "simulation.py".')

        #check prior function
        if not self.prior:
            raise IOError('Please provide a valid prior information ' + 
                          'for each variable parameter. \n ' + 
                          '"See file priors.py"')

        #check distance function
        if not self.distance:
            raise IOError('Please provide a valid distance function. ' + 
                          '\n See file "distances.py"')

    def BuildFirstPSystem(self, output=True):
        """
        Build the first particle system, storing the parameter values 
        satisfying distance threshold. 

        input:  output (optional) -> choose to write output data file, 
                                     default is True
                 
        output:	[0] array containing first particle system
	        collumns -> [ param_to_fit[0], param_to_fit[1], .., 
                              param_to_fit[n], distance_from_dataset1 ] 	

		[1] integer -> total number of draws necessary to build 
                               the first particle system
        """

        print('Building first particle system:')

        theta = []

        #check if there are previous partial results
        temp_files = [self.params['file_root'] + '0_p' + str(part) + '.dat'
                      for part in range(int(self.params['split_output'][0]))]

        file_list = os.listdir(os.getcwd())
        if temp_files[0] in file_list:
           for partial_calc in range(int(self.params['split_output'][0])):
               if temp_files[partial_calc] in file_list:
                   local_data = np.loadtxt(temp_files[partial_calc])
                   for line in local_data:
                       theta.append(line) 
               else:
                   begin_int = partial_calc
                   print('Found ' + str(len(theta)) + ' values in t = 0')
                   print('Calculations will begin in particle system' + \
                         ' t = 0 part ' + str(begin_int))
                   break 
        else:
            begin_int = 0

        # check if there are left over from previous runs
        if not 'begin_int' in locals():
            raise UnboundLocalError('Erase intermediate files from ' + 
                                    'previous attempts!')

        for iteration in range(begin_int, 
                                int(self.params['split_output'][0])):
            time_ini = time.time()
            args = [self.params for item in 
                    range(self.params['Mini'] // \
                                    int(self.params['split_output'][0]))]

            if self.params['ncores'] > 0:
                pool = Pool(processes=self.params['ncores'])
                p = pool.map_async(SetDistanceFromSimulation, args)
                try:
                     dist = p.get(0xFFFF)
                except KeyboardInterrupt:
                    print('Interruputed by the user!')
                    sys.exit()

                pool.close()
                pool.join()

            else:
                dist = [SetDistanceFromSimulation(obj) for obj in args]

            time_end = time.time() - time_ini

            total_time = time_end/self.M

            #initiate variables to store total number of draws
            #and surviving parameters 
            theta_local = []    

            for line in dist:
                theta_t = [par for par in line[2]]
                theta_t = theta_t + list(line[0])            
                theta_t.append(str(self.params['Mini']/self.M))
                theta_t.append(total_time)
                theta_local.append(theta_t) 
                theta.append(theta_t)           

            if output:
                ftemp = open(self.params['file_root'] + '0_p' + 
                             str(iteration) + '.dat', 'w')   
                for line in theta_local:
                    for element in line:
                        ftemp.write(str(element) + '    ')
                    ftemp.write('\n')
                ftemp.close()    

        #choose smaller distance 
        d1 = np.array([ np.sqrt(sum(line[j]**2 
                         for j in range(len(self.params['param_to_fit']),
                                         len(self.params['param_to_fit']) 
                                         + self.params['dist_dim']))) 
                                         for line in theta])

        d1B = list(d1) 
        d1.sort()

        indx = [d1B.index(item) for item in d1]

        theta_new = [theta[elem] for elem in indx[:self.M]]    

        #write particle system to file
        op = open(self.params['file_root'] + '0.dat', 'w')
        for item in self.params['param_to_fit']:
            op.write(item  + '    ')

        for i2 in range(self.params['dist_dim']):
            op.write('distance' + str(i2 + 1) + '    ')     
 
        op.write('NDraws    time       ')
        for i2 in range(self.params['dist_dim']):
            op.write('dist_threshold' + str(i2 + 1) + '    ')
        op.write('\n')

        for line in theta_new:
            for elem in line:
                op.write(str(elem)  + '    ')
            for i3 in range(self.params['dist_dim']):
                indx = self.params['dist_dim'] + i3
                op.write(str(max(np.array(theta_new)[:,-indx])) + '    ')
            op.write('\n')
        op.close()

        #determine initial weights
        W1 = [1.0/self.M for i2 in range(self.M)]

        op2 = open(self.params['file_root'] + '0weights.dat', 'w')
        for item in W1:
            op2.write(str(item) + '\n')
        op2.close()

        #erase temporary files 
        if output:
            for iteration in range(int(self.params['split_output'][0])):
                os.remove(self.params['file_root'] + '0_p' + 
                          str(iteration) + '.dat')
 
        return np.array(theta_new), self.params['Mini']

    def BuildPSystem(self, previous_particle_system, W, t, output=True):
        """
        Build particle system. 

        input: 	previous_particle_system -> model parameters surviving 
                                            previous distance threshold
	                                    collumns -> [param_to_fit[0], 
                                                         param_to_fit[1], ..,
                                                         param_to_fit[n], 
                                                       distance_from_dataset1] 	

                W -> vector of weights
                t -> particle system index

        output: updated particle system
                collumns -> [surviving_model_parameters, distance, 
                             number_necessary_draws, computational_time (s), 
                             distance_threshold]

        """
        
        #calculate weighted covariance matrix from previous particle system
        indx = len(self.params['param_to_fit'])
        ds = DescrStatsW(previous_particle_system[:,:indx], weights=W)
        cov1 = ds.cov
    
        var = {}
        var['params'] = self.params
        var['W'] = W
        var['previous_particle_system'] = previous_particle_system
        var['previous_cov_matrix'] = cov1
        surv_param = []

        #check if there are previous partial results
        temp_files = [self.params['file_root'] + str(t) + '_p' + 
                                                 str(part) + '.dat'
                      for part in range(int(self.params['split_output'][0]))]

        file_list = os.listdir(os.getcwd())
        if temp_files[0] in file_list:
            for partial_calc in range(int(self.params['split_output'][0])):
                if temp_files[partial_calc] in file_list:
                    local_data = np.loadtxt(temp_files[partial_calc])
                    for line in local_data:
                        surv_param.append(line) 
                else:
                    begin_int = partial_calc
                    print ('Found ' + str(len(surv_param)) + \
                          ' values in t = ' + str(t))
                    print ('Calculations will begin in particle system t = ' + \
                          str(t) + ' part ' + str(begin_int))
                    break 
        else:
            begin_int = 0
            
        #run sampler in separate chuncks 
        for iteration in range(begin_int, int(self.params['split_output'][0])):   
            args = [var for j in range(self.params['M']//int(self.params['split_output'][0]))]

            if self.params['ncores'] > 0:
                pool = Pool(self.params['ncores'])
                p = pool.map_async(SelectParamInnerLoop, args)
                try:
                    surv_param_local = p.get(0xFFFF)
                except KeyboardInterrupt:
                    print('Interruputed by the user!')
                    sys.exit()

                pool.close()
                pool.join() 

            else:
                surv_param_local = [SelectParamInnerLoop(item) for item in args]

            for item in surv_param_local:
                surv_param.append(item)

            if output:
                ftemp = open(self.params['file_root'] + str(t) + '_p' + 
                             str(iteration) + '.dat', 'w')  
                for line in surv_param_local:
                    for element in line:
                        ftemp.write(str(element) + '    ')
                    ftemp.write('\n')
                ftemp.close() 

        #begin writing output file
        op = open(self.params['file_root'] + str(t) + '.dat' , 'w')
        for item in self.params['param_to_fit']:
            op.write(item  + '    ' )
        for jj in range(self.params['dist_dim']):
            op.write('distance' + str(jj) + '    ')     
        op.write('NDraws    time    ')
        for kk in range(self.params['dist_dim']): 
            op.write('dist_threshold' + str(kk + 1) + '    ')
        op.write('\n')
        for line in surv_param:
            for elem in line:
                op.write(str(elem) + '    ')
            op.write('\n')
        op.close()

        if output:
            for iteration in range(int(self.params['split_output'][0])):   
                os.remove(self.params['file_root'] + str(t) + '_p' + 
                          str(iteration) + '.dat')

        return np.array(surv_param)


    def UpdateWeights(self, W, previous_particle_system, 
                            current_particle_system, output=True):   
        """
        Update weights given new particle system.

        input: 	W ->  vector of current weights

                previous_particle_system -> model parameters surviving 
                                            previous distance threshold
	       	                            collumns -> [param_to_fit[0], 
                                                         param_to_fit[1], .., 
                                                         param_to_fit[n], 
                                                       distance_from_dataset1] 	

                current_particle_system -> model parameters surviving 
                                           current distance threshold
	       	                           collumns -> [param_to_fit[0], 
                                                        param_to_fit[1], .., 
                                                        param_to_fit[n], 
                                                       distance_from_dataset1] 	

        output:   vector of updated weights                       
        """

        print('updating weights')

        #calculate weighted covariance matrix from previous particle system
        indx = len(self.params['param_to_fit'])
        ds = DescrStatsW(previous_particle_system[:,:indx], weights=W)
        cov1 = ds.cov

        new_weights = []

        #determine prior distributions
        distributions = [self.prior[par]['func'](self.params['prior'][par], 
                                                 func=True) 
                                       for par in self.params['param_to_fit']]

        for i4 in range(len(current_particle_system)):
            nominator = np.prod([distributions[i2].pdf(
                                   current_particle_system[i4][i2]
                                                      )  
                                   for i2 in range(
                                      len(self.params['param_to_fit'])
                                                   )
                                 ])

            indx = len(self.params['param_to_fit'])
            denominator = sum(W[i3]*multivariate_normal.pdf(
                              current_particle_system[i4][:indx], 
                              previous_particle_system[i3][:indx], cov=cov1) 
                              for i3 in range(len(W)))

            new_weights.append(nominator/denominator)

        final_weights = [item/sum(new_weights) for item in new_weights]

        if output == True:
            op = open(self.params['file_root'] + str(self.T) + 
                      'weights.dat', 'w')
            for item in final_weights:
                op.write(str(item) + '\n')
            op.close()

        return final_weights
              


    def fullABC(self, build_first_system=False, nruns=-9):
        """
        Run complete ABC sampler algorithm. 

        input:	root_file_name -> root of file name to be used 
                                  in all runs (string)

		build_first_system (optional) -> boolean (read or generate 
                                                 first particle system). 
                                                 Default is False.

                nruns (optional) -> int (number of ABC iterations to run)
                                    if < 0 use convergence criteria
                                    elif > 0 run only through this number
                                             of iterations
                                    default is -9
   	
        output:	particle systems and corresponding 
                weights written in data files.
        """

        #determine initial weights
        W1 = [1.0/self.M for i2 in range(self.M)]

        if build_first_system == True:
            #build first particle system
            sys0 = self.BuildFirstPSystem()

        #read first particle system from file
        op = open(self.params['file_root'] + '0.dat', 'r')
        lin = op.readlines()
        op.close()

        t1 = np.atleast_2d([elem.split() for elem in lin[1:]])

        sys1 = np.array([np.array([float(line[i1]) 
                         for i1 in range(len(self.params['param_to_fit']) 
                                        + self.params['dist_dim'])]) 
                              for line in t1])

        #determine number of draws in previous particle system generation
        K =  sum(int(float(line[len(self.params['param_to_fit']) + 
                 self.params['dist_dim']])) for line in t1)

        if self.params['screen']:
            print('number of draws PS0 = ' + str(K))

        #initiate iteration counter
        t = 0
        K = self.M

        if nruns < 0:
            while float(self.M)/K > self.delta:      

                t = t + 1

                self.T = t 

                sys_new = self.BuildPSystem(sys1, W1, t)
        
                W2 = self.UpdateWeights(W1, sys1, sys_new)

                K = sum(sys_new[:, (len(self.params['param_to_fit']) + 
                        self.params['dist_dim'])])

                del sys1, W1

                sys1 = sys_new
                W1 = W2

                del sys_new, W2 

                print(' finished PS ' + str(t) + ',    convergence = ' + \
                       str(float(self.M)/K))
           
            self.T = t

        elif nruns > 0:
            for iterations in range(nruns):

                t = t + 1

                self.T = t 

                sys_new = self.BuildPSystem(sys1, W1, t)
        
                W2 = self.UpdateWeights(W1, sys1, sys_new)

                K = sum(sys_new[:, (len(self.params['param_to_fit']) + 
                            self.params['dist_dim'])])

                del sys1, W1

                sys1 = sys_new
                W1 = W2

                del sys_new, W2 

                print(' finished PS ' + str(t) + ',    convergence = ' + 
                       str(float(self.M)/K))
           
                self.T = t
        
    def  ContinueStoppedRun(self, t, nruns=-9):
        """
        Continue ABC sampler algorithm from a specific time-step (run). 

        input: 	t -> index of last completed particle system (int)

                nruns (optional) -> int (number of ABC iterations to run)
                                    if < 0 use convergence criteria
                                    elif > 0 run only through this number
                                             of iterations
                                    default is -9

	output:	subsequent particle systems and corresponding weights 
                written in data files.
        """
 

        op = open(self.params['file_root'] + str(t) + '.dat', 'r')
        lin = op.readlines()
        op.close()

        t1 = [elem.split() for elem in lin]
        
        sys1 = np.array([np.array([float(line[i1]) 
                          for i1 in range(len(self.params['param_to_fit']) + \
                                           self.params['dist_dim'])]) 
                                            for line in t1[1:]])
        

        #determine number of draws in previous particle system generation
        K =  int(sum(float(line[list(t1[0]).index('NDraws')]) 
                     for line in t1[1:]))
        print('number of draws PS' + str(t) + ' = ' + str(K))
        
        if t > 0:        
            W1 = np.loadtxt(self.params['file_root'] + str(t) + 'weights.dat')
        elif t == 0:
            W1 = [1.0/self.M for i2 in range(self.M)]
    

        if nruns < 0:
            while float(self.M)/K > self.delta:

                t = t + 1

                self.T = t

                sys_new = self.BuildPSystem(sys1, W1, t)
        
                W2 = self.UpdateWeights(W1, sys1, sys_new)

 
                K = sum(sys_new[:, len(self.params['param_to_fit' ]) + 
                                       self.params['dist_dim']])

                del sys1, W1

                sys1 = sys_new
                W1 = W2

                del sys_new, W2 


                print(' finished PS' + str(t) + ',    convergence = ' + \
                       str(float(self.M)/K))
        
                  
            self.T = t

        elif nruns > 0:
            for iterations in range(nruns):

                t = t + 1

                self.T = t

                sys_new = self.BuildPSystem(sys1, W1, t)
        
                W2 = self.UpdateWeights(W1, sys1, sys_new)

 
                K = sum(sys_new[:, len(self.params['param_to_fit' ]) + 
                                       self.params['dist_dim']])

                del sys1, W1

                sys1 = sys_new
                W1 = W2

                del sys_new, W2 


                print(' finished PS' + str(t) + ',    convergence = ' + \
                       str(float(self.M)/K))
        
                  
            self.T = t




def main():
  print(__doc__)



if __name__=='__main__':
  main()    

    




