#!/usr/bin/env python

"""
Approximate Bayesian Computation library.
"""

__author__ = "E. E. O. Ishida, S. D. P. Vitenti, M. Penna-Lima"
__maintainer__ = "E. E. O. Ishida"
__copyright__ = "Copyright 2015"
__version__ = "0.1"
__email__ = "emilleishida@gmail.com"
__status__ = "Prototype"
__license__ = "GPL"



from scipy.stats import multivariate_normal
import numpy
from scipy.interpolate import Rbf
import time
from scipy.stats.mstats import mquantiles
from scipy.stats import norm
from statsmodels.stats.weightstats import DescrStatsW
import sys


################################################


class ABC( object ):

    def __init__( self, dataset1=False, params=False, simulation_func=False, distance_func=False, prior_func=False ):

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
	s		-	Scale parameter
	M		-	Number of particles in each particle system.
        epsilon1	-	Distance threshold for first particle system.
	qthreshold	-	Quantile for choosing subsequent distance thresholds.
        T		-	Total number of particle system at convergence
	

	Method attributes
	------------------------

	pars		-	Set of parameters drawn from prior distribution.
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

        self.data	 	= dataset1    			#set of 2-dimensional arrays of "real" catalog.         
        self.simulation		= simulation_func		#function which performs the simulation
        self.params		= params			#all parameters 
        self.distance		= distance_func                 #distance function
        self.prior		= prior_func                    #prior distribution function

        #check length of observed data
        if isinstance( self.data, bool):
            raise IOError( 'No real data catalog provided.' )

        if len( self.data ) < 2 :
            raise IOError( 'Real data catalog too short. Table must contain more than 1 element.' )

        try:
            len( self.data[0] )
        except TypeError:
            raise TypeError( 'Real data catalog must be at least 2 dimensional.' )
                

        #check minimum keywords in params
        self.min_keys = ['simulation_params', 'param_to_fit', 'prior_par', 'param_lim', 'M', 'epsilon1', 'qthreshold', 'delta','s', 'file_root']  
        for item in self.min_keys:
            if item not in self.params.keys():
                raise IOError( 'Keyword ' + str( item ) + '  is missing from inputed dictionary of parameters (params)!' ) 

        #check simulation function
        if not self.simulation:
            raise IOError( 'Please, provide a valid simulation function. \n See file "simulation.py". ' )

        #check prior function
        if not self.prior:
            raise IOError( 'Please provide a valid prior distribution function. \n "See file priors.py"' )

        #check distance function
        if not self.distance:
            raise IOError( 'Please provide a valid distance function. \n See file "distances.py"' )


           
        #define parameters separately so inputed values are not changed       
        self.delta		= params['delta']		#convergence criteria 
        self.s 		        = params['s']			#scale parameter (scalar)
        self.M		        = params['M']			#number of elements in each particle system
        self.epsilon1    	= params['epsilon1']		#distance threshold for first particle system
        self.qthreshold		= params['qthreshold']		#quantile to define the distance threshold for subsequent particle system
    


    def DrawAllParams( self ):
        """
        Draw complete set of  parameters from prior.

        output: array of parameters sampled from the prior
        """

        pars = []
        for j in range( len( self.params[ 'param_to_fit' ] ) ):

            
            p1 = self.prior(  self.params[ 'prior_par' ][ j ], self.params[ 'param_lim' ][ j ] )   

            pars.append( p1 )

        return numpy.array( pars )

   
    def SetDistanceFromSimulation( self ):
        """
        Draw cosmological parameter values from prior, generate simulation and calculate distance from a given comparison  catalog. 
        In the context of Ishida et al., 2015 that translates into calculating the distance from the real to the simulated cataloges.

        output: scalar (distance between dataset1 and dataset2)
    
        """

        #draw parameters from prior
        ParTry = self.DrawAllParams()

        #update cosmological parameter values in dictionary
        for i1 in range( len( self.params[ 'param_to_fit' ] ) ):
            self.params[ self.params[ 'param_to_fit' ][ i1 ] ] = ParTry[ i1 ]

        #generate simulation
        DataSimul = self.simulation( self.params['simulation_params'] )
             


        #calculate distance
        dist = self.distance( self.data, DataSimul, self.s )

        if dist > 0:
            return dist
        else:
            print 'dist = ' + str(dist)
            print 'DataSimul = ' + str( DataSimul )
            print 'ParTry = '+str(ParTry)
            op1=open('simulation.dat','w')
            for line in DataSimul:
                for item in line:
                    op1.write( str(item) + '    ')
                op1.write( '\n' )
            op1.close() 
            sys.exit()
    

    def BuildFirstPSystem( self, output=True, filename='PS_0.dat' ):
        """
        Build the first particle system, storing the parameter values satisfying distance threshold. 

        input:	output (optional) -> boolean (choose to write output data file, default is True)
    		filename (optional) -> string (name for output file, default is 'PS_0.dat')
                 

        output:	[0] array containing first particle system
		    collumns -> [ param_to_fit[0], param_to_fit[1], .., param_to_fit[n], distance_from_dataset1 ] 	

		[1] integer -> total number of draws necessary to build the first particle system
        """

        print 'Building first particle system:'

        #initiate variables to store total number of draws and surviving parameters
      
        theta = []    
   
        while len( theta ) < self.M:
 
            dist = 10**10
            time0 = time.time()
            K = 0

            while dist > self.epsilon1: 
                K = K + 1 

                time1 = time.time()
                dist = self.SetDistanceFromSimulation()
 
            theta_t = [ self.params[ item ] for item in self.params['param_to_fit'] ] 

            if dist > 0:

                total_time = time.time() - time1

                theta_t.append( dist )
                theta_t.append( str( K ) )
                theta_t.append( total_time )
                theta.append( theta_t  )
        
                print '        particle index = ' + str( len( theta) ) + ',  number of draws = ' + str( K ) + ',   distance=' + str( dist )

        

        #write particle system to file
        if output:
            op = open( filename, 'w' )
            for item in self.params[ 'param_to_fit' ]:
                op.write( item  + '    '  )
            op.write(  'distance     NDraws    time        dist_threshold\n' )
            for line in theta:
                for elem in line:
                    op.write( str( elem )  + '    ' )
                op.write( str( self.epsilon1 ) + '\n' )
            op.close()
 
        return numpy.array( theta ), K


    def SelectParamInnerLoop( self, previous_particle_system, W, previous_cov_matrix ):
        """
        Draw model parameters based on previous particle system and return those satisfying distance threshold.

        input:	previous_particle_system -> model parameters surviving previous distance threshold
		collumns -> [ param_to_fit[0], param_to_fit[1], .., param_to_fit[n], distance_from_dataset1 ] 	

		W -> vector of weights 

		previous_cov_matrix -> covariance matrix based on previous particle system results
		cosmological_parameters -> dictionary of necessary cosmological parameters
		keys must include: [H0, Omegab, Omegam, Tgamma0, ns, sigma8, w]	

		epsilon -> distance threshold to be satisfied	

        output: vector -> [ surviving_model_parameters, distance, number_necessary_draws, computational_time (s), distance_threshold ]
        """

        dist = 10**10
        K = 0

        #time marker
        time_start = time.time()

        #determine distance threshold
        epsilon = mquantiles( previous_particle_system[:,len(self.params['param_to_fit'])], prob=self.qthreshold )[0]

        while dist > epsilon:
 
            #update counter
            K = K + 1 
  
            #draw model parameters to serve as mean 
            index_theta0 = numpy.random.choice( xrange( len( W ) ), p=W )
            theta0 = numpy.atleast_2d( previous_particle_system[ index_theta0 ][: len( self.params['param_to_fit']) ] )
  

            #initialize boolean parmeter vector
            theta_t = [ False for i in xrange( len( self.params['param_to_fit']) ) ]
            
            cont = 0
            while False in theta_t:

                cont = cont + 1
          
                #draw model parameter values for simulation
                mvn = multivariate_normal.rvs( theta0, previous_cov_matrix )

                try:
                    len( mvn )
                    theta_t_try = list( mvn )
                except TypeError:
                    theta_t_try = [ mvn ]

                theta_t = []
                
                for k1 in xrange( len( self.params['param_to_fit'] ) ):

                    if  theta_t_try[ k1 ] >= self.params[ 'param_lim' ][ k1 ][0] and theta_t_try[ k1 ] < self.params[ 'param_lim' ][ k1 ][1]:
                        theta_t.append( True )

                    else:
                        theta_t.append( False )          

            #update parameter values in dictionary
            for i1 in range( len( self.params[ 'param_to_fit' ] ) ):
                self.params[ self.params[ 'param_to_fit' ][ i1 ] ] = theta_t_try[ i1 ]

            #generate simulation
            DataSimul = self.simulation( self.params['simulation_params'] )
      
            #calculate distance
            dist = self.distance( self.data, DataSimul, self.s )
 

        theta_t_try.append( dist )
        theta_t_try.append( K )    
        theta_t_try.append( time.time() - time_start )
        theta_t_try.append( epsilon )

        return theta_t_try


    def BuildPSystem( self, previous_particle_system, W, t, filename='system.dat'  ):
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
        ds = DescrStatsW( previous_particle_system[:,:len(self.params['param_to_fit'])], weights=W )
        cov1 = ds.cov


        if filename != None:
            #begin writing output file
            op = open( filename, 'w' )
            for item in self.params[ 'param_to_fit' ]:
                op.write( item  + '    '  )
            op.write(  'distance     NDraws    time    dist_threshold\n' )

        particle_system = []

        for j in xrange( self.M ):

            surv_param = self.SelectParamInnerLoop( previous_particle_system,  W, cov1 )

            particle_system.append( surv_param )

            print 'J = ' + str( j ) + ',    K = ' + str( surv_param[ len( self.params['param_to_fit' ]) + 1  ] ) + ',    distance=' + str( surv_param[len( self.params['param_to_fit' ])  ] )
            
            if str( surv_param[len( self.params['param_to_fit' ])  ] ) == 'nan':
                raise ValueError( '"nan" is not an acceptable distance value!' )

            if filename != None:
                for elem in surv_param:
                    op.write( str( elem ) + '    ' )
                op.write( '\n' )

        if filename != None:
            op.close()

        return numpy.array( particle_system )


    def UpdateWeights( self, W, previous_particle_system, current_particle_system, output=True, filename='weights.dat' ):   
        """
        Update weights given new particle system.

        input: 	W ->  vector of current weights

                previous_particle_system -> model parameters surviving previous distance threshold
	       	collumns -> [ param_to_fit[0], param_to_fit[1], .., param_to_fit[n], distance_from_dataset1 ] 	

                current_particle_system -> model parameters surviving current distance threshold
	       	collumns -> [ param_to_fit[0], param_to_fit[1], .., param_to_fit[n], distance_from_dataset1 ] 	

        output:   vector of updated weights      
                 
        """


        #calculate weighted covariance matrix from previous particle system
        ds = DescrStatsW( previous_particle_system[:,:len(self.params['param_to_fit'])], weights=W )
        cov1 = ds.cov

        new_weights = []

        #determine prior distributions
        distributions = [ self.prior( self.params['prior_par'][ i1 ], self.params['param_lim'][ i1 ], func=True) for i1 in xrange( len( self.params[ 'param_to_fit' ] ) ) ]

        for i4 in range( len( current_particle_system ) ):
              
            nominator = numpy.prod( [ distributions[ i2 ].pdf( current_particle_system[ i4 ][ i2 ] ) for i2 in xrange( len(  self.params[ 'param_to_fit' ] ) ) ] )

            denominator = sum( W[ i3 ]*multivariate_normal.pdf( current_particle_system[ i4 ][:len(self.params['param_to_fit'])], previous_particle_system[ i3 ][:len(self.params['param_to_fit'])], cov=cov1 ) for i3 in xrange( len( W ) ) )

            new_weights.append( nominator/denominator )

        final_weights = [ item/sum( new_weights ) for item in new_weights ]

        if output == True:
            op = open( filename, 'w')
            for item in final_weights:
                op.write( str( item ) + '\n' )
            op.close()

        return final_weights
              


    def fullABC( self, root_file_name, build_first_system=False ):
        """
        Run complete ABC sampler algorithm. 

        input:	root_file_name -> root of file name to be used in all runs (string)

		build_first_system (optional) -> boolean  (read or generate first particle system). Default is False.
   	
        output:	particle systems and corresponding weights written in data files.
        """
 

        if build_first_system == True:

            #build first particle system
            sys0 = self.BuildFirstPSystem( filename=root_file_name+'0.dat' )

            op = open( root_file_name+'0.dat', 'r' )
            lin = op.readlines()
            op.close()

            t1 = [ elem.split() for elem in lin[1:] ]

            sys1 = numpy.array([ numpy.array([ float( line[ i1 ] ) for i1 in xrange( len( self.params['param_to_fit' ] ) + 1 ) ]) for line in t1 ])

            
        else:
            #read first particle system from file
            op = open( root_file_name+'0.dat', 'r' )
            lin = op.readlines()
            op.close()

            t1 = [ elem.split() for elem in lin[1:] ]

            sys1 = numpy.array([ numpy.array([ float( line[ i1 ] ) for i1 in xrange( len( self.params['param_to_fit' ] ) + 1 ) ]) for line in t1 ])

        #determine number of draws in previous particle system generation
        K =  sum( int( line[ len( self.params['param_to_fit' ] ) + 1  ] ) for line in t1 )
        print 'K = ' + str( K )

        #determine initial weights
        W1 = [ 1.0/self.M for i2 in xrange( self.M ) ]

        #initiate iteration counter
        t = 0

        while float( self.M )/K > self.delta:

            t = t + 1

            sys_new = self.BuildPSystem( sys1, W1, t, filename=root_file_name + str( t ) + '.dat'  )
        
            W2 = self.UpdateWeights( W1, sys1, sys_new, filename=root_file_name + str( t ) + 'weights.dat' )

 
            K = sum( sys_new[:, len( self.params['param_to_fit' ] ) + 1 ] )

            del sys1, W1

            sys1 = sys_new
            W1 = W2

            del sys_new, W2 

            print ' T = ' + str( t ) + ',    convergence = ' + str( float( self.M )/K )

        self.T = t
        

        
    def  ContinueStoppedRun( self, t, root_file_name ):
        """
        Continue ABC sampler algorithm from a specific time-step (run). 

        input: 	t -> index of last completed particle system (int)

		root_file_name -> root of file name to be used in all subsequent runs (string)
	
	output:	subsequent particle systems and corresponding weights written in data files.
        """
 

        op = open( root_file_name+str(t)+'.dat', 'r' )
        lin = op.readlines()
        op.close()

        t1 = [ elem.split() for elem in lin[1:] ]
        
        sys1 = numpy.array([ numpy.array([ float( line[ i1 ] ) for i1 in xrange( len( self.params['param_to_fit' ] ) + 1 ) ]) for line in t1 ])
        
        #determine number of draws in previous particle system generation
        K =  sum( int( line[ len( self.params['param_to_fit' ] ) + 1 ] ) for line in t1 )
        print 'K = ' + str( K )
        
                
        W1 = numpy.loadtxt( root_file_name+str(t)+'weights.dat' )

        while float( self.M )/K > self.delta:

            t = t + 1

            sys_new = self.BuildPSystem( sys1, W1, t, filename=root_file_name + '_' + str( t ) + '.dat'  )
        
            W2 = self.UpdateWeights( W1, sys1, sys_new, filename=root_file_name + '_' + str( t ) + 'weights.dat' )

 
            K = sum( sys_new[:, len( self.params['param_to_fit' ] ) + 1 ] )

            del sys1, W1

            sys1 = sys_new
            W1 = W2

            del sys_new, W2 

            print ' T = ' + str( t ) + ',    convergence = ' + str( float( self.M )/K )
        
                  
        self.T = t



def main():
  print(__doc__)



if __name__=='__main__':
  main()    

    




