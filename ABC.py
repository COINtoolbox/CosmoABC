#!/usr/bin/python 

import sys
import numpy as np
import math


from scipy import stats
from scipy.stats.mstats import mquantiles
import multiprocessing as mp
#from pathos.pp import ParallelPythonPool as Pool
#import  pathos.multiprocessing as mp
import time as time




class StoreInfo( object ):
	
	#def __init__(self,n_iter, Ncpu, N, prior,dist):
	def __init__(self,n_iter, Ncpu, N):
	#def __init__(self):
		
		#self.prior=prior
		
		#elf.kernel=Kernel
		
		#self.rho=dist
		
		
		############################
		# Variable to store data 
		self.n_iter= n_iter 
		self.Ncpu = Ncpu
		self.N = N
		self.index = None
		self.sdata_old = None
		self.sdata = None
		self.epsilon = None
		self.mean = None
		self.weighted_mean = None
		self.sdata_c = None
		self.median = None
		self.weighted_cov = None
		self.cov = None
		self.weights = None
		

class Distributions( object ):
	
	def __init__(self,prior = "Normal" , mean = None, cov = None, bounds = None):
		
		self.prior=prior
		
		self.mean = mean
		
		self.cov = cov
		
		self.size = len(mean)
		
		self.bounds = bounds
		
		
		if ( self.prior == "Uniform" ):
			
			self.rv=stats.uniform( self.bounds[0], self.bounds[1] )
		
		elif ( self.prior == "Normal" ):
			
			self.rv=stats.multivariate_normal( self.mean , self.cov  )
			
		else:
			sys.exit("Stop: Invalid Prior name\n")
			
			
		
		
	
	def gen(self):
		
		if (self.prior == "Normal" ):
			
			flag = np.array([ False for i in xrange( self.size ) ])
			
			while (flag.all() == False ): 
				
				p = self.rv.rvs()
				
				flag=((p>=self.bounds[0]) & (p<=self.bounds[1]))
				
		else:
			
			p=self.rv.rvs()
		
		return p
	
	
	def pdf( self, x ):
		
		return self.rv.pdf(x)
	
	def set_cov( self , cov ):
		
		self.cov = cov
		
	def set_mean( self , mean):
		
		self.mean = mean



class RHO( object ):
	
	
	#def __init__( self , data_obs ,  model , SummaryStatistics , dist, *args):
	def __init__( self , data_obs ,  model , SummaryStatistics , dist, *args):
		
		self.model = model 
		
		self.sumstats = SummaryStatistics
				
		self.data_obs = data_obs
		
		self.sum_stats_obs = self.sumstats ( data_obs )
		
		#print self.sum_stats_obs
		
		self.dist = dist
		
		self.args=args
		
		#self.Xobs = data_obs[0]
		
		#self.Yobs = data_obs[1]
		
		#self.Ebos = data_obs[2]
	
	#def __call__(self,params):
	def call(self,params):
		
		simul_data=self.model ( params  ,self.args)
		#print simul_data
		
		sum_stats_model = self.sumstats ( simul_data )
		
		#print sum_stats_model 
		
		#distance =  self.diff( self.sum_stats_obs , sum_stats_model ) 
		return self.dist ( self.sum_stats_obs , sum_stats_model )
		#return  np.sum((self.sum_stats_obs  - sum_stats_model)**2)
		
		


		
def kernel(  mean, cov, bounds ):
     
    
    
    n=len( mean )
    
   # U,S,V=linalg.svd(hyper_par_cov)
    
   # new_mean=numpy.dot(numpy.transpose(U),hyper_par)
   
    #print mean
    
    #print cov
    
    #print bounds
    

    flag = np.array([ False for i in xrange( n ) ])
    
    #print "***********************************"
    
    while (flag.all() == False ):
		
          L = np.linalg.cholesky( cov )
          
          #print "L:", L
          
          x = np.random.normal( size=n )
          #print x
          
          params =np.array( mean ) + np.dot(L, x )
          
          #print "candidate:", params
          
          flag=((params>=bounds[0]) & (params<=bounds[1]))
          
         # print params,flag
    #print "***********************************"
    #print "generated parameters=", numpy.array( params[0] )
    return np.array(params)
    

#def first_inner_loop(storeinfo,prior,rho):
def first_inner_loop(args):
	
		storeinfo,prior,rho = args
		
		#print storeinfo.Ncpu
		
		#print storeinfo.N
		
		if ( storeinfo.Ncpu!=1 ): 
			
			size = 2*int(storeinfo.N/storeinfo.Ncpu)
			
		else:
			
			size=2*storeinfo.N
		
		#print size
		
		np.random.seed()
		
		results = []
		
		for i in xrange( size ):
						
			#p_proposed = storeinfo.prior.gen( )
			
			p_proposed = prior.gen( )
			
			#d=storeinfo.rho.call ( p_proposed )
			
			#d=storeinfo.rho.call ( p_proposed )
			
			d=rho.call ( p_proposed )

                        d2 = sum( abs( d[1] - d[2] ) )
			
			results.append(  list( p_proposed ) + [ d[0], d2 ] )
			
			
		
		return results # returns a list
		
def inner_loop(args):
		
		#storeinfo , tol = args
		
		storeinfo,prior,rho , tol = args
		
		np.random.seed()		
		
		index_list=[ ]
		
		results = [ ]
		
		d_tol= tol

                #print 'd_tol = ', d_tol
		
		if ( storeinfo.Ncpu!=1 ): 
			
			size = storeinfo.N/storeinfo.Ncpu
			
		else:
			
			size=storeinfo.N
					
		
		k=0
		
		while (k<size):
		
		
			
			index = weighted_values(range( storeinfo.N ), storeinfo.weights,1) 
			
			#print "Sampled index", index
			
			#print "prior bounds"
			
			#print storeinfo.prior.bounds
			
			#print prior.bounds
			
			#print "choosen triplet=" , storeinfo.sdata[index]
			
			
			#p_proposed = kernel( storeinfo.sdata[index] , storeinfo.weighted_cov, storeinfo.prior.bounds )
			
			p_proposed = kernel( storeinfo.sdata[index] , storeinfo.weighted_cov, prior.bounds )
			
			#print "proposed point:", p_proposed
			
			#d = storeinfo.rho( p_proposed )
			
			d=rho.call ( p_proposed )
			
			#print 'd = ', d
			
			if ( d[0] <= d_tol[0] ):

                                d2 = sum( abs( d[1] - d[2] ) )                

                                #print 'd2 = ', d2

                                if ( d2 <= d_tol[1] ):
				
					print "Particle=%d" % k
				
					k+=1
				
					index_list.append( index )

                                        #print 'index = ', index  
				
					print " Dist1 =%.4f \t Dist2=%.4f \t  epsilon1=%.4f  \t epsilon2 = %.4f"  % (  d[0] , d2, d_tol[0], d_tol[1] )
				
					print "Accepted Point: " + str (p_proposed)  + '\n'
				
					results.append( list( p_proposed ) + [ d[0], d2 ] )
		
		return results , index_list # both objects are lists



		
def DataToFiles(name,data):
	
		
		
		f=open(name,"w")
		
		for i in data:
			#f.write(" ".join([ "%.4\t" for j in i ]  ) + "\n" )
			
			f.write(" ".join([ "%.5f\t" % j for j in i ]  ) + "\n")
		
		f.close()


def weighted_values( values, probabilities, size):
    bins = np.cumsum( probabilities)
    return values[ np.digitize(np.random.random_sample(size), bins) ]



def norm_pdf_multivariate(x, mu, sigma):
  size = len(x)
  if size == len(mu) and (size, size) == sigma.shape:
    det = np.linalg.det(sigma)
    if det == 0:
        raise NameError("The covariance matrix can't be singular")

    norm_const = 1.0/ ( math.pow((2*np.pi),float(size)/2) * math.pow(det,1.0/2) )
    x_mu = np.matrix(x - mu)
    inv = np.matrix(sigma).I        
    result = math.pow(math.e, -0.5 * (x_mu * inv * x_mu.T))
    return norm_const * result
  else:
    raise NameError("The dimensions of the input don't match")
    
 



class ABC(object):
	
	#def __init__(self,SI):
	def __init__(self,SI,Prior,Rho):
		
		self.SI=SI
		self.prior=Prior
		self.rho=Rho
		self.results = None
		
	
	def sampler(self, db ):
		
		q=0.750 #0.75
		
		start_time = time.time()
		
		self.SI.epsilon = [ ]
		
		if (self.SI.Ncpu!=1): 
			
			p = mp.Pool ( self.SI.Ncpu )
		
		for i in xrange ( self.SI.n_iter ):
			
			if ( i == 0 ):
				
				print "First iteration loop!"
				
				if ( self.SI.Ncpu == 1 ):
					
					print "%d CPU being used!" % self.SI.Ncpu
					
					#self.results = first_inner_loop ( self.SI )
					
					self.results = first_inner_loop ( (self.SI, self.prior, self.rho) )
				
					self.SI.sdata = np.array ( self.results )
			
				else:
					
					print "%d CPU being used!" % self.SI.Ncpu
				
					#self.results = p.map ( first_inner_loop , [ self.SI for j in xrange ( self.SI.Ncpu ) ] )
					
					self.results = p.map ( first_inner_loop , [ (self.SI, self.prior, self.rho ) for j in xrange ( self.SI.Ncpu ) ] )
				
					#print self.results
				
					self.SI.sdata = np.vstack ( [ self.results[j] for j in xrange(self.SI.Ncpu) ]  )
				
					#del p
				
				self.SI.weights = np.ones ( self.SI.N ) / self.SI.N	 # set initial weights = 1/Nparticles
			
				indx = np.argsort( self.SI.sdata[:,-2] ) # sorting coputing distances
				
				#print "epsilons:" ,self.SI.sdata[indx,-1]
			
				#print "indx=", indx
			
				self.SI.epsilon.append( [ mquantiles( self.SI.sdata[indx,j1] ,[q] )[0] for j1 in [-2,-1] ] )
				
				print "espilon choosed:", [ mquantiles( self.SI.sdata[indx,j1] ,[q] )[0] for j1 in [-2,-1] ] 
			
			#	print "sdata before ordering"
			
			#	print self.SI.sdata
			
				DataToFiles( db + "_" + str(i) + ".dat", self.SI.sdata[indx])
			
				self.SI.sdata = self.SI.sdata[indx][:self.SI.N,:-2] # sort data according to the computed distances
					
			#	print "sdata after ordering"
			#	print self.SI.sdata
				
				self.SI.mean = np.mean( self.SI.sdata , axis = 0 )
								
				self.SI.weighted_mean =  np.mean( np.transpose( np.transpose (self.SI.sdata) * self.SI.weights ) , axis = 0) # data weighted mean
				
				self.SI.sdata_c = self.SI.sdata[:] - self.SI.weighted_mean[:] # data  rescaled with weighted mean, so it has zero mean 
				
				self.SI.meadian = np.median( self.SI.sdata, axis = 0 )
				
				self.SI.cov = np.cov( self.SI.sdata, rowvar = 0 )
				
				self.SI.weighted_cov = 2*np.dot( np.transpose( self.SI.sdata_c ), np.dot( np.diag( self.SI.weights ) , self.SI.sdata_c ) ) # compute new covariance matrix

				
				print "Data mean"
				
				print self.SI.mean
				
				#print "Data covariance"
					
				#print self.SI.cov
				
				#print "Data weighted covariance"
					
				#print self.SI.weighted_cov
				
				#print "Data median"
				
				#print self.SI.median
				
				#print "Data weighted mean"
				
				#print self.SI.weighted_mean
				
				#print "Data weighted covariance"
				
				del indx
				
			else:
				
				print "Iteration: %d" % i
				
				self.SI.sdata_old = self.SI.sdata [ : ]
				
				if ( self.SI.Ncpu == 1 ):
					
					
					#results =  inner_loop ( (self.SI , self.SI.epsilon[i-1] ) ) # to complete
					
					self.results =  inner_loop ( (self.SI , self.prior, self.rho, self.SI.epsilon[i-1] ) ) # to complete
					
					self.SI.sdata = np.array ( self.results[0] )
					
					self.SI.index = np.array ( self.results[1] )#, dtype=int )
				
				else:
					
					#self.results = p.map (inner_loop , [ (self.SI , self.SI.epsilon[i-1] ) for j in xrange ( self.SI.Ncpu ) ] )
					
					self.results = p.map (inner_loop , [ (self.SI , self.prior, self.rho, self.SI.epsilon[i-1] ) for j in xrange ( self.SI.Ncpu ) ] )
					
					self.SI.sdata = np.vstack ( [ self.results[j][0] for j in xrange(self.SI.Ncpu) ]  )
					
					self.SI.index = np.hstack ( [ self.results[j][1] for j in xrange(self.SI.Ncpu) ]  )
					
					
				DataToFiles( db + "_" + str( i ) + ".dat", self.SI.sdata )
				
				#print "epsilons:", self.SI.sdata[:,-1]
				
				self.SI.epsilon.append( [ mquantiles( self.SI.sdata[:,j1] ,[q] )[0] for j1 in [-2,-1] ] )
				
				print "espilon choosed:", [ mquantiles( self.SI.sdata[:,j1] ,[q] )[0] for j1 in [-2,-1] ] 
				
				self.SI.sdata=self.SI.sdata[:,:-2]
				
				denominator = np.array( [ sum( [self.SI.weights[n1] * norm_pdf_multivariate ( self.SI.sdata[n2], \
				
				self.SI.sdata_old[n1],self.SI.cov ) for n1 in xrange( self.SI.N ) ] ) for n2  in xrange( self.SI.N ) ])
				
				numerator= np.array( [ self.prior.pdf(self.SI.sdata[n1])  for n1 in xrange(self.SI.N) ] )  
				
				self.SI.weights= numerator / denominator
				
				self.SI.weights/= self.SI.weights.sum()
				
				
				
				self.SI.mean = np.mean( self.SI.sdata , axis = 0 )
								
				self.SI.weighted_mean =  np.mean( np.transpose( np.transpose ( self.SI.sdata ) * self.SI.weights ), axis = 0) # data weighted mean
				
				self.SI.sdata_c = self.SI.sdata[:] - self.SI.weighted_mean[:] # data  rescaled with weighted mean, so it has zero mean 
				
				self.SI.meadian = np.median( self.SI.sdata, axis = 0 )
				
				self.SI.cov = np.cov( self.SI.sdata, rowvar = 0 )
				
				self.SI.weighted_cov = 2*np.dot( np.transpose( self.SI.sdata_c ), np.dot( np.diag( self.SI.weights ) , self.SI.sdata_c ) ) # compute new covariance matrix

				
				#print "Data mean"
				
				#print self.SI.mean
				
				#print "Data covariance"
					
				#print self.SI.cov
				
				#print "Data median"
				
				#print self.SI.median
				
				#print "Data weighted mean"
				
				#print self.SI.weighted_mean
				
				#print "Data weighted covariance"
				
				self.SI.weighted_cov
				
				
				
				
				
		
		
		if self.SI.Ncpu!=1:
			p.close()
			p.join()
		else: 
			pass
		print "Done!"
		print "time=%.4f seconds" % (time.time()-start_time)	
			
		#p=mp.Pool(2)
		#p = mp.ProcessingPool(2)
		#result = p.map(first_inner_loop, [ self.x for i in xrange(2)])
		#p.close()
		#p.join()
		
		#return  result
	
		


