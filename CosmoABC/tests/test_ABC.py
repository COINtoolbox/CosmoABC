from unittest import TestCase

from CosmoABC.distances import *
from CosmoABC.priors import *
from CosmoABC.ABC_sampler import *
from CosmoABC.plots import *
import numpy
from statsmodels.stats.weightstats import DescrStatsW

def ysim( v ):

    l1 = numpy.random.normal( loc=v['mu'], scale=0.1, size=100 )
    
    return numpy.atleast_2d( l1 ).T 


mu = 2.5
v1 = {'mu': mu }
data = ysim( v1 )

params = {}
params['param_to_fit']=['mu' ]							
params['param_lim']=[[2.0, 3.0]]	
params['prior_par'] = [[2.0, 3.0]]
params['simulation_params'] = v1


params['s']=0.15					
params['epsilon1'] = 0.5			
params['M'] = 10				
params['delta'] =1.0				
params['qthreshold'] = 0.75			

params['file_root'] = 'example_PS'			#root to output file name for subsequent particle systems

success = 0

#initiate ABC sampler
sampler_ABC = ABC( dataset1=data, params=params, simulation_func=ysim, prior_func=[ flat_prior ], distance_func=distance_GRBF) 

success = success + 1

print 'Passed test  ' + str( success ) 

#draw parameters
r1 = sampler_ABC.DrawAllParams()

success = success + 1

print 'Passed test  ' + str( success ) 

#set distance
r2 = sampler_ABC.SetDistanceFromSimulation()

success = success + 1

print 'Passed test  ' + str( success )

#Build First Particle system
r3 = sampler_ABC.BuildFirstPSystem( output=False )

success = success + 1

print 'Passed test  ' + str( success )

#determine weights
W = [ 1.0/params['M'] for i in xrange( params['M'] ) ]

r3B=numpy.array([ [float( line[0] ), float( line[1] ), float( line[2]) ] for line in r3[0] ])

#determine previous covariance matrix
ds = DescrStatsW( r3B[:,:len(params['param_to_fit'])], weights=W )
cov1 = ds.cov

#Select Parameters for subsequent loop
r4 = sampler_ABC.SelectParamInnerLoop( r3B, W, cov1 )

success = success + 1

print 'Passed test  ' + str( success )

#build second particle system
r5 = sampler_ABC.BuildPSystem( r3B, W, 1, filename=None  )

success = success + 1

print 'Passed test  ' + str( success )

#update weights
r6 = sampler_ABC.UpdateWeights( W, r3B, r5, output=False )

success = success + 1

print 'Passed test  ' + str( success )

print 'All main ABC functions survived initial tests!!!'
