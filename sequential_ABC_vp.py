from math import *
from gi.repository import GObject
import matplotlib.pyplot as plt
from gi.repository import NumCosmo as Nc
from gi.repository import NumCosmoMath as Ncm

from scipy.stats.mstats import mquantiles
from scipy.stats import norm,rv_discrete


import numpy
import random 
import os 

from NCountSimul_true_mass_nw_vp import *
from multiprocessing import Pool

import pylab as plt
import time


starttime=time.time()

############################################################
### Fiducial Cosmological parameters: ref. arXiv:1203.5775 (table 5. wCDM CMB+BAO+H0+SNeIa+SPTcl), Tgamma0 and ns are not given.

zmin = 0.3              #minimum redshift
zmax = 1.32             #maximum redshift
H0 = 71.15               #Hubble parameter
Omegam = 0.262           #Dark matter density
Omegab = 0.0439          #Baryon density
Tgamma0 = 2.725          #Radiation temperature today
ns = 0.97                #spectral index 
sigma8 = 0.807           #sigma8
w = -1.01                #Dark energy equation of state     

#mass bin
#dm = [5*10**13, 10**14, 10**14.25, 10**14.5, 10**14.75,  10**15, 10**15.25,  10**15.5, 10**15.75 ]
dm = [10**14.3, 10**14.5, 10**14.7,  10**14.9, 10**15.1, 10**15.3,  10**15.5, 10**15.7 ]

#quantile list
quant_list = [ 0.02, 0.09, 0.25, 0.5, 0.75, 0.91, 0.98]

#sky area
area = 2500

#output directory
path1 = 'all_res_SPT'


#######
#priors over Om (dark matter), w and sigma8


#time_steps = 10

Ncpu=1

time_steps =3 # number of time steps to take

Ninit=100 # initial number of samples 

N=50      # particle sample size after the first iteration

												
epsilon_ini = 1e20		#Starting tolerance

seed = 100  # seed for ncount

##############################################################

path1="teste_3p"



#create output directory
if not os.path.exists( path1 ):
    os.makedirs( path1 )
    
    
    
CosmoParams=ChooseParamsInput()

CosmoParams.params={"H0":H0,"Ob":Omegab,"Om":Omegam,"OL":1.-Omegam-Omegab,"Tgamma":Tgamma0,"ns":ns,"sigma8":0.8,"w":-1.0}

CosmoParams.keys=["Om","w", "sigma8"]

CosmoParams.keys_values=numpy.array( [  0.25 , -1.01,  0.7 ] )

CosmoParams.keys_cov=numpy.diag( [ 0.5 , 0.5, 0.5 ] )**2.

CosmoParams.keys_bounds=numpy.array( [ [0.0 , -3.0, 0.3 ] , [1.0-Omegab, 0.0, 1.0] ] )



#CosmoParams.keys=["Om"]

#CosmoParams.keys_values=numpy.array( [  0.25  ] )

#CosmoParams.keys_cov=numpy.diag( [ 0.5  ] )**2.

#CosmoParams.keys_bounds=numpy.array( [ [0.0  ] , [1.0-Omegab] ] )

CosmoParams.prior_dist="normal"

Nparams=len( CosmoParams.keys )
    
   

epsilon = [ ]

def worker( args ):
    
    return choose_surv_par (*args)

if Ncpu!=1:
  p=Pool(Ncpu)
else:
  pass

for ii in xrange(time_steps):
    print ii
    if (ii==0):
    
        #new simulation object
        ncount = NCountSimul (zmin, zmax, log ( dm[0] ), log ( 10**16 ), area )

        #Generate fiducial data
        data_fid = numpy.array( ncount.simulation( zmax, seed, CosmoParams.params)[1] )

        #Calculate summary statistics for fiducial data
        summ_fid = summary_quantile( data_fid, dm, quant_list )
        
        if ( Ncpu == 1 ):
           
           par_surv, indx,par_cov = choose_surv_par( summ_fid, dm, quant_list, epsilon_ini, Ninit, Ncpu,CosmoParams, zmin, zmax, area, seed )
        
        else: 
           job_args=[ (summ_fid, dm, quant_list, epsilon_ini, N, Ncpu , CosmoParams,  zmin, zmax, area,  seed) for i in xrange(Ncpu)  ]
           
           results = p.map(worker, job_args)
           #p.close()
          # p.join() 
        
           par_surv=[]
           indx=[]
           par_cov=numpy.array(results[0][-1])
           
           for i in xrange(Ncpu):
              
              
              par_surv.append(results[i][0])
           
              indx.append(results[i][1])
           
           par_surv = numpy.vstack( par_surv )
           indx = numpy.vstack( indx ).flatten()
        
       # print "Data generated"
       # print par_surv
       # print "indexes"
       # print indx
       # print "covariance"
       # print par_cov
        
        par_surv=par_surv[ par_surv[:,-1].argsort()  ]
      
        par_surv=par_surv[xrange(N)]
        
        epsilon.append( mquantiles( [ par_surv[ j1 , -1] for j1 in xrange( len( par_surv ) ) ], [0.75] )[0]  )
        
        CosmoParams.sdata=par_surv[:]
        
        CosmoParams.sdata_weights = numpy.ones( N )/N
        
        
      
      #write results to file
        #write results to file
        op0 = open(os.path.join( path1 , 'sequential_ABC_SPT_res_0.dat'), 'w')
        op0.write( '#' ) 
        for par in CosmoParams.keys:
            op0.write('#' + par + '    ' )
        op0.write( 'distance    epsilon \n' )
        for elem in par_surv:
           for item in elem:
              op0.write( str( item ) + '    ' )
           op0.write( str( epsilon_ini ) + '\n' )
        op0.close()

      
    else:
        
        if ( Ncpu == 1 ): 
           
           par_surv, indx, par_cov = choose_surv_par( summ_fid, dm, quant_list, epsilon[ii-1], N,Ncpu, CosmoParams,  zmin, zmax, area, seed )
        
        else:
           
           job_args=[ (summ_fid, dm, quant_list, epsilon[ii-1], N, Ncpu , CosmoParams,  zmin, zmax, area,  seed) for i in xrange(Ncpu)  ]
        
           results = p.map(worker, job_args)
          # p.close()
         #  p.join()
        
           par_surv=[]
           indx=[]
           par_cov=numpy.array(results[0][-1])
           
           for i in xrange(Ncpu):
              
              
              par_surv.append(results[i][0])
           
              indx.append(results[i][1])
           
           par_surv = numpy.vstack( par_surv )
           indx = numpy.vstack( indx ).flatten()
              
        
        CosmoParams.sdata=par_surv[:]
        
        # store epsilon
        epsilon.append( mquantiles( [ par_surv[ j1 , -1 ] for j1 in xrange( len( par_surv ) ) ], [0.75] )[0]  )
        
        
        #update the weights of the last simulation given the mean and standard deviation from the previous set of simulations
      
        denominator=numpy.array([ sum([CosmoParams.sdata_weights[i]*norm_pdf_multivariate( par_surv[n,:Nparams] , par_surv[ int(indx[ i ]) , :Nparams ],par_cov ) for i in range(N)]) for n in range(N)]	)
      
        numerator = numpy.array([ norm_pdf_multivariate(par_surv[i,:Nparams],CosmoParams.keys_values,CosmoParams.keys_cov ) for i in xrange(N)]) 
    
        weights_new= numerator / denominator
      
        weights=weights_new/weights_new.sum()
        
        CosmoParams.sdata_weights=weights[:]
        
        op1 = open(os.path.join( path1 , 'sequential_ABC_SPT_res_' + str( ii ) + '.dat'), 'w')
        op1.write( '#' ) 
        for par in CosmoParams.keys:
            op1.write('#' + par + '    ' )
        op1.write( 'distance    epsilon\n ' )
        for elem in par_surv:
           for item in elem:
              op1.write( str( item ) + '    ' )
           op1.write( str( epsilon[-2] ) )
           op1.write( '\n' )
        op1.close()

        del op1

      
     

if Ncpu!=1:
  p.join()
  p.close()
else: 
   pass
print "Done!"
print "time=%.4f seconds" % (time.time()-starttime)
