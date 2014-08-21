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

from NCountSimul_true_mass_nw  import *

import pylab as plt

################################################################




#make random changes

############################################################
### Fiducial Cosmological parameters: ref. arXiv:1203.5775 (table 5. wCDM CMB+BAO+H0+SNeIa+SPTcl), Tgamma0 and ns are not given.

zmin = 0.3              #minimum redshift
zmax = 1.32             #maximum redshift
H0 = 71.15 #Hubble parameter
Omegam = 0.262 #Dark matter density
Omegab = 0.0439 #Baryon density
Tgamma0 = 2.725 #Radiation temperature today
ns = 0.97 #spectral index
sigma8 = 0.807 #sigma8
w = -1.01 #Dark energy equation of state 

#sky area
area = 2500

time_steps =10 # number of time steps to take

Ninit=1000 # initial number of samples 

N=250      # particle sample size after the first iteration

seed=100

path1="teste_3p_NO_weight_dist"



#quantile list
quant_list = [ 0.02, 0.09, 0.25, 0.5, 0.75, 0.91, 0.98]

epsilon_ini = [1e20, 1e20]		#Starting tolerance

mock_data = "Mock_Data.dat"  # path and name to mock data file

# Load mock data
data_fid = numpy.loadtxt( mock_data )

#keep number of fiducial data set
nobjs_fid = len( data_fid )

#mass bin
#dm = [5*10**13, 10**14, 10**14.25, 10**14.5, 10**14.75,  10**15, 10**15.25,  10**15.5, 10**15.75 ]
dm = [10**14.3, 10**14.5, 10**14.7,  10**14.9, 10**15.1, 10**15.3,  10**15.5, 10**15.7 ]

#Calculate summary statistics for fiducial data
summ_fid = summary_quantile( data_fid, dm, quant_list )
##############################################################
#create output directory
if not os.path.exists( path1 ):
    os.makedirs( path1 )
    
    
    
CosmoParams=ChooseParamsInput()

CosmoParams.params={"H0":H0,"Ob":Omegab,"Om":Omegam,"OL":1.-Omegam-Omegab,"Tgamma":Tgamma0,"ns":ns,"sigma8":0.8,"w":-1.0}

CosmoParams.keys=["Om","w","sigma8"]

CosmoParams.keys_values=numpy.array( [  0.25 , -1.01,  0.7 ] )

CosmoParams.keys_cov=numpy.diag( [1.0 , 1.0, 1.0 ] )**2.

CosmoParams.keys_bounds=numpy.array( [ [0.0 , -3.0, 0.3 ] , [1.0-Omegab, 0.0, 1.0] ] )




CosmoParams.prior_dist="normal"

Nparams=len( CosmoParams.keys )
    
#new simulation object
ncount = NCountSimul (zmin, zmax, log ( dm[0] ), log ( 10**16 ), area )

epsilon = [ ]


for ii in xrange(time_steps):
    print ii
    if (ii==0):
    
        
        
        par_surv, indx,par_cov = choose_surv_par( summ_fid, dm, quant_list, epsilon_ini, Ninit, CosmoParams, zmin, zmax, area, ncount, seed, nobjs_fid )
        
        # print "Data generated"
        # print par_surv
        # print "indexes"
        # print indx
        # print "covariance"
        # print par_cov

        #sort distances according to the sum of the summary statistics distance and population tests
        final_dist = numpy.array([ sum( [ par_surv[ i ][ l1 ]/numpy.mean(par_surv[:, l1])  for l1 in xrange( -2, 0) ]) for i in range( len( par_surv ) ) ])
      
        #par_surv=par_surv[ par_surv[:,-2].argsort()  ]
      
        param_sorted = final_dist.argsort() 
      
        #par_surv=par_surv[xrange(N)]
        par_surv = numpy.array([ par_surv[ i3 ] for i3 in param_sorted[:N] ])
                 

        epsilon.append( [ mquantiles( [ par_surv[ j1 , j2] for j1 in xrange( len( par_surv ) ) ], [0.75] )[0] for j2 in range(-2,0) ]  )
        
        CosmoParams.sdata=par_surv[:]
        
        CosmoParams.sdata_weights = numpy.ones( N )/N
        
        
        #write results to file
        op0 = open(os.path.join( path1 , 'sequential_ABC_SPT_res_0.dat'), 'w')
        op0.write( '#' ) 
        for par in CosmoParams.keys:
            op0.write('#' + par + '    ' )
        op0.write( 'distancia1    distancia2 \n' )
        for elem in par_surv:
           for item in elem:
              op0.write( str( item ) + '    ' )
           op0.write('\n' )
        op0.close()

      
    else:
    
        
        par_surv, indx, par_cov = choose_surv_par( summ_fid, dm, quant_list, epsilon[ii-1], N, CosmoParams,  zmin, zmax, area, ncount, seed, nobjs_fid )
        CosmoParams.sdata=par_surv[:]
        
        # store epsilon
        epsilon.append( [ mquantiles( [ par_surv[ j1 , j2 ] for j1 in xrange( len( par_surv ) ) ], [0.75] )[0] for j2 in range(-2,0) ] )
        
        
        #update the weights of the last simulation given the mean and standard deviation from the previous set of simulations
        ###########
           ####updated number of epsilon parameters     
        denominator=numpy.array([ sum([CosmoParams.sdata_weights[i]*norm_pdf_multivariate( par_surv[n,:Nparams] , par_surv[ indx[ i ] , :Nparams ],par_cov ) for i in xrange(N)]) for n in xrange(N)]	)
      
        numerator = numpy.array([ norm_pdf_multivariate(par_surv[i,:Nparams],CosmoParams.keys_values,CosmoParams.keys_cov ) for i in xrange(N)]) 
    
        weights_new= numerator / denominator
      
        weights=weights_new/weights_new.sum()
        
        CosmoParams.sdata_weights=weights[:]
        
        op1 = open(os.path.join( path1 , 'sequential_ABC_SPT_res_' + str( ii ) + '.dat'), 'w')
        op1.write( '#' ) 
        for par in CosmoParams.keys:
            op1.write('#' + par + '    ' )
        op1.write( 'distancia1    distancia2\n ' )
        for elem in par_surv:
           for item in elem:
              op1.write( str( item ) + '    ' )
           op1.write( '\n' )
        op1.close()

        del op1

      
      

