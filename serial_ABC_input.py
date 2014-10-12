from math import *
from gi.repository import GObject
import matplotlib.pyplot as plt
from gi.repository import NumCosmo as Nc
from gi.repository import NumCosmoMath as Ncm

from scipy.stats.mstats import mquantiles
from scipy.stats import norm
from serial_ABC_functions import *

import numpy
import random 
import os 
import shutil

from multiprocessing import Pool
import time

starttime=time.time()
############################################################
### Fiducial Cosmological parameters: ref. arXiv:1203.5775 (table 5. wCDM CMB+BAO+H0+SNeIa+SPTcl), Tgamma0 and ns are not given.
zmin = 0.3              #minimum redshift
zmax = 1.32             #maximum redshift
H0 = 71.15               #Hubble parameter
Omegam = 0.262           #Dark matter density
Omegab = 0.044          #Baryon density
Tgamma0 = 2.725          #Radiation temperature today
ns = 0.97                #spectral index 
sigma8 = 0.807           #sigma8
w = -1.01                #Dark energy equation of state     

#choose observable
observable = 'SZ'

mass_min = 10**14.3
mass_max = 10**16

#quantile list
quant_list = [ 0.02, 0.09, 0.25, 0.5, 0.75, 0.91, 0.98]


#quantile threshold
choose_quantile = 0.75

#sky area
area = 2500

#######
# path and name to data file
mock_data = "dataset_ABC_seed123.dat"

N=200      # particle sample size after the first iteration

												
epsilon_ini = [1e20, 1e20]		#Starting tolerance

#if seed == False, use random 
seed = False  # seed for ncount

var_tol = 0.075

#options are 'Om', 'w' and 'sigma8'
parm_to_fit=['Om']

##############################################################
output_param_file_root = 'SMC_ABC_' + observable + '_'  
output_covariance_file = 'covariance_evol_' + observable + '.dat'

CosmoParams=ChooseParamsInput()

CosmoParams.params={"H0":H0,"Ob":Omegab,"Om":Omegam,"OL":1.-Omegam-Omegab,"Tgamma":Tgamma0,"ns":ns,"sigma8":sigma8,"w":w}

if len( parm_to_fit ) == 3:
    CosmoParams.keys=["Om","w","sigma8"]
    CosmoParams.keys_values=numpy.array( [  0.25 , -0.8,  0.7 ] )
    CosmoParams.keys_cov=numpy.diag( [1.0 , 1.0, 1.0 ] )**2.
    CosmoParams.keys_bounds=numpy.array( [ [0.0 , -10.0, 0.6 ] , [1.0-Omegab, 1.0,  1.0] ] )
    CosmoParams.desired_variance = [ var_tol,  var_tol, var_tol]
    CosmoParams.desired_median = [var_tol, var_tol, var_tol]
    CosmoParams.desired_mean = [var_tol, var_tol, var_tol]

elif ( 'Om' in parm_to_fit ) and ( 'w' in parm_to_fit ):
    CosmoParams.keys=["Om","w"]
    CosmoParams.keys_values=numpy.array( [  0.25 , -0.8 ] )
    CosmoParams.keys_cov=numpy.diag( [1.0 , 1.0] )**2.
    CosmoParams.keys_bounds=numpy.array( [ [0.0 , -10.0] , [1.0-Omegab, 1.0 ] ] )
    CosmoParams.desired_variance = [ var_tol,  var_tol ]
    CosmoParams.desired_median = [var_tol, var_tol ]
    CosmoParams.desired_mean = [var_tol, var_tol ]

elif ( 'Om' in parm_to_fit ) and ( 'sigma8' in parm_to_fit ):
    CosmoParams.keys=["Om","sigma8"]
    CosmoParams.keys_values=numpy.array( [  0.25 , 0.7 ] )
    CosmoParams.keys_cov=numpy.diag( [1.0 , 1.0] )**2.
    CosmoParams.keys_bounds=numpy.array( [ [0.0 , 0.6] , [1.0-Omegab, 1.0 ] ] )
    CosmoParams.desired_variance = [ var_tol,  var_tol ]
    CosmoParams.desired_median = [var_tol, var_tol ]
    CosmoParams.desired_mean = [var_tol, var_tol ]

elif ( 'w' in parm_to_fit ) and ( 'sigma8' in parm_to_fit ):
    CosmoParams.keys=["w","sigma8"]
    CosmoParams.keys_values=numpy.array( [  -0.8, 0.7 ] )
    CosmoParams.keys_cov=numpy.diag( [1.0 , 1.0] )**2.
    CosmoParams.keys_bounds=numpy.array( [ [-10.0, 0.6] , [1.0, 1.0 ] ] )
    CosmoParams.desired_variance = [ var_tol,  var_tol ]
    CosmoParams.desired_median = [var_tol, var_tol ]
    CosmoParams.desired_mean = [var_tol, var_tol ]     

elif len( parm_to_fit ) == 1 and 'Om' in parm_to_fit:
    CosmoParams.keys=["Om"]
    CosmoParams.keys_values=numpy.array( [ 0.25 ] )
    CosmoParams.keys_cov=numpy.diag( [1.0] )**2.
    CosmoParams.keys_bounds=numpy.array( [ [0.0] , [1.0-Omegab] ] )
    CosmoParams.desired_variance = [ var_tol ]
    CosmoParams.desired_median = [var_tol ]
    CosmoParams.desired_mean = [var_tol ]     

elif len( parm_to_fit ) == 1 and 'w' in parm_to_fit:
    CosmoParams.keys=["w"]
    CosmoParams.keys_values=numpy.array( [ -0.8 ] )
    CosmoParams.keys_cov=numpy.diag( [1.0] )**2.
    CosmoParams.keys_bounds=numpy.array( [ [-10.0] , [1.0] ] )
    CosmoParams.desired_variance = [ var_tol ]
    CosmoParams.desired_median = [var_tol ]
    CosmoParams.desired_mean = [var_tol ]     

elif len( parm_to_fit ) == 1 and 'sigma8' in parm_to_fit:
    CosmoParams.keys=["sigma8"]
    CosmoParams.keys_values=numpy.array( [ 0.7 ] )
    CosmoParams.keys_cov=numpy.diag( [1.0] )**2.
    CosmoParams.keys_bounds=numpy.array( [ [0.6] , [1.0] ] )
    CosmoParams.desired_variance = [ var_tol ]
    CosmoParams.desired_median = [var_tol ]
    CosmoParams.desired_mean = [var_tol ]     


CosmoParams.prior_dist="normal"

Nparams=len( CosmoParams.keys )

#read "observed data"
data_fid = numpy.loadtxt( mock_data )
   
op3 = open( output_covariance_file, 'w')
for item in CosmoParams.keys:
    op3.write( item + '    ' )

op3.write( '\n' )   

epsilon = [ ]

var_flag = 0
cont = 0




while var_flag < Nparams:

    print 'cont  = ' + str( cont )

    if ( cont == 0 ):

        if observable == 'true_mass':
            dm_choose = [ numpy.log( mass_min ) ]
        else:
            dm_choose = [ min( data_fid[:,1] ) ]

        mq = mquantiles( data_fid[:,1], prob=quant_list )
        for it in mq:
            dm_choose.append( it )

        if observable == 'true_mass':
            dm_choose.append( numpy.log( mass_max ) )
        else:
            dm_choose.append( max( data_fid[:,1] ) ) 



        #if observable == 'true_mass':
        dm = dm_choose

        #keep number of fiducial data set
        nobjs_fid = len( data_fid )

        #Calculate summary statistics for fiducial data
        summ_fid = summary_quantile( data_fid, dm, quant_list )
                
        par_surv, indx,par_cov = choose_surv_par (summ_fid, dm, quant_list, epsilon_ini, 10*N, CosmoParams,  zmin, zmax, area,  seed, nobjs_fid, [numpy.log(mass_min), numpy.log(mass_max)], observable)
        
   
        cont = cont + 1 
      
        #sort distances according to the sum of the summary statistics distance and population tests
        final_dist = numpy.array([ sum( [ par_surv[ i ][ l1 ]/numpy.mean(par_surv[:, l1])  for l1 in xrange( -2, 0) ]) for i in range( len( par_surv ) ) ])
      
        #par_surv=par_surv[ par_surv[:,-2].argsort()  ]
      
        param_sorted = final_dist.argsort() 
      
        #par_surv=par_surv[xrange(N)]

        CosmoParams.sdata_weights = numpy.ones( N )/N
        par_surv = numpy.array([ par_surv[ i3 ] for i3 in param_sorted[:N] ])
       

        epsilon.append( [ mquantiles( [ par_surv[ j1 , j2] for j1 in xrange( len( par_surv ) ) ], [ choose_quantile ] )[0] for j2 in range(-2,0) ]  )
        
        CosmoParams.sdata=par_surv[:]

        op4 = open('weights_' + str( cont ) + '.dat', 'w' )
        op4.write( str( cont ) + '    ' )
        for iii in CosmoParams.sdata_weights:
            op4.write( str( iii ) + '    ' )
        op4.close()
        
        CosmoParams.variance = [ numpy.std( CosmoParams.sdata[:, i] ) for i in range( Nparams )]
        CosmoParams.median =  [ numpy.median( CosmoParams.sdata[:, i] ) for i in range( Nparams )]
        CosmoParams.mean =  [ numpy.mean( CosmoParams.sdata[:, i] ) for i in range( Nparams )]
        print 'covariances = ' + str( CosmoParams.variance )
        print 'mean = ' + str( CosmoParams.mean )
        print 'median = ' + str( CosmoParams.median )

        for elem in CosmoParams.variance:
            op3.write( str( elem ) + '    ' )
        op3.write( '\n' )     
   

        #write results to file
        op0 = open(os.path.join( output_param_file_root + '0.dat'), 'w')
        for par in CosmoParams.keys:
            op0.write( par + '    ' )
        op0.write( 'distancia1    distancia2 \n' )
        for elem in par_surv:
           for item in elem:
              op0.write( str( item ) + '    ' )
           op0.write('\n' )
        op0.close()

      
    else:
    
        print 'cont = ' + str( cont )

        par_surv_old = par_surv
        
           
        par_surv, indx,par_cov = choose_surv_par (summ_fid, dm, quant_list, epsilon[ -1 ], N, CosmoParams,  zmin, zmax, area,  seed, nobjs_fid, [numpy.log(mass_min), numpy.log(mass_max)], observable)
        
        CosmoParams.sdata=par_surv[:]
        
        #check if desired variance was achieved
        new_cov = [ numpy.std( CosmoParams.sdata[:, i2] ) for i2 in range( Nparams ) ]
        new_median = [ numpy.median( CosmoParams.sdata[:, i2] ) for i2 in range( Nparams ) ]
        new_mean = [ numpy.mean( CosmoParams.sdata[:, i2] ) for i2 in range( Nparams ) ]
        CosmoParams.variance_diff = [ abs( new_cov[ i2 ] - CosmoParams.variance[ i2 ] )/CosmoParams.variance[ i2 ] for i2 in range( Nparams )]
        CosmoParams.median_diff = [ abs( new_median[ i2 ] - CosmoParams.median[ i2 ] )/CosmoParams.median[ i2 ] for i2 in range( Nparams )]
        CosmoParams.mean_diff = [ abs( new_mean[ i2 ] - CosmoParams.mean[ i2 ] )/CosmoParams.mean[ i2 ] for i2 in range( Nparams )]

        var_flag = 0
        
        
        
        for j4 in range( Nparams ):
            if CosmoParams.variance_diff[ j4 ] < CosmoParams.desired_variance[ j4 ] and CosmoParams.mean_diff[ j4 ] < CosmoParams.desired_median[ j4 ]:
                var_flag = var_flag + 1
        
        print 'cov_diff = ' + str( CosmoParams.variance_diff )
        print 'median_diff  = ' + str( CosmoParams.median_diff )
        print 'mean_diff  = ' + str( CosmoParams.mean_diff )
 

        for elem in CosmoParams.variance:
            op3.write( str( elem ) + '    ' )
        op3.write( '\n' )
        
        
        if var_flag < Nparams:

            CosmoParams.variance = new_cov               

            CosmoParams.median = new_median
            # store epsilon
            epsilon.append( [ mquantiles( [ par_surv[ j1 , j2 ] for j1 in xrange( len( par_surv ) ) ], [ choose_quantile ] )[0] for j2 in range(-2,0) ] )
        
        
            #update the weights of the last simulation given the mean and standard deviation from the previous set of simulations  
            denominator=numpy.array([ sum([CosmoParams.sdata_weights[i]*norm_pdf_multivariate( par_surv[n,:Nparams] , par_surv_old[ i  , :Nparams ], par_cov ) for i in xrange(N)]) for n in xrange(N)]	)
      
            numerator = numpy.array([ norm_pdf_multivariate(par_surv[i,:Nparams],CosmoParams.keys_values,CosmoParams.keys_cov ) for i in xrange(N)]) 
    
            weights_new= numerator / denominator 
      
            weights=weights_new/weights_new.sum()
        
            CosmoParams.sdata_weights=weights[:]
            
        op4 = open('weights_' + str( cont ) + '.dat', 'w' )    
        op4.write( str( cont ) + '    ' )
        for iii in CosmoParams.sdata_weights:
            op4.write( str( iii ) + '    ' )
        op4.close() 
        
        op1 = open(os.path.join( output_param_file_root + str( cont ) + '.dat'), 'w')
        for par in CosmoParams.keys:
            op1.write( par + '    ' )
        op1.write( 'distancia1    distancia2    epsilon1    epsilon2\n ' )
        for elem in par_surv:
           for item in elem:
              op1.write( str( item ) + '    ' )
           op1.write( str( epsilon[-2][0] ) + '    ' + str( epsilon[-2][1] ) + '\n' )
        op1.close()

        del op1
        
        cont = cont + 1 
            
      
op3.close()      


print "Done!"
print "time=%.4f seconds" % (time.time()-starttime)
