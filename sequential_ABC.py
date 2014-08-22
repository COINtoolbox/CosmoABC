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
import shutil


import pylab as plt


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
observable = 'true_mass'


#mass bin
#dm = [5*10**13, 10**14, 10**14.25, 10**14.5, 10**14.75,  10**15, 10**15.25,  10**15.5, 10**15.75 ]
dm_choose = [10**14.3, 10**14.5, 10**14.7,  10**14.9, 10**15.1, 10**15.3,  10**15.5, 10**15.7 ]

mass_min = 10**14.3
mass_max = 10**16

#quantile list
quant_list = [ 0.02, 0.09, 0.25, 0.5, 0.75, 0.91, 0.98]

#sky area
area = 2500

#######
# path and name to mock data file
mock_data = "/home/emille/Dropbox/WGC/ABC/data/Mock_Data.dat"  


Ninit=2000 # initial number of samples 

N=1000      # particle sample size after the first iteration

												
epsilon_ini = [1e20, 1e20]		#Starting tolerance

#if seed == False, use random 
seed = False  # seed for ncount

Nparams=2

##############################################################


# Load mock data
data_fid = numpy.loadtxt( mock_data )
 
#path to results directory
path1="/home/emille/Dropbox/WGC/ABC/RESULTS/" + observable + '/' + str( Nparams ) + 'p/om_w/'
output_param_file_root = 'SMC_ABC_' + observable + '_'  
output_covariance_file = 'covariance_evol_' + observable + '.dat'

#create output directory
if not os.path.exists( path1 ):
    os.makedirs( path1 )

#copy input file to directory for this run
shutil.copy2( 'sequential_ABC.py', path1 + 'input.py' )   

if observable == 'true_mass':
    from NCountSimul_true_mass_nw import *
    shutil.copy2( 'NCountSimul_true_mass_nw.py', path1 + 'functions.py' )   
else:
    from NCountSimul_SPT import *
    shutil.copy2( 'NCountSimul_SPT.py', path1 + 'functions.py' )  


   
CosmoParams=ChooseParamsInput()

CosmoParams.params={"H0":H0,"Ob":Omegab,"Om":Omegam,"OL":1.-Omegam-Omegab,"Tgamma":Tgamma0,"ns":ns,"sigma8":0.8,"w":-1.0}

#CosmoParams.keys=["Om","w","sigma8"]
CosmoParams.keys=["Om","sigma8"]

#CosmoParams.keys_values=numpy.array( [  0.25 , -1.01,  0.7 ] )
CosmoParams.keys_values=numpy.array( [  0.25 ,  0.7 ] )

#CosmoParams.keys_cov=numpy.diag( [1.0 , 1.0, 1.0 ] )**2.
CosmoParams.keys_cov=numpy.diag( [1.0 , 1.0 ] )**2.


#CosmoParams.keys_bounds=numpy.array( [ [0.0 , -10.0, 0.3 ] , [1.0-Omegab, 1.0, 1.0] ] )
CosmoParams.keys_bounds=numpy.array( [ [0.0 , 0.3 ] , [1.0-Omegab,  1.0] ] )


#desired final variance (percentage in relation to the initial variance)
#CosmoParams.desired_variance = [ 0.1, 0.1, 0.1]
#CosmoParams.desired_median = [0.1,0.1,0.1]

CosmoParams.desired_variance = [ 0.1,  0.1]
CosmoParams.desired_median = [0.1,0.1]


#CosmoParams.keys=["Om"]

#CosmoParams.keys_values=numpy.array( [  0.25  ] )
#CosmoParams.keys_cov=numpy.diag( [ 0.5  ] )**2.

#CosmoParams.keys_bounds=numpy.array( [ [0.0  ] , [1.0-Omegab] ] )

CosmoParams.prior_dist="normal"

Nparams=len( CosmoParams.keys )


   
op3 = open( path1 + output_covariance_file, 'w')
for item in CosmoParams.keys:
    op3.write( item + '    ' )

op3.write( '\n' )   

epsilon = [ ]

var_flag = 0
#for ii in xrange(time_steps):
cont = 0
while var_flag < Nparams:

    print 'cont  = ' + str( cont )

    if ( cont == 0 ):
    
        #new simulation object
        ncount = NCountSimul (zmin, zmax, log ( mass_min ), log ( mass_max ), area )
  

        #Generate fiducial data
        data_fid = numpy.array( ncount.simulation( zmax, seed, CosmoParams.params)[1] )


        if observable == 'true_mass':
            dm = dm_choose
        else:
            dm = mquantiles( data_fid[:,1], prob=quant_list[1:-1] )


        #keep number of fiducial data set
        nobjs_fid = len( data_fid )

        #Calculate summary statistics for fiducial data
        summ_fid = summary_quantile( data_fid, dm, quant_list )
        
        par_surv, indx,par_cov = choose_surv_par( summ_fid, dm, quant_list, epsilon_ini, Ninit, CosmoParams, zmin, zmax, area, ncount, seed, nobjs_fid )
        
        cont = cont + 1 
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
        
        CosmoParams.variance = [ numpy.std( CosmoParams.sdata[:, i] ) for i in range( Nparams )]
        CosmoParams.median =  [ numpy.median( CosmoParams.sdata[:, i] ) for i in range( Nparams )]
        print 'covariances = ' + str( CosmoParams.variance )

        for elem in CosmoParams.variance:
            op3.write( str( elem ) + '    ' )
        op3.write( '\n' )     
   

        #write results to file
        op0 = open(os.path.join( path1 , output_param_file_root + '0.dat'), 'w')
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
    
        print 'cont = ' + str( cont )

        par_surv, indx, par_cov = choose_surv_par( summ_fid, dm, quant_list, epsilon[-1], N, CosmoParams,  zmin, zmax, area, ncount, seed, nobjs_fid )
        CosmoParams.sdata=par_surv[:]
        
        #check if desired variance was achieved
        new_cov = [ numpy.std( CosmoParams.sdata[:, i2] ) for i2 in range( Nparams ) ]
        new_median = [ numpy.median( CosmoParams.sdata[:, i2] ) for i2 in range( Nparams ) ]
        CosmoParams.variance_diff = [ abs( new_cov[ i2 ] - CosmoParams.variance[ i2 ] )/CosmoParams.variance[ i2 ] for i2 in range( Nparams )]
        CosmoParams.median_diff = [ abs( new_median[ i2 ] - CosmoParams.median[ i2 ] )/CosmoParams.median[ i2 ] for i2 in range( Nparams )]

        var_flag = 0
        for j4 in range( Nparams ):
            if CosmoParams.variance_diff[ j4 ] < CosmoParams.desired_variance[ j4 ] and CosmoParams.median_diff < CosmoParams.desired_median[ j4 ]:
                var_flag = var_flag + 1
        
        print 'cov_diff = ' + str( CosmoParams.variance_diff )
        print 'median_diff  = ' + str( CosmoParams.median_diff )
 

        for elem in CosmoParams.variance:
            op3.write( str( elem ) + '    ' )
        op3.write( '\n' )
        
        
        if var_flag < Nparams:

            CosmoParams.variance = new_cov               
            CosmoPArams.median = new_median

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

        
            op1 = open(os.path.join( path1 , output_param_file_root + str( cont ) + '.dat'), 'w')
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
            cont = cont + 1 
            
      
op3.close()      

