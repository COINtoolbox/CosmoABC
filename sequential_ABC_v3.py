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

from NCountSimul_true_mass_nw import *

import pylab as plt

############################################################
### Fiducial Cosmological parameters: ref. arXiv:1203.5775 (table 5. wCDM CMB+BAO+H0+SNeIa+SPTcl), Tgamma0 and ns are not given.

z_min = 0.3              #minimum redshift
z_max = 1.32             #maximum redshift
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

#time_steps = 10
time_steps = 10 # number of time steps to take

Ninit=50 # initial number of samples 

N=10       # particle sample size after the first iteration

												
epsilon_ini = 1e20		#Starting tolerance

prior_dist = 'normal'

seed = 100  # seed for ncount

##############################################################

path1= '/home/emille/Dropbox/WGC/ABC/github/test_1p_shift/'

#create output directory
if not os.path.exists( path1 ):
    os.makedirs( path1 )

CosmoParams={"H0":H0,"Ob":Omegab,"Om":Omegam,"OL":1-Omegam-Omegab,"Tgamma":Tgamma0,"ns":ns,"sigma8":sigma8,"w":w}

######################################################
#*****************************
#***************************** uncomment this for the 3 parameter test

#Keys=["Om","w", "sigma8"]

#Nparams=len(Keys)

#Keys_values = numpy.array( [  0.25 , -1.01,  0.7 ] )

#Keys_cov0 = numpy.diag( [ 1. , 1., 1. ] )**2.

#Keys_bounds=numpy.array( [ [0.0 , -3.0, 0.3 ] , [1.0-Omegab, 0.0, 1.0] ] )

#########################################################
#*****************************
#***************************** uncomment this for the 2 parameter test

#Keys=["Om","sigma8"]

#Nparams=len(Keys)

#Keys_values = numpy.array( [  0.25 ,  0.7 ] )

#Keys_cov0 = numpy.diag( [ 1. , 1. ] )**2.

#Keys_bounds=numpy.array( [ [0.0 , 0.3 ] , [1.0-Omegab, 1.0] ] )

##################################################################
#*****************************
#***************************** uncomment this for the 1 parameter test

Keys=["Om"]

Nparams=len(Keys)

Keys_values = numpy.array( [  0.5 ] )

Keys_cov0 = numpy.diag( [ 1. ] )**2.

Keys_bounds=numpy.array( [ [0.0  ] , [1.0-Omegab ] ] )

##################################################################

epsilon_list = [ ]


for ii in xrange(time_steps):
    print ii
    if (ii==0):
    
        #new simulation object
        ncount = NCountSimul (z_min, z_max, log ( dm[0] ), log ( 10**16 ), area )

        #Generate fiducial data
        data_fid = numpy.array( ncount.simulation( z_max, seed, CosmoParams)[1] )

        #Calculate summary statistics for fiducial data
        summ_fid = summary_quantile( data_fid, dm, quant_list )
      
        par_surv =  choose_surv_par( summ_fid, dm, quant_list, epsilon_ini, Ninit, CosmoParams, Keys, Keys_cov0,Keys_bounds, prior_dist, z_min, z_max, area, ncount , seed )
      
        par_surv=par_surv[ par_surv[:,-1].argsort()  ]
      
        par_surv=par_surv[xrange(N)]
      
        epsilon_list.append( mquantiles( [ par_surv[ j1 , -1] for j1 in range( len( par_surv ) ) ], [0.75] )[0]  )
        weights_matrix = numpy.array([ [1.0/N if j==i else 0 for j in range(N)  ] for i in range(N)])
      
        weights = [ 1.0/N for j1 in xrange(N) ]

      
        #write results to file
        op0 = open( path1 + '/sequential_ABC_SPT_res_0.dat', 'w')
        op0.write( '#' ) 
        for par in Keys:
            op0.write( par + '    ' )
        op0.write( 'distance    epsilon \n' )
        for elem in par_surv:
           for item in elem:
              op0.write( str( item ) + '    ' )
           op0.write( str( epsilon_ini ) + '\n' )
        op0.close()

      
    else:
      
        #choose one of the instances of simulation set
        indx = weighted_values(range(N), weights,1)

        mean_prev=par_surv[:,:Nparams].mean(axis=0)

      
        par_surv3 = numpy.array([ [par_surv[i][j] - mean_prev[j] for j in range(Nparams) ] for i in range(len (par_surv))]) 

        data_cov = numpy.dot(numpy.transpose(par_surv3), numpy.dot(weights_matrix,par_surv3))

        print 'new_cov = ' + str( data_cov )

        op2 = open( path1 + '/covariance_' + str( ii ) + '.dat','w' )
        for line in data_cov:
            for item in line:
                op2.write( str( item ) + '    ' )
            op2.write( '\n' )
        op2.close()

        del op2
   
        if (Keys_cov0.shape == ( 1, 1 ) ):   
        
            Keys_cov =  2*( data_cov )
      
        else:
           
            Keys_cov=  2*data_cov
      
      
        par_list = []
        indx_list =[]
        for k1 in xrange( N ):

            d1 = 10*epsilon_list[-1] 

            
            while d1 >= epsilon_list[-1]:

                indx = weighted_values(range(N), weights,1)

                 
                #########################################################
                #draw from a multivariate normal distribution

                flag1=numpy.array([ False for i in xrange(Nparams) ])
    
                while (flag1.all() == False ):
                    L = numpy.linalg.cholesky( Keys_cov )
                    norm = numpy.random.normal(size=1*Nparams).reshape(Nparams, 1)
                    new_par = par_surv[ indx ][:-1] + numpy.dot(L, norm).reshape(1,Nparams)
 
                    flag1=((new_par[0]>=Keys_bounds[0]) & (new_par[0]<=Keys_bounds[1])) 
                ###########################################################
                #construct the sampling
 
                for i in xrange( len( Keys ) ):
                    CosmoParams[ Keys[i] ] = new_par[0][ i ]

                if "Om" in Keys:
                    CosmoParams["OL"]=1-CosmoParams["Om"]-CosmoParams["Ob"]
        
                print new_par[0]

                #calculate summary statistics distance to fiducial data  
                d = set_distances( summ_fid, dm, quant_list,  z_min, z_max, area, ncount, seed, CosmoParams  )
                d1 = sum( d )


                print 'tries = ' + str( k1 )

                if d1 <=  epsilon_list[-1]:

                    indx_list.append( indx )
               
                    print '        dist = ' + str( d1 ) + ', epsilon = ' + str( epsilon_list[-1] )

                    result = list( new_par[0] ) + [ d1 ]

                    #add distance to parameter list
                    par_list.append( result )

                ###############################################################
 
        par_surv = numpy.array( par_list )
        epsilon_list.append( mquantiles( [ par_list[ j1 ][ -1 ] for j1 in range( len( par_list ) ) ], [0.75] )[0]  )
      
        #update the weights of the last simulation given the mean and standard deviation from the previous set of simulations
      
        denominator=numpy.array([ sum([weights[i]*norm_pdf_multivariate( par_surv[n,:Nparams] , par_surv[ indx_list[ i ] , :Nparams ],Keys_cov) for i in xrange(N)]) for n in xrange(N)]	)
      
        numerator = numpy.array([ norm_pdf_multivariate(par_surv[i,:Nparams],Keys_values,Keys_cov0 ) for i in xrange(N)]) 
    
        weights_new= numerator / denominator
      
        weights=weights_new/weights_new.sum()

        weights_matrix = numpy.array([ [ weights[ j ] if j==i else 0 for j in range(N)  ] for i in range(N)])
   

        #write results to file
        op1 = open(path1 + '/sequential_ABC_SPT_res_' + str( ii ) + '.dat', 'w')
        op1.write( '#' ) 
        for par in Keys:
            op1.write( par + '    ' )
        op1.write( 'distance    epsilon\n ' )
        for elem in par_surv:
           for item in elem:
              op1.write( str( item ) + '    ' )
           op1.write( str( epsilon_list[-2] ) )
           op1.write( '\n' )
        op1.close()

        del op1


