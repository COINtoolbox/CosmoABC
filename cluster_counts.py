#!/usr/bin/python

import numpy as np
import random
import math

from scipy.stats import anderson_ksamp
from scipy.stats.mstats import mquantiles

from ABC import ABC, Distributions, StoreInfo, RHO
from NCount_Simul import NCountSimul

## Fiducial Cosmological parameters: ref. arXiv:1203.5775 (table 5. wCDM CMB+BAO+H0+SNeIa+SPTcl), Tgamma0 and ns are not given.

zmin = 0.3              #minimum redshift
zmax = 1.32             #maximum redshift
H0 = 70.0             #Hubble parameter
Omegam = 0.262  +0.05         #Dark matter density
Omegab = 0.0439          #Baryon density
Tgamma0 = 2.725          #Radiation temperature today
ns = 0.97                #spectral index 
sigma8 = 0.807  +0.05         #sigma8
w = -1.0                #Dark energy equation of state     


CosmoParams={"H0":H0,"Ob":Omegab,"Om":Omegam,"OL":1.-Omegam-Omegab,"Tgamma":Tgamma0,"ns":ns,"sigma8":sigma8,"w":w}


 
############ Load observational data ##########################

data_obs = np.loadtxt ( "Mock_Data_True_Mass.dat" )

##############################################################




#############################################################################################################################

# Cluster mass stuff!
dm_choose = [10**14.3, 10**14.5, 10**14.7,  10**14.9, 10**15.1, 10**15.3,  10**15.5, 10**15.7 ]

mass_min = 10**14.3
mass_max = 10**16
###############################################################################################################################

#quantile list
quant_list = [ 0.02, 0.09, 0.25, 0.5, 0.75, 0.91, 0.98]

###########################################################################################################################

#sky area
area = 2500

observable = 'true_mass'

if observable == 'true_mass':
    from NCount_Simul import NCountSimul
    
else:
    from NCount_Simul_SPT import NCountSimul


#choose observable
#observable = 'true_mass'

#if observable == 'true_mass':
#    from NCountSimul_true_mass_nw import *
    
#else:
#    from NCountSimul_SPT import *
     

nclusters = NCountSimul( zmin, zmax, np.log( mass_min ), np.log( mass_max ), area,False)



def model (p,*args):
	
	CosmoParams[ "Om" ] = p[ 0 ] 
	CosmoParams[ "Ode "] = 1.-CosmoParams[ "Om" ]- CosmoParams[ "Ob" ]
	CosmoParams[ "sigma8" ] = p[1]
	
	
	
	return np.array( nclusters.simulation( zmax,  CosmoParams )[1] )
	


def SumStats ( data ):
	
	mass = np.array( [ np.log( item ) for item in dm_choose ] )
	
	
	data_bin = np.array([ [ item[0] for item in data  if item[1] >= mass[ i ] ] for i in range (len(mass) -1) ])
	
	res = np.array( [ mquantiles( elem, prob=quant_list ) if len( elem ) > 0  else [ 0 for jj in quant_list]  for elem in data_bin ] )
	
	if sum( [ len( data_bin[ k ] ) for k in range( len( data_bin ) ) ] ) > 0:
		
		pop = np.array ( [ float( len( data_bin[ i ] ) )/sum( [ len( data_bin[ k ] ) \
		for k in range( len( data_bin ) ) ] )  for i in range( len( data_bin ) ) ] )
		
	else:
		
		pop = np. array( [ 0 for i in range( len( data_bin ) ) ] )
	
	return  res ,  pop 
	
	
	
def dist( summary_fid , summary_sim  ):
	
	#distance = np.sum( np.sqrt( (summary_fid[0] - summary_sim[0])**2. )   )
	distance = np.sum([np.sqrt(np.sum((summary_fid[0][j]-summary_sim[0][j])**2)) for j in xrange(len(summary_fid[0]))])
	
	#distance =  np.sum( [ np.sqrt( sum( [ ( summary_fid[0][ i ][ j ] -\
	 # summary_sim[0][ i ][ j ] ) ** 2  for j in range( len( summary_fid[0][ i ] ) ) ] ) ) for i in range( len( summary_fid[0] ) ) ] )
	  
	return distance, summary_fid[1], summary_sim[1]
	  



	#d = dist=RHO ( data , model , SumStats , dist , (False,) )

d = RHO ( data_obs , model , SumStats ,  dist )

Prior=Distributions( "Normal" , np.array( [ 0.3 , 0.7 ] ), np.diag([1,1]), np.array( [[0.0,0.4],[1-Omegab,1.0]]) )
	
#P=StoreInfo(10,2,100,Prior, d)

N_iter= 20

N=250

Ncpu = 2

P=StoreInfo(N_iter,Ncpu,N)
sc = ABC(P,Prior,d)
	
out=sc.sampler("SimpleTest")
