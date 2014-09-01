import numpy as np

#from NCount_Simul import *
from NCount_Simul_SPT import *



def generate_data(name):
   
   
   
   ############################################################
   ### Fiducial Cosmological parameters: ref. arXiv:1203.5775 (table 5. wCDM CMB+BAO+H0+SNeIa+SPTcl), Tgamma0 and ns are not given.

   zmin = 0.3              #minimum redshift
   zmax = 1.32             #maximum redshift
   H0 = 70.0             #Hubble parameter
   Omegam = 0.262           #Dark matter density
   Omegab = 0.0439          #Baryon density
   Tgamma0 = 2.725          #Radiation temperature today
   ns = 0.97                #spectral index 
   sigma8 = 0.807           #sigma8
   w = -1.0                #Dark energy equation of state     

#mass bin
#dm = [5*10**13, 10**14, 10**14.25, 10**14.5, 10**14.75,  10**15, 10**15.25,  10**15.5, 10**15.75 ]
   dm = [10**14.3, 10**14.5, 10**14.7,  10**14.9, 10**15.1, 10**15.3,  10**15.5, 10**15.7 ]
   
   #sky area
   area = 2500
   
   #seed
   
   seed = 100 
   

   CP={"H0":H0,"Ob":Omegab,"Om":Omegam,"OL":1.-Omegam-Omegab,"Tgamma":Tgamma0,"ns":ns,"sigma8":sigma8,"w":w}

   ncount = NCountSimul (zmin, zmax, np.log ( dm[0] ), np.log ( 10**16 ), area, seed )
   
   data_fid = np.array( ncount.simulation( zmax,  CP )[1] )
   
   f=open(name+".params","w")
   
   f.write("Parameters used\n")
   f.write("zmin=%.3f\n" % zmin)
   f.write("zmax=%.3f" % zmax)
   f.write("H0=%.2f" % H0)
   f.write("Omega_baryons=%.3f" % Omegab)
   f.write("Omega_matter=%.3f" % Omegam)
   f.write("Tgamma0=%.3f" % Tgamma0)
   f.write("zmin=%.3f" % ns)
   f.write("sigma8=%.3f" % sigma8)
   f.write("w=%.3f" % w)
   
   
   np.savetxt(name+".dat",data_fid,fmt="%.4f\t %.4f")
   
   return data_fid
   

