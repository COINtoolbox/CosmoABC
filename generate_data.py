from math import *
from gi.repository import GObject
import matplotlib.pyplot as plt
from gi.repository import NumCosmo as Nc
from gi.repository import NumCosmoMath as Ncm

from NCountSimul_true_mass_nw_v4 import *


def generate_data(name):
   
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
   
   #sky area
   area = 2500
   
   #seed
   
   seed = 100 
   
   CosmoParams=ChooseParamsInput()

   CosmoParams.params={"H0":H0,"Ob":Omegab,"Om":Omegam,"OL":1.-Omegam-Omegab,"Tgamma":Tgamma0,"ns":ns,"sigma8":0.8,"w":-1.0}

   ncount = NCountSimul (zmin, zmax, log ( dm[0] ), log ( 10**16 ), area )
   
   data_fid = numpy.array( ncount.simulation( zmax, seed, CosmoParams.params)[1] )
   
   f=open(name+".params","w")
   
   f.write("Parameters used\n")
   f.write("zmin=%.3f\n" % zmin)
   f.write("zmax=%.3f\n" % zmax)
   f.write("H0=%.2f\n" % H0)
   f.write("Omega_baryons=%.3f\n" % Omegab)
   f.write("Omega_matter=%.3f\n" % Omegam)
   f.write("Tgamma0=%.3f\n" % Tgamma0)
   f.write("zmin=%.3f\n" % ns)
   f.write("sigma8=%.3f\n" % sigma8)
   f.write("w=%.3f\n" % w)
   f.write("\n")
   f.write("Mass bins\n")
   for i in dm:
     
      f.write("%.4e\t" % i)

   f.close()
   
   
   numpy.savetxt(name+".dat",data_fid,fmt="%.4f\t %.4f")
   
   return data_fid
   

