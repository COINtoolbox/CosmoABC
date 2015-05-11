"""
Functions for catalog simulations using NumCosmo. 
"""

import numpy as np
import random
import math
from numpy import pi

from gi.repository import NumCosmo as Nc
from gi.repository import NumCosmoMath as Ncm

############################################################################
###	From NumCosmo: for simulation	####################################
############################################################################


class NCountSimul:

  def __init__(self, z_min, z_max, z_sigma, lnM_min, 
                lnM_max, area, observable):
    Ncm.cfg_init()
    self.cosmo = Nc.HICosmo.new_from_name(Nc.HICosmo, "NcHICosmoDEXcdm")
    dist = Nc.Distance.new(z_max*1.5)
    wp =  Nc.Window.new_from_name("NcWindowTophat")
    tf = Nc.TransferFunc.new_from_name("NcTransferFuncEH")
    vp = Nc.MatterVar.new(Nc.MatterVarStrategy.FFT, wp, tf)
    gf = Nc.GrowthFunc.new()
    mulf = Nc.MultiplicityFunc.new_from_name("NcMultiplicityFuncTinkerCrit{'Delta':<500.0>}")
    mf = Nc.MassFunction.new(dist, vp, gf, mulf)
    
    if observable == 'SZ':
        cluster_m = Nc.ClusterMass.new_from_name("NcClusterMassBenson{'M0':<3e14>, 'z0':<0.6>, 'signif-obs-min':<5.0>,'Asz':<6.24>, 'Bsz':<1.33>, 'Csz':<0.83>, 'Dsz':<0.24>}")
    elif observable == 'true_mass':
        cluster_m = Nc.ClusterMass.new_from_name("NcClusterMassNodist{'lnM-min':<% 20.15g>, 'lnM-max':<% 20.15g>}" % (lnM_min, lnM_max))
    else:
        raise NameError('Invalid observable choice. Should be ' + 
                        '"true_mass" or "SZ"')
        
    cluster_z = Nc.ClusterRedshift.new_from_name("NcClusterPhotozGaussGlobal{'pz-min':<%f>, 'pz-max':<%f>, 'z-bias':<0.0>, 'sigma0':<%f>}" % (z_min, z_max, z_sigma))
    cad = Nc.ClusterAbundance.new(mf, None, cluster_z, cluster_m)

    self.ncdata = Nc.DataClusterNCount.new(cad)
    self.mset = Ncm.MSet()
    self.mset.set(self.cosmo)
    self.mset.set(cluster_m)

    self.rng = Ncm.RNG.pool_get("example_ca_sampling");
    self.ncdata.init_from_sampling(self.mset, cluster_z, cluster_m, 
                                   area * (pi/180.0) ** 2, self.rng)

    del dist
    del vp
    del gf
    del mulf
    del mf
    del cad
    del cluster_z
    del cluster_m
    
  def simulation(self, z_max, CP, seed):
    self.cosmo.props.H0      = CP["H0"]
    self.cosmo.props.Omegab  = CP["Ob"]
    self.cosmo.props.Omegac  = CP["Om"]
    self.cosmo.props.Omegax  = CP["OL"]
    self.cosmo.props.Tgamma0 = CP["Tgamma"]
    self.cosmo.props.ns      = CP["ns"]
    self.cosmo.props.sigma8  = CP["sigma8"]
    self.cosmo.props.w       = CP["w"]
    

    if seed == False:
        self.rng.set_random_seed(False)
    else:   
        Ncm.RNG.set_seed(self.rng, seed)
    
    self.ncdata.resample(self.mset, self.rng)

    lnM_true = self.ncdata.get_lnM_true()
    z_true = self.ncdata.get_z_true()
    lnM_obs = self.ncdata.get_lnM_obs()
    z_obs = self.ncdata.get_z_obs()

    nobjects = self.ncdata.get_len()

    header = ['z_true', 'Mass_true', 'z_obs', 'Mass_obs']

    sim = []
    for i in range(nobjects):
        sim.append([z_obs.get (i, 0), lnM_obs.get (i, 0)])

    del lnM_true
    del z_true
    del lnM_obs
    del z_obs
    
    return header, sim
    
class ChooseParamsInput(object):
    params=None
    keys=None
    keys_values=None
    keys_cov=None
    keys_bounds=None
    sdata=None
    sdata_weights=None
    prior_dist=None
    
    def set_de(self):
       if "Om" in self.keys:
         self.params["OL"] = 1. - self.params["Om"] - self.params["Ob"]
    
    def update_keys(self,x):
       for i in xrange(len( self.keys )):
          self.params[self.keys[i]] = x[i]



def numcosmo_sim_cluster(simul_params, save=False):
        """
        Perform simulation using NumCosmo library given a set of input 
        cosmological parameters.

        input:	simul_params -> NumCosmo object containing 
                                cosmological parameters

		save (optional) -> boolean (save simulation in 
                                   output file. Default is False)	

        output: set of 2-dimensional arrays of simulated catalog:
		collumns -> [ redshift, observable]
        """    
 
        #prepara simulation object
        ncount=NCountSimul(simul_params['zmin'], simul_params['zmax'], 
                           simul_params['zsigma'], simul_params['MinMass'], 
                           simul_params['MaxMass'], simul_params['area'], 
                           simul_params['observable'])

        data_simul = []
        
        #generate simulation
        data_simul = np.array(ncount.simulation(simul_params['zmax'], 
                              simul_params, simul_params['seed'])[1])
        
       
        if len(data_simul) == 0:
            data_simul = np.array([[0,0]])

       
        if save == True:
            op1 = open('simulation.dat', 'w')
            for line in data_simul:
                for item in line:
                    op1.write(str( item ) + '    ')
                op1.write('\n')
            op1.close()
     
        return data_simul

############################################################################


def main():
  print(__doc__)

if __name__=='__main__':
  main()      
