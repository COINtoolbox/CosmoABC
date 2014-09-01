#!/usr/bin/python
from gi.repository import NumCosmo as Nc
from gi.repository import NumCosmoMath as Ncm
from numpy import pi
import numpy as np
import random
import math





class NCountSimul:

  def __init__ (self, z_min, z_max, lnM_min, lnM_max, area,seed):
    Ncm.cfg_init ()
    self.cosmo = Nc.HICosmo.new_from_name (Nc.HICosmo, "NcHICosmoDEXcdm")
    dist = Nc.Distance.new (z_max * 1.5)
    wp =  Nc.Window.new_from_name ("NcWindowTophat")
    tf = Nc.TransferFunc.new_from_name ("NcTransferFuncEH")
    vp = Nc.MatterVar.new (Nc.MatterVarStrategy.FFT, wp, tf)
    gf = Nc.GrowthFunc.new ()
    mulf = Nc.MultiplicityFunc.new_from_name ("NcMultiplicityFuncTinkerCrit{'Delta':<500.0>}")
    mf = Nc.MassFunction.new (dist, vp, gf, mulf)
    
    cluster_m = Nc.ClusterMass.new_from_name ("NcClusterMassBenson{'M0':<3e14>, 'z0':<0.6>, 'signif-obs-min':<5.0>, 'Asz':<6.24>, 'Bsz':<1.33>, 'Csz':<0.83>, 'Dsz':<0.24>}")
    cluster_z = Nc.ClusterRedshift.new_from_name ("NcClusterPhotozGaussGlobal{'pz-min':<%f>, 'pz-max':<%f>, 'z-bias':<0.0>, 'sigma0':<0.05>}" % (z_min, z_max))
    cad = Nc.ClusterAbundance.new (mf, None, cluster_z, cluster_m)

    self.ncdata = Nc.DataClusterNCount.new (cad)

    self.mset = Ncm.MSet ()
    self.mset.set (self.cosmo)
    self.mset.set (cluster_m)

    self.rng = Ncm.RNG.pool_get ("example_ca_sampling");

    self.ncdata.init_from_sampling (self.mset, cluster_z, cluster_m, area * (pi / 180.0)**2, self.rng)
    
    if seed == False:
        self.rng.set_random_seed( False )
    else:   
        Ncm.RNG.set_seed ( self.rng , seed )
    self.ncdata.resample (self.mset, self.rng)

    del dist
    del vp
    del gf
    del mulf
    del mf
    del cad
    del cluster_z
    del cluster_m

  def simulation (self, z_max,  CP):
    self.cosmo.props.H0      = CP["H0"]
    self.cosmo.props.Omegab  = CP["Ob"]
    self.cosmo.props.Omegac  = CP["Om"]
    self.cosmo.props.Omegax  = CP["OL"]
    self.cosmo.props.Tgamma0 = CP["Tgamma"]
    self.cosmo.props.ns      = CP["ns"]
    self.cosmo.props.sigma8  = CP["sigma8"]
    self.cosmo.props.w       = CP["w"]        

    

    lnM_true = self.ncdata.get_lnM_true ()
    z_true = self.ncdata.get_z_true ()
    lnM_obs = self.ncdata.get_lnM_obs ()
    z_obs = self.ncdata.get_z_obs ()

    nobjects = self.ncdata.get_len ()

    header = ['z_true', 'Mass_true', 'z_obs', 'Mass_obs']

    sim = []
    for i in range (nobjects):
        #sim.append ([ z_true.get(i), lnM_true.get (i), z_obs.get (i, 0), lnM_obs.get (i, 0)])
        sim.append ([ z_obs.get (i, 0), lnM_obs.get (i, 0)])

    del lnM_true
    del z_true
    del lnM_obs
    del z_obs
    
    return header, np.array( sim ) 


