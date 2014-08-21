from gi.repository import NumCosmo as Nc
from gi.repository import NumCosmoMath as Ncm
from numpy import pi
import numpy
import random

from scipy.stats import anderson_ksamp
from scipy.stats.mstats import mquantiles

class NCountSimul:

  def __init__ (self, z_min, z_max, lnM_min, lnM_max, area):
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

    del dist
    del vp
    del gf
    del mulf
    del mf
    del cad
    del cluster_z
    del cluster_m

  def simulation (self, z_max, seed, CP):
    self.cosmo.props.H0      = CP["H0"]
    self.cosmo.props.Omegab  = CP["Ob"]
    self.cosmo.props.Omegac  = CP["Om"]
    self.cosmo.props.Omegax  = CP["OL"]
    self.cosmo.props.Tgamma0 = CP["Tgamma"]
    self.cosmo.props.ns      = CP["ns"]
    self.cosmo.props.sigma8  = CP["sigma8"]
    self.cosmo.props.w       = CP["w"]        

    Ncm.RNG.set_seed ( self.rng , seed )
    self.ncdata.resample (self.mset, self.rng)

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
    
    return header, sim 

def summary_anderson( fid_data, sim_data, mass_bin ):


    #convert mass bins to natural log
    mass = [ numpy.log( item ) for item in mass_bin ]

    #separate fiducial data
    fid_data_bin = numpy.array([ [ item[0] for item in fid_data  if item[1] >= mass[ i ] ] for i in range (len(mass) -1) ])
    
    #separate simulated data
    sim_data_bin = numpy.array([ [ item[0] for item in sim_data  if item[1] >= mass[ i ] ] for i in range (len(mass) -1) ])

    for ii in range( len( fid_data_bin ) ):
        print 'sim = ' + str( len( sim_data_bin[ ii ] ) ) + '   fid = ' + str( len( fid_data_bin[ ii ] ) )

    #perform anderson darling test
    res = numpy.array( [ [ anderson_ksamp( [ fid_data_bin[ i ], sim_data_bin[ i ] ] )[ j ] for j in [0,2] ]  for i in range( len( fid_data_bin ) ) if ( len( fid_data_bin[ i ] ) > 0 and len( sim_data_bin[ i ] ) > 0 ) ])

    return res 


def summary_quantile( data, mass_bin, q_list ):
    #Summary statistics by as the quantiles in q_list.
    #Calculates the quantiles for a distribution of redshift for observed masses above a mass threshold given by mass_bin.
    #
    #input: data -> list of list of floats 
    #               [ redshit, log_obs_Mass ]
    #       mass_bin -> list of floats
    #               [ mass_threshold_in_base_10 ]   
    #       q_list -> list of float
    #               [ quantiles ]
    #
    #output: res -> list of list of float
    #               [ calculated_quantiles_for_each_mass_threshold ]
    #        pop -> list of float
    #               [ fraction_of_number_of_data_for_each_mass_threshold ]    


    #convert mass bins to natural log
    mass = [ numpy.log( item ) for item in mass_bin ]

    #separate fiducial data
    data_bin = numpy.array([ [ item[0] for item in data  if item[1] >= mass[ i ] ] for i in range (len(mass) -1) ])

    #for item in data_bin:
    #    print ' bin = ' + str( len( item ) )
  
    #calculate quantile 
    res = [ mquantiles( elem, prob=q_list ) if len( elem ) > 0  else [ 0 for jj in q_list]  for elem in data_bin ] 
 
    if sum( [ len( data_bin[ k ] ) for k in range( len( data_bin ) ) ] ) > 0:    
        pop = [ float( len( data_bin[ i ] ) )/sum( [ len( data_bin[ k ] ) for k in range( len( data_bin ) ) ] )  for i in range( len( data_bin ) ) ]   
    else:
        pop = [ 0 for i in range( len( data_bin ) ) ] 

    return res, pop

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
         self.params["OL"]=1.-self.params["Om"]-self.params["Ob"]
    
    def update_keys(self,x):
       for i in xrange( len( self.keys ) ):
          self.params[ self.keys[ i ] ]=x[ i ]


def deviation_quantile( summary_fid, summary_sim ):
    #Calculates the distance between 2 data sets given the quantiles calculated by summary_quantile function.
    #Distances are calculated for each mass threshold, e.g. mass1, as
    #    sqrt( sum_in_i( ( fraction_data_mass1_fid * quantile[ i ]_mass1_fid - fraction_data_mass1_sim * quantile[ i ]_mass1_sim ) ** 2 ))
    #
    #
    #input: summmary_fid, summary_sim -> list of list of float 
    #                                    [ calculated_quantiles_for_each_mass_threshold ] 
    #                                    outputs from summary_quantile for 2 distinct data sets
    #
    #output: distance -> list of float
    #                    [ distances_for_each_mass_bin ]             

  
    distance = [ numpy.sqrt( sum( [ ( summary_fid[1][i] * summary_fid[0][ i ][ j ] - summary_sim[1][i] * summary_sim[0][ i ][ j ] ) ** 2  for j in range( len( summary_fid[0][ i ] ) ) ] ) ) for i in range( len( summary_fid[0] ) ) ]
    
    return distance

def choose_par( hyper_par, dist, Ob ):
    """
    Sample cosmological parameter from prior. 
    Current constraints are:
        0 <= om < Ob
        0 <= om + ol + Ob < 1
        w <= 0

    input:
        hyper_par : list of list of float
            hyper_par[0][0], hyper_par[0][1] : om1, om2 -> float, float: parameters that describe prior for matter energy density
            hyper_par[1][0], hyper_par[1][1] : w1, w2   -> float, float: that describe prior for dark energy equation of state
            hyper_par[2][0], hyper_par[2][1] : sig81, sig82 -> float, float:parameters that describe prior for sigma8
        dist     -> string: prior distribution:  'normal' -> gaussian 
                                                 'flat'   -> top-hat
        Ob :  float
            baryon energy density parameter

    output:
        om_try -> float: dark matter density 
        w_try  -> float: equation of state parameter
        sig8_try -> float: sigma8
    """


    om_try = 0.0
    w_try  = 0.0
    sig8_try = 0.0

    flag = False
    while flag == False:
        
        if dist == 'normal':
            om_try = random.normalvariate( hyper_par[0][0], hyper_par[0][1] )
            w_try = random.normalvariate( hyper_par[1][0], hyper_par[1][1] ) 
            sig8_try  = random.normalvariate( hyper_par[2][0], hyper_par[2][1] )
            

        elif dist == 'flat':
            om_try = random.uniform(  hyper_par[0][0] - hyper_par[0][1], hyper_par[0][0] + hyper_par[0][1] )
            w_try = random.uniform(  hyper_par[1][0] - hyper_par[1][1], hyper_par[1][0] + hyper_par[1][1] )
            sig8_try  = random.uniform(  hyper_par[2][0] - hyper_par[2][1],  hyper_par[2][0] + hyper_par[2][1] ) 

        if om_try <= Ob or om_try > 1 - Ob or w_try > 0.0 or sig8_try > 1.0 or sig8_try < 0.2:
            flag = False 
        else:
            flag = True


    return om_try, w_try, sig8_try



def set_distances( summary_fid, mass_bin, quantile_list, omX, wX, sigX, zmin, zmax, area, ncount1, H0=71.15, Omegab=0.0439, Tgamma0=2.725, ns=0.97 ):
    """
    Calculate summary statistics difference between the fiducial data and a specific model defined by the inputed cosmological parameters.

    input: 
           summary_fid: list of float outputed by function summary_quantile
                        [ fraction_of_number_of_data_for_each_mass_threshold ]    

           mass_bin: list of floats
                     [ mass_threshold_in_base_10 ]  
 
           quantile_list: list of float
                          [ quantiles ]
     
           omX: float
                model dark matter density
 
           olX: float
                model dark energy density

           wX:  float
                model dark energy equation of state parameter

    output:  
           difference: float
                       summary statistic difference between the fiducial data and inputed cosmological model as outputed by function deviation_quantile 
        
    """

    #simulate instance of data
    data_sim = numpy.array( ncount1.simulation( zmax, H0, Omegab, omX, 1-omX-Omegab, Tgamma0, ns, sigX, wX )[1] )

    #calculate summary statistics for simulated data
    summ_sim = summary_quantile( data_sim, mass_bin, quantile_list )
    
    #calculate deviation for every mass threshold
    difference = deviation_quantile( summary_fid, summ_sim ) 

    del summ_sim
    del data_sim
  
    if difference == 0.0:
        return 1000
    else:
        return difference

def choose_surv_par( summary_fid, mass_bin, quantile_list, tolerance, n_tries, hyper_par, dist, zmin, zmax, area, ncount1, Ob  ):
    """
    Select model parameters surviving summary statistics distance threshold.

    input:
           summary_fid: list of float outputed by function summary_quantile
                        [ fraction_of_number_of_data_for_each_mass_threshold ]    

           mass_bin: list of floats
                     [ mass_threshold_in_base_10 ]  
 
           quantile_list: list of float
                          [ quantiles ]

           tolerance: float
                      maximum threshold of summary statistic distance. Parameter values are kept if the calculated model distance is lower.

           n_tries: int
                    number of surviving parameter values

           hyper_par : list of list of float
            hyper_par[0][0], hyper_par[0][1] : om1, om2 -> float, float: parameters that describe prior for matter energy density
            hyper_par[2][0], hyper_par[1][1] : w1, w2   -> float, float: that describe prior for dark energy equation of state
            hyper_par[1][0], hyper_par[2][1] : sig81, sig82 -> float, float:parameters that describe prior for sigma8

           dist: string
                 prior distribution:  'normal' -> gaussian 
                                      'flat'   -> top-hat
           Ob :  float
            baryon energy density parameter 

    output:
           result: list of list of float
                   parameter and distance surviving threshold requirement
                   [[ dark_matter_density, dark_energ_density, equation_of_state_parameter, distance ]]  

           data: list of list of float   
                 simulated data for surviving model 
                 [[ redshift, mass]]
    """

    #list to store results
    result = []
    data_all = []

    while len( result ) < n_tries :

        print 'tries = ' + str( len( result ) )

        #choose model parameters
        par_list = list( choose_par( hyper_par, dist, Ob )  )

        print par_list

        #calculate summary statistics distance to fiducial data  
        d = set_distances( summary_fid, mass_bin, quantile_list, par_list[0], par_list[1], par_list[2], zmin, zmax, area, ncount1  )
        d1 = sum( d )

        if d1 <=  tolerance:

            
            print '        dist = ' + str( d1 )

            #add distance to parameter list
            par_list.append( d1 )

            #append parameters if tolerance is satisfied
            result.append( par_list  )
    
    return result

