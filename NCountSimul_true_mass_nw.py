from gi.repository import NumCosmo as Nc
from gi.repository import NumCosmoMath as Ncm
from numpy import pi
import numpy
import random
import math

from scipy.stats import anderson_ksamp
from scipy.stats.mstats import mquantiles
from numpy import linalg 


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
    
    cluster_m = Nc.ClusterMass.new_from_name ("NcClusterMassNodist{'lnM-min':<% 20.15g>, 'lnM-max':<% 20.15g>}" % (lnM_min, lnM_max))
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



def choose_par( hyper_par, hyper_par_cov, bounds,dist ):
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
    size=len(hyper_par)
    
   # U,S,V=linalg.svd(hyper_par_cov)
    
   # new_mean=numpy.dot(numpy.transpose(U),hyper_par)
    

    flag=numpy.array([ False for i in xrange(size) ])
    
    while (flag.all() == False ):
    
 
        
        if dist == 'normal':
           # y=numpy.random.multivariate_normal(new_mean,numpy.diag(S))
            
           # params=numpy.dot(U,y)
            params=numpy.random.multivariate_normal(hyper_par,hyper_par_cov)
            
            flag=((params>=bounds[0]) & (params<=bounds[1]))

    return params
    


def set_distances( summary_fid, mass_bin, quantile_list, zmin, zmax, area, ncount1, seed, CP):
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
    
    data_sim = numpy.array( ncount1.simulation( zmax, seed, CP )[1] )

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



def choose_surv_par( summary_fid, mass_bin, quantile_list, tolerance, n_tries, CP, keys, hyper_par_cov, bounds,dist, zmin, zmax, area, ncount1, seed ):
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
    
    
    hyper_par = numpy.array([CP[iten] for iten in keys])
    
    
    result = []
    data_all = []

    while len( result ) < n_tries :

        print 'tries = ' + str( len( result ) )

        #choose model parameters
        par_list = list( choose_par( hyper_par,hyper_par_cov,bounds, dist )  )
        
        for i in xrange( len( keys ) ):
           CP[keys[i]]=par_list[i]

        if "Om" in keys:
           CP["OL"]=1-CP["Om"]
        
        print par_list

        #calculate summary statistics distance to fiducial data  
        d = set_distances( summary_fid, mass_bin, quantile_list,  zmin, zmax, area, ncount1, seed, CP  )
        d1 = sum( d )

        if d1 <=  tolerance:

            
            print '        dist = ' + str( d1 )

            #add distance to parameter list
            par_list.append( d1 )

            #append parameters if tolerance is satisfied
            result.append( par_list  )
    
    return numpy.array(result)


def weighted_values(values, probabilities, size):
    bins = numpy.cumsum(probabilities)
    return values[numpy.digitize(numpy.random.random_sample(size), bins)]


def norm_pdf_multivariate(x, mu, sigma):
  size = len(x)
  if size == len(mu) and (size, size) == sigma.shape:
    det = linalg.det(sigma)
    if det == 0:
        raise NameError("The covariance matrix can't be singular")

    norm_const = 1.0/ ( math.pow((2*numpy.pi),float(size)/2) * math.pow(det,1.0/2) )
    x_mu = numpy.matrix(x - mu)
    inv = numpy.matrix(sigma).I        
    result = math.pow(math.e, -0.5 * (x_mu * inv * x_mu.T))
    return norm_const * result
  else:
    raise NameError("The dimensions of the input don't match")




