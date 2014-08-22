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

    if seed == False:
        self.rng.set_random_seed( False )
    else:   
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

  
    #distance = [ numpy.sqrt( sum( [ ( summary_fid[1][i] * summary_fid[0][ i ][ j ] - summary_sim[1][i] * summary_sim[0][ i ][ j ] ) ** 2  for j in range( len( summary_fid[0][ i ] ) ) ] ) ) for i in range( len( summary_fid[0] ) ) ]
    
    print '****no weight in distance **' 
    distance = [ numpy.sqrt( sum( [ ( summary_fid[0][ i ][ j ] -  summary_sim[0][ i ][ j ] ) ** 2  for j in range( len( summary_fid[0][ i ] ) ) ] ) ) for i in range( len( summary_fid[0] ) ) ]
    
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
    Nparams=len(hyper_par)
    
   # U,S,V=linalg.svd(hyper_par_cov)
    
   # new_mean=numpy.dot(numpy.transpose(U),hyper_par)
    

    flag = numpy.array([ False for i in xrange( Nparams ) ])
    
    print "***********************************"
    while (flag.all() == False ):
          L = numpy.linalg.cholesky( hyper_par_cov )
          norm = numpy.random.normal(size=1*Nparams).reshape(Nparams, 1)
          params =numpy.array( hyper_par) + numpy.dot(L, norm).reshape(1,Nparams)
          #print norm,hyper_par,params
          flag=((params[0]>=bounds[0]) & (params[0]<=bounds[1]))
         # print params,flag
    #print "***********************************"
    #print "generated parameters=", numpy.array( params[0] )
    return numpy.array(params[0])
    


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
    data_size = len( data_sim )

    #calculate summary statistics for simulated data
    summ_sim = summary_quantile( data_sim, mass_bin, quantile_list )
    
    #calculate deviation for every mass threshold
    difference = deviation_quantile( summary_fid, summ_sim ) 

    del summ_sim
    del data_sim 
  
    if difference == 0.0:
        return 1000
    else:
 
        ##################
        # return the dimension of simulated data set   
        return difference, data_size


def weighted_values(values, probabilities, size):
    bins = numpy.cumsum(probabilities)
    return values[numpy.digitize(numpy.random.random_sample(size), bins)]
    
    

def choose_surv_par( summary_fid, mass_bin, quantile_list, tolerance, n_tries, CP, zmin, zmax, area, ncount, seed, ndata_fid ):
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

    
    
    Nparams=len(CP.keys)
    
    if ( CP.sdata is None ):
       
       flag = True 
       
       par_mean = CP.keys_values[:]
       
       par_surv = par_mean[:]
       
       par_cov =  CP.keys_cov[:] 
      
    else:
       
       flag = False 
       
       par_mean = CP.sdata[:,:Nparams].mean(axis=0) # compute the mean of the previous data
       
       par_surv = CP.sdata[:,:Nparams]-par_mean[:]  # center previous data
       
    #   print "###################"
    #   print par_surv
    #   print "###################"
       
       weights = CP.sdata_weights[:]
       
     
       par_cov = 2*numpy.dot( numpy.transpose( par_surv ), numpy.dot( numpy.diag( weights ) , par_surv ) ) # compute new covariance matrix
       
       print par_cov

   # 2*numpy.dot( numpy.transpose( par_surv[:,:Nparams] ),numpy.dot( numpy.diag( weights ) , par_surv[:,:Nparams] ) )

    bounds = CP.keys_bounds[:]

    indx_list= []
    
    result = []
    
    par_list = []
        
    for k in xrange( n_tries ) : 
     
        d1 = 10*tolerance[0]
        d2 = 10*tolerance[1]
       # print "d1=",d1,"tolerance=",tolerance
        while ( d1 >= tolerance[0] ) or d2 >= tolerance[1]:
        
        
           #########################################################
           #draw from a multivariate normal distribution

           if ( flag == True ):
              
              indx=numpy.array( [ 0 ] )
              
              new_par = choose_par( par_surv, par_cov, bounds, CP.prior_dist ) 
           
           else:
              
              #choose one of the instances of simulation set
              indx = weighted_values(range( n_tries ), weights,1)
              print "Choosed iten=",indx,CP.sdata[indx,:Nparams]
              new_par = choose_par( CP.sdata[indx,:Nparams], par_cov, bounds, CP.prior_dist ) 
           
           #print new_par
           
           CP.update_keys( new_par ) # update dictionary with containing the cosmological parameters
        
           CP.set_de()  # update OL key
        
           #print CP.params
           
           d = set_distances( summary_fid, mass_bin, quantile_list, zmin, zmax, area, ncount, seed, CP.params  )

           ##################
           # set dimension of output from set_distances
           d1 = sum( d[0] )

           
           if ( d1 <=  tolerance[0] ) and d[1] > 0:

               d2 = max( abs( 1 - float( d[1] )/float( ndata_fid ) ), abs( 1 - float( ndata_fid  )/float( d[1] ) ) )
 
               if  d2 <= tolerance[1]:

                    print 'tries = ' + str( len( par_list ) )

                    indx_list.append( indx )
                    
                    print '        dist = ' + str( d1 ) + ', epsilon1 = ' + str( tolerance[0] ) + ',   epsilon2 = ' + str( tolerance[1] )

                    result = list( new_par ) + [ d1, d2 ]   
                    print "Accpeted point:", result

                    #add distance to parameter list
                    par_list.append( result )

            ###############################################################
 
     
     
    del par_surv
    del new_par
    del par_mean
    del d1
    del flag
    del bounds
       
   
    return numpy.array( par_list ),numpy.array( indx_list ),par_cov



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




