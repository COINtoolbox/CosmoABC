"""
Functions for distance determination.
Created by Emille Ishida in 2015.
""" 

import numpy as np
import sys
from scipy.stats.mstats import mquantiles


###############################################################
### GRBF distance    -  Appendix B from Ishida et al., 2015

def GRBF(vec):
    """
    Gaussian Radial Basis Function (GRBF). 
    Equation B.1 of Ishida et al., 2015.

    input:  vec -> vector of input parameters
                   vec[0] -> array of observabed values for 1 object
                   vec[1] -> dictionary of parameters

    output: scalar -> GRBF at obs values
    """

    try:
        core_exp = 0
        for obj in vec[1]['sim']:
            d = vec[0] - obj
            term  = np.exp(-(np.dot(np.dot(d, vec[1]['inv_cov']), d.T))/2.0)

            core_exp = core_exp + vec[1]['const'] * term
  
        return np.log(core_exp)

    except KeyboardInterrupt:
        pass

def logf(var):
    """
    Log of GRBF between 2 data sets. 
    Equation B.3 of Ishida et al, 2015.

    input: var -> dictionary of parameters

    output: scalar  
    """
    args = [[line, var] for line in var['dataset1']]

    result = [GRBF(item) for item in args]

    output = sum(result) - len(var['sim'])

    return output


def norm_GRBF(var):
    """
    GRBF relating the observed data to itself. 
    
    input: var -> dictionary of parameters

    output: scalar -> normalization constant
    """
    var['sim'] = var['dataset1']

    return logf(var)


def prep_GRBF(var):
    """
    Calculates all parameters needed by the GRBF definition. 

    input: var -> dictionary of parameters

    output: var -> updated dictionary of parameters enabling 
                   GRBF calculations
    """

    var['cov'] = (var['s'] ** 2) * np.cov(var['dataset1'].T)
    var['inv_cov'] = np.linalg.inv(var['cov'])
    det = np.exp(np.linalg.slogdet(var['cov'])[1])
    var['const'] = 1.0/(2 * np.pi * np.sqrt(det))

    var['norm_GRBF'] = norm_GRBF(var)

    return var

def distance_GRBF(dataset2, var):
    """
    Distance based on GRBF as defined in equation B.4 of Ishida et al, 2015.

    input: dataset2 -> array of simulated catalogue
 
    output: scalar -> GRBF distance
    """

    var['sim'] = dataset2

    return [-2 * logf(var) + 2 * var['norm_GRBF']]

###############################################################
### Quantile distances - section 4.2 of Ishida et al., 2015

def summ_quantiles( dataset1, Parameters ):
    """
    Updates the dictionary of input parameters with the summary based on 
    quantiles for the observed catalogue. 
    ** This is an auxiliary function in order to optimize parallelization**

    input: dataset1 -> observed catalogue
   
    output: Parameters -> updated dictionary of input parameters where the 
                          keyword 'extra' correspond to the summary calculated 
                          for the observed data set.
    """

    
    Parameters['dist_dim'] = len(dataset1[0]) + 1
 
    qlist = np.arange(0.05, 1.0, 0.95/Parameters['quantile_nodes'])

    Parameters['extra'] = []
    for i1 in range(len( dataset1[0])):

        Parameters['extra'].append(np.array([mquantiles(dataset1[:,i1], 
                                   prob=item) for item in qlist]))
   
    return Parameters


def distance_quantiles(dataset2, Parameters):
    """
    Compute the distance between two samples based on quantiles of 
    observed features.

    input: dataset2 -> simulated catalog
           Parameters -> dictionary of input parameters	
              
    output: distance -> scalar
    """


    qlist = np.arange(0.05, 1.0, 0.95/Parameters['quantile_nodes'])


    qd = []
    for j1 in range(len(dataset2[0])):
        qd.append(np.array([mquantiles(dataset2[:, j1 ], prob=item) 
                              for item in qlist]))
        

    l1 = len(Parameters['dataset1'])
    l2 = len(dataset2)

    d = []
    for j2 in range(len(dataset2[0])):
        d.append(np.linalg.norm(Parameters['extra'][ j2] - qd[j2]))

    
    d.append(max(abs(1-float(l1)/l2), abs(1-float(l2)/l1)))

    return d

###############################################################

def main():
  print(__doc__)

if __name__=='__main__':
  main()


