"""
Functions for distance determination.

***
In order to insert a new distance function, add it to this file.
***
""" 


import numpy as np

from scipy.interpolate import Rbf
from scipy.stats.mstats import mquantiles

def SumGRBF(dataset1, dataset2, Parameters):
    """
    Compute the sum of Gaussian Radial Basis Function (GRBF).  

    input: data1 -> first data matrix ( collumns -> features, lines -> objects )
   	   data2 -> second data matrix, fixed basis (collumns -> features, 
                                                     lines -> objects )
 	   Parameters -> dictionary of input parameters

    output: distance -> scalar
    """

    #GRBF must be one in all points in the fixed basis sample
    one = [1 for j in range(len(dataset2))]

    if len(dataset2[0]) == 1:

        rbf1 = Rbf(dataset2,  one, function=Parameters['kernel_func'], 
                   smooth=Parameters['s'])
    
        #sum GRBF for all elements in data1 and centralize in the number of 
        #existing points in data2
        sum1 = abs(np.nansum([np.log(rbf1(line[0])) for line in dataset1 ]) 
                   - len( dataset2 ))
    
    if len(dataset2[0]) == 2:

        rbf1 = Rbf(dataset2[:,0], dataset2[:,1], one, 
                   function=Parameters['kernel_func'], smooth=Parameters['s'])
    
        #sum GRBF for all elements in data1 and centralize in the number of 
        #existing points in data2
        sum1 = abs(np.nansum([np.log(rbf1(line[0], line[1])) 
                   for line in dataset1]) - len(dataset2))
           
    return sum1 
        


def distance_grbf(dataset2, Parameters):
    """
    Compute the distance between two catalogues based on Gaussian Radial 
    Basis function.
    
    input: dataset2 -> simulated catalogue
           Parameters -> dictionary of input parameters
              
    output: distance -> scalar
    """

    if sum(dataset2[0]) == 0 or len(dataset2) > 5*len(Parameters['dataset1']) or dataset2.shape[0] == 1:
        return 10**10
    else:
 
        j1 = SumGRBF(Parameters['dataset1'], dataset2, Parameters)

        #distance based on the logarithm of GRBF and normalized by 
        #results from data1
        if j1 > 0 and str( j1 ) != 'nan' :
            d = abs( -2 * j1 + 2 * Parameters['extra'])
            return [d]

        else:
            op1 = open('dataset2_error_source.dat', 'w')
            for line in dataset2:
                for item in line:
                    op1.write(str(item) + '    ')
                op1.write('\n')
            op1.close() 
               
            raise ValueError('ERROR in function SumGRBF!!  Corresponding simulation is stored in file "dataset2_error.dat"')



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
    for i1 in xrange(len( dataset1[0])):

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


    qlist = np.arange( 0.05, 1.0, 0.95/Parameters['quantile_nodes'] )


    qd = []
    for j1 in xrange( len( dataset2[0] ) ):
        qd.append( np.array([ mquantiles( dataset2[:, j1 ], prob=item ) 
                              for item in qlist ]) )
        

    l1 = len( Parameters['dataset1'] )
    l2 = len( dataset2 )

    d = []
    for j2 in xrange( len( dataset2[0] ) ):
        d.append(np.linalg.norm(Parameters['extra'][ j2] - qd[j2]))

    
    d.append(max(abs(1-float(l1)/l2), abs(1-float(l2)/l1)))

    return d


def main():
  print(__doc__)

if __name__=='__main__':
  main()


