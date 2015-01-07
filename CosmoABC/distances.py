"""
Functions for distance determination.

***
In order to insert a new distance function, add it to this file.
***
""" 


import numpy
from scipy.interpolate import Rbf

def SumGRBF( dataset1, dataset2, s1, f1='gaussian'):
    """
    Compute the sum of Gaussian Radial Basis Function (GRBF) for all points in data1 having data2 as a fixed reference basis.
    Following equation 18 notation of Ishida et al., 2015: \mathcal{D} -> data1 and Y-> data2.  

    input:	data1 -> first data matrix ( collumns -> features, lines -> objects )
   		data2 -> second data matrix, fixed basis (collumns -> features, lines -> objects )
 		s -> epsilon in scipy notation (scalar)
                f1 (optional) -> kernel function implemented in scipy.interpolate.Rbf rotine (string)

    output: scalar
    """

    #GRBF must be one in all points in the fixed basis sample
    one = [ 1 for j in range( len( dataset2 ) )  ]

    rbf1 = Rbf( dataset2[:,0], dataset2[:,1], one, function=f1, epsilon=s1 )
    
    #sum GRBF for all elements in data1 and centralize in the number of existing points in data2
    sum1 = abs( numpy.nansum( [ numpy.log( rbf1( line[0], line[1] ) ) for line in dataset1 ]  ) - len( dataset2 ))
    
    return sum1 
        


def distance_GRBF( dataset1, dataset2, s1, f='gaussian' ):
    """
    Compute the distance between two samples (equation 19).
    Following equation 19 notation of Ishida et al., 2015: \mathcal{D} -> data1 and Y-> data2.  

    input:	data1 -> first data matrix ( collumns -> features, lines -> objects )
   		data2 -> second data matrix, fixed basis (collumns -> features, lines -> objects )
 		s -> scale parameter (scalar)
                f (optional) -> kernel function implemented for GRBF rotine (string)
              
    output: scalar
    """

    if sum( dataset2[0] ) == 0:
        return 10**10
     
  
    else:
        #distance based on the logarithm of GRBF and normalized by results from data1
        if SumGRBF( dataset1, dataset2, s1 ) > 0:
            d = abs( -2 * self.SumGRBF( dataset1, dataset2, s1, f1=f ) + 2 * self.SumGRBF( dataset1, dataset1, s1, f1=f ))
            return d

        else:
            op1 = open( 'dataset2_error_source.dat', 'w' )
            for line in dataset2:
                for item in line:
                    op1.write( str(item) + '    ' )
                op1.write('\n')
            op1.close() 
               
            raise ValueError('ERROR in function SumGRBF!!  Corresponding simulation is stored in file "dataset2_error.dat" ' )



def main():
  print(__doc__)

if __name__=='__main__':
  main()


