"""
Functions for plotting.
"""


import datetime
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
from weighted_gaussian_kde import *
import matplotlib.gridspec as gridspec
import numpy


def plot_1D( T, file_output, Parameters):
    """
    Make 1-dimensional plot for ABC results. 

    input:  T 		-> number of total particle systems (int)  
	    file_output -> name of output figure file (str)
            Parameters 	-> input parameter dictionary - same as used for ABC sampler (dict)

    output: summary pdf file with multiple pages for each generated particle system. 

    *******
    If you want to change the range of parameters for plotting, change values for the key 'param_lim' in the input dictionary parameter! 
    *******
    """

    if file_output[-3:] != 'pdf':
        raise NameError( 'Name file for figure output must be a pdf!' )       


    sampling = numpy.array( [ [ i for i in numpy.arange( Parameters['param_lim'][ j ][0], Parameters['param_lim'][ j ][1], (Parameters['param_lim'][ j ][1]-Parameters['param_lim'][ j ][0])/1000) ] for j in range( len( Parameters['param_to_fit'] ) )] )

    with PdfPages( file_output ) as pdf:

        epsilon_ev = []
        time_ev = []
        ndraws_ev = []

        for i in range( T ):
           
            op1 = open( Parameters['file_root'] + str( i ) + '.dat', 'r' )
            lin1 = op1.readlines()
            op1.close()

            d = [ elem.split() for elem in lin1 ]
            d1 = numpy.array([ [ float( item ) for item in line ] for line in d[1:] ])
     

            epsilon_ev.append( float( d[1][ d[0].index( 'dist_threshold' ) ] ) )
            time_ev.append( sum( float( line[ d[0].index('time') ] ) for line in d1[1:]) )
            ndraws_ev.append( sum( float( line[ d[0].index('NDraws') ] ) for line in d1[1:]) )

            if i > 0:
                w1 = numpy.loadtxt( Parameters['file_root'] + str( i ) + 'weights.dat' )
            else:
                w1 = numpy.array([ 1.0/len(d1) for k in range( len(d1) ) ])


            kde1 = gaussian_kde(d1[:,0] , weights=w1)
            y1 = kde1( sampling[0] )

            #### Plot posteriors
            plt.figure()
            plt.plot( sampling[0], y1, color='blue')
            plt.xlabel( Parameters['param_to_fit'][0] )
            plt.ylabel( 'density', fontsize=8 )
            plt.tick_params(axis='both', which='major', labelsize=8)
            plt.xlim( Parameters['param_lim'][0][0]-0.1, Parameters['param_lim'][0][1]+0.1)
            pdf.savefig()
            plt.close()

    
        convergence = numpy.array( [ float(Parameters['M'])/item  for item in ndraws_ev ] )
    
        #Plot epsilon evolution
        plt.figure()
        plt.scatter( range( T ), numpy.array(epsilon_ev), color='red', marker='*', label='distance threshold' )
        plt.legend()
        plt.xlabel( 'Particle System' )
        plt.ylabel( r'$\epsilon$' )
        pdf.savefig()
        plt.close()

        plt.figure()
        plt.scatter( range( T ), numpy.array(time_ev), color='blue', marker='o', label='time' )
        plt.legend(loc='upper left')
        plt.xlabel( 'Particle System' )
        plt.ylabel( 'time' )
        pdf.savefig()
        plt.close()

        plt.figure()
        plt.scatter( range( T ), convergence/max( convergence ), color='green', marker='x', label='convergence' )
        plt.legend()
        plt.xlabel( 'Particle System' )
        plt.ylabel( 'convergence criteria(M/K)' )
        pdf.savefig()
        plt.close()



def plot_2D( T, file_output, Parameters):
    """
    Make 2-dimensional plot for ABC results. 

    input:  T 		-> number of total particle systems (int)  
	    file_output -> name of output figure file (str)
            Parameters 	-> input parameter dictionary - same as used for ABC sampler (dict)

    output: summary pdf file with multiple pages for each generated particle system. 

    *******
    If you want to change the range of parameters for plotting, change values for the key 'param_lim' in the input dictionary parameter! 
    *******
    """

    if file_output[-3:] != 'pdf':
        raise NameError( 'Name file for figure output must be a pdf!' )       


    sampling = numpy.array( [ [ i for i in numpy.arange( Parameters['param_lim'][ j ][0], Parameters['param_lim'][ j ][1], (Parameters['param_lim'][ j ][1]-Parameters['param_lim'][ j ][0])/1000) ] for j in range( len( Parameters['param_to_fit'] ) )] )
   
    sampling2 = numpy.array( [ [ i for i in numpy.arange( Parameters['param_lim'][ j ][0], Parameters['param_lim'][ j ][1], (Parameters['param_lim'][ j ][1]-Parameters['param_lim'][ j ][0])/100) ] for j in range( len( Parameters['param_to_fit'] ) )] )

    



    with PdfPages( file_output ) as pdf:

        epsilon_ev = []
        time_ev = []
        ndraws_ev = []

        for i in range( T ):
           
            op1 = open( Parameters['file_root'] + str( i ) + '.dat', 'r' )
            lin1 = op1.readlines()
            op1.close()

            d = [ elem.split() for elem in lin1 ]
            d1 = numpy.array([ [ float( item ) for item in line ] for line in d[1:] ])
     

            epsilon_ev.append( float( d[1][ d[0].index( 'dist_threshold' ) ] ) )
            time_ev.append( sum( float( line[ d[0].index('time') ] ) for line in d1[1:]) )
            ndraws_ev.append( sum( float( line[ d[0].index('NDraws') ] ) for line in d1[1:]) )

            if i > 0:
                w1 = numpy.loadtxt( Parameters['file_root'] + str( i ) + 'weights.dat' )
            else:
                w1 = numpy.array([ 1.0/len(d1) for k in range( len(d1) ) ])


            kde1 = gaussian_kde(d1[:,0] , weights=w1)
            y1 = kde1( sampling[0] )

            kde2 = gaussian_kde( d1[:,1], weights=w1)
            y2 = kde2( sampling[1] )
 

            d2 = numpy.array( d1[:,:len(Parameters['param_to_fit'])] ) 
            kde3 = gaussian_kde( d2.transpose() )
            xx, yy = numpy.meshgrid( sampling2[0], sampling2[1])
            y3 = kde3(( numpy.ravel(xx), numpy.ravel(yy)))
            zz = numpy.reshape( y3, xx.shape)    

            kwargs = dict(extent=(Parameters['param_lim'][0][0], Parameters['param_lim'][0][1], Parameters['param_lim'][1][0], Parameters['param_lim'][1][1]), cmap='hot', origin='lower')

            #### Plot posteriors
            f = plt.figure()
            gs0 = gridspec.GridSpec( 3, 1, left=0.1, right=0.95, wspace=0.2, hspace=0.5 )
            gs1 = gridspec.GridSpecFromSubplotSpec(2,2, subplot_spec=gs0[:-1], wspace=0.0, hspace=0.0 )
        
            axA = plt.Subplot(f, gs1[:,:])
            f.add_subplot( axA )

            gs2 = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=gs0[-1], wspace=0.25)
            ax1 = plt.Subplot( f, gs2[0] )
            f.add_subplot( ax1 )

            ax2 = plt.Subplot( f, gs2[1] )
            f.add_subplot( ax2 )

            axA.set_title( 'Particle System t = ' + str( i + 1) )
            axA.imshow( zz.T, **kwargs)
            axA.set_xlabel( Parameters['param_to_fit'][0] )
            axA.set_ylabel( Parameters['param_to_fit'][1] )    
            axA.set_aspect('auto')     
            axA.tick_params(axis='both', which='major', labelsize=10)

            ax1.plot( sampling[0], y1, color='blue')
            ax1.set_xlabel( Parameters['param_to_fit'][0] )
            ax1.set_ylabel( 'density', fontsize=8 )
            ax1.tick_params(axis='both', which='major', labelsize=8)
            ax1.set_xlim( Parameters['param_lim'][0][0]-0.1, Parameters['param_lim'][0][1]+0.1)

            ax2.plot( sampling[1], y2, color='blue')
            ax2.set_xlabel( Parameters['param_to_fit'][1] )
            ax2.set_ylabel( 'density', fontsize = 8 )
            ax2.tick_params(axis='both', which='major', labelsize=8)
            ax2.set_xlim( Parameters['param_lim'][1][0]-0.1, Parameters['param_lim'][1][1]+0.1) 

            pdf.savefig()
            plt.close()

    
        convergence = numpy.array( [ float(Parameters['M'])/item  for item in ndraws_ev ] )
    
        #Plot epsilon evolution
        plt.figure()
        plt.scatter( range( T ), numpy.array(epsilon_ev), color='red', marker='*', label='distance threshold' )
        plt.legend()
        plt.xlabel( 'Particle System' )
        plt.ylabel( r'$\epsilon$' )
        pdf.savefig()
        plt.close()

        plt.figure()
        plt.scatter( range( T ), numpy.array(time_ev), color='blue', marker='o', label='time' )
        plt.legend(loc='upper left')
        plt.xlabel( 'Particle System' )
        plt.ylabel( 'time' )
        pdf.savefig()
        plt.close()

        plt.figure()
        plt.scatter( range( T ), convergence/max( convergence ), color='green', marker='x', label='convergence' )
        plt.legend()
        plt.xlabel( 'Particle System' )
        plt.ylabel( 'convergence criteria(M/K)' )
        pdf.savefig()
        plt.close()

