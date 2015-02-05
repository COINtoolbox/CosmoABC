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

    y0 = [ Parameters['prior_func'][0]( Parameters['prior_par'][0], Parameters['param_lim'][0] ) for x in xrange( Parameters['M'] ) ]
    w0 = [ 1.0/len(y0) for i in y0 ]

    kde0 = gaussian_kde(y0 , weights=w0)
    y00 = kde0( sampling[0] )

    epsilon_ev = []
    time_ev = []
    ndraws_ev = []

    marker = ['*','+','x','o','^','d','<']
    color = ['red', 'green','blue','brown', 'gray', 'purple']

    with PdfPages( file_output ) as pdf:

   
        plt.figure()
        plt.title( 'Particle System: 0' , fontsize=15 )
        plt.plot( sampling[0], y00, color='blue')
        plt.xlabel( Parameters['param_to_fit'][0] )
        plt.ylabel( 'density', fontsize=12 )
        plt.tick_params(axis='both', which='major', labelsize=12)
        plt.xlim( Parameters['param_lim'][0][0], Parameters['param_lim'][0][1])
        pdf.savefig()
        plt.close()

        
        print 'Finished plotting particle system T=0'

        for i in range( T ):
           
            op1 = open( Parameters['file_root'] + str( i ) + '.dat', 'r' )
            lin1 = op1.readlines()
            op1.close()

            d = [ elem.split() for elem in lin1 ]
            d1 = numpy.array([ [ float( item ) for item in line ] for line in d[1:] ])
     

            epsilon_ev.append( [ float( d[1][ d[0].index( 'dist_threshold' + str( jj + 1 ) ) ] ) for jj in xrange( len( Parameters['epsilon1'] ) ) ] )
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
            plt.title( 'Particle System: ' + str(i), fontsize=15 )
            plt.plot( sampling[0], y1, color='blue')
            plt.xlabel( Parameters['param_to_fit'][0] )
            plt.ylabel( 'density', fontsize=12 )
            plt.tick_params(axis='both', which='major', labelsize=12)
            plt.xlim( Parameters['param_lim'][0][0], Parameters['param_lim'][0][1])
            pdf.savefig()
            plt.close()

    
        convergence = numpy.array( [ float(Parameters['M'])/item  for item in ndraws_ev ] )
    
        #Plot epsilon evolution
        plt.figure()
        plt.title('Distance threshold evolution')
        for kk in xrange( len(Parameters['epsilon1'] ) ):
            plt.scatter( range( 1,  T + 1 ), numpy.array(epsilon_ev)[:,kk]/max(numpy.array(epsilon_ev)[:,kk]), color=color[ kk ], marker=marker[ kk ], label='distance threshold' + str( kk ) )
        plt.legend()
        plt.xlabel( 'Particle System' )
        plt.ylabel( r'$\epsilon$' )
        pdf.savefig()
        plt.close()

        plt.figure()
        plt.title( 'Computational time evolution' )
        plt.scatter( range( 1, T + 1 ), numpy.array(time_ev), color=color[ len( Parameters['epsilon1'] ) ], marker=marker[len( Parameters['epsilon1'] )], label='time' )
        plt.legend(loc='upper left')
        plt.xlabel( 'Particle System' )
        plt.ylabel( 'time (s)' )
        pdf.savefig()
        plt.close()

        plt.figure()
        plt.title('Convergence criteria evolution')
        plt.scatter( range(1, T + 1), convergence, color=color[len( Parameters['epsilon1'] )+1], marker=marker[len( Parameters['epsilon1'] )+1], label='convergence' )
        plt.legend()
        plt.xlabel( 'Particle System' )
        plt.ylabel( 'convergence criteria(M/K)' )
        pdf.savefig()
        plt.close()

        print 'Finished plotting particle system T=' + str( i + 1 )



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

    y0 = numpy.array([ [ Parameters['prior_func'][0]( Parameters['prior_par'][0], Parameters['param_lim'][0] ), Parameters['prior_func'][1]( Parameters['prior_par'][1], Parameters['param_lim'][1] ) ] for x in xrange( Parameters['M'] ) ])
    
    w0 = [ 1.0/len(y0) for i in y0 ]

    kde01 = gaussian_kde(y0[:,0] , weights=w0)
    y01 = kde01( sampling[0] )

    kde02 = gaussian_kde(y0[:,1], weights=w0)
    y02 = kde02( sampling[1] )

    d20 = numpy.array( y0 ) 
    kde30 = gaussian_kde( d20.transpose() )
    xx, yy = numpy.meshgrid( sampling2[0], sampling2[1])
    y30 = kde30(( numpy.ravel(xx), numpy.ravel(yy)))
    zz0 = numpy.reshape( y30, xx.shape) 

    kwargs = dict(extent=(Parameters['param_lim'][0][0], Parameters['param_lim'][0][1], Parameters['param_lim'][1][0], Parameters['param_lim'][1][1]), cmap='hot', origin='lower')

    epsilon_ev = []
    time_ev = []
    ndraws_ev = []
   
    marker = ['*','+','x','o','^','d','<']
    color = ['red', 'green','blue','brown', 'gray', 'purple']

    with PdfPages( file_output ) as pdf:


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

        axA.set_title( 'Particle System t = 0' )
        axA.imshow( zz0, **kwargs)
        axA.set_xlabel( Parameters['param_to_fit'][0] )
        axA.set_ylabel( Parameters['param_to_fit'][1] )    
        axA.set_aspect('auto')     
        axA.tick_params(axis='both', which='major', labelsize=10)

        ax1.plot( sampling[0], y01, color='blue')
        ax1.set_xlabel( Parameters['param_to_fit'][0] )
        ax1.set_ylabel( 'density', fontsize=8 )
        ax1.tick_params(axis='both', which='major', labelsize=8)
        ax1.set_xlim( Parameters['param_lim'][0][0], Parameters['param_lim'][0][1])

        ax2.plot( sampling[1], y02, color='blue')
        ax2.set_xlabel( Parameters['param_to_fit'][1] )
        ax2.set_ylabel( 'density', fontsize = 8 )
        ax2.tick_params(axis='both', which='major', labelsize=8)
        ax2.set_xlim( Parameters['param_lim'][1][0], Parameters['param_lim'][1][1]) 

        pdf.savefig()
        plt.close()

        print 'Finished plotting particle system T=0' 

        for i in range( T ):
           
            op1 = open( Parameters['file_root'] + str( i ) + '.dat', 'r' )
            lin1 = op1.readlines()
            op1.close()

            d = [ elem.split() for elem in lin1 ]
            d1 = numpy.array([ [ float( item ) for item in line ] for line in d[1:] ])
     

            epsilon_ev.append( [ float( d[1][ d[0].index( 'dist_threshold' + str( jj + 1) ) ] ) for jj in xrange( len( Parameters['epsilon1'] ) )] )
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
            y3 = kde3(( numpy.ravel(xx), numpy.ravel(yy)))
            zz = numpy.reshape( y3, xx.shape)    


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
            axA.imshow( zz, **kwargs)
            axA.set_xlabel( Parameters['param_to_fit'][0] )
            axA.set_ylabel( Parameters['param_to_fit'][1] )    
            axA.set_aspect('auto')     
            axA.tick_params(axis='both', which='major', labelsize=10)

            ax1.plot( sampling[0], y1, color='blue')
            ax1.set_xlabel( Parameters['param_to_fit'][0] )
            ax1.set_ylabel( 'density', fontsize=8 )
            ax1.tick_params(axis='both', which='major', labelsize=8)
            ax1.set_xlim( Parameters['param_lim'][0][0], Parameters['param_lim'][0][1])

            ax2.plot( sampling[1], y2, color='blue')
            ax2.set_xlabel( Parameters['param_to_fit'][1] )
            ax2.set_ylabel( 'density', fontsize = 8 )
            ax2.tick_params(axis='both', which='major', labelsize=8)
            ax2.set_xlim( Parameters['param_lim'][1][0], Parameters['param_lim'][1][1]) 

            pdf.savefig()
            plt.close()
           
            print 'Finished plotting particle system T=' +str( i+1 )

    
        convergence = numpy.array( [ float(Parameters['M'])/item  for item in ndraws_ev ] )
    
        #Plot epsilon evolution
        plt.figure()
        plt.title('Distance threshold evolution')
        for jj in xrange( len( Parameters['epsilon1'] ) ):
            plt.scatter( range( 1, T + 1), numpy.array(epsilon_ev)[:,jj]/max(numpy.array(epsilon_ev)[:,jj]), color=color[ jj ], marker=marker[ jj ], label='distance threshold' + str( jj ) )
        plt.legend()
        plt.xlabel( 'Particle System' )
        plt.ylabel( r'$\epsilon$' )
        pdf.savefig()
        plt.close()

        plt.figure()
        plt.title('Computational time evolution')
        plt.scatter( range( 1, T + 1 ), numpy.array(time_ev), color=color[len( Parameters['epsilon1'] )], marker=marker[len( Parameters['epsilon1'] )], label='time' )
        plt.legend(loc='upper left')
        plt.xlabel( 'Particle System' )
        plt.ylabel( 'time' )
        pdf.savefig()
        plt.close()

        plt.figure()
        plt.title('Convergene criteria evolution')
        plt.scatter( range( 1, T + 1 ), convergence, color=color[ len( Parameters['epsilon1'] ) + 1], marker=marker[len( Parameters['epsilon1'] )+1], label='convergence' )
        plt.legend()
        plt.xlabel( 'Particle System' )
        plt.ylabel( 'convergence criteria(M/K)' )
        pdf.savefig()
        plt.close()


def plot_3D( T, file_output, Parameters):
    """
    Make 3-dimensional plot for ABC results. 

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


    #sampling for 1D plots
    sampling = numpy.array( [ [ i for i in numpy.arange( Parameters['param_lim'][ j ][0], Parameters['param_lim'][ j ][1], (Parameters['param_lim'][ j ][1]-Parameters['param_lim'][ j ][0])/1000) ] for j in range( len( Parameters['param_to_fit'] ) )] )
   

    #sampling for 2D plots
    sampling2 = numpy.array( [ [ i for i in numpy.arange( Parameters['param_lim'][ j ][0], Parameters['param_lim'][ j ][1], (Parameters['param_lim'][ j ][1]-Parameters['param_lim'][ j ][0])/100) ] for j in range( len( Parameters['param_to_fit'] ) )] )

    #define variables for 2D plots 
    xx12, yy12 = numpy.meshgrid( sampling2[0], sampling2[1])
    xx13, yy13 = numpy.meshgrid( sampling2[0], sampling2[2])
    xx23, yy23 = numpy.meshgrid( sampling2[1], sampling2[2])
    
    kwargs12 = dict(extent=(Parameters['param_lim'][0][0], Parameters['param_lim'][0][1], Parameters['param_lim'][1][0], Parameters['param_lim'][1][1]), cmap='hot', origin='lower') 
    kwargs13 = dict(extent=(Parameters['param_lim'][0][0], Parameters['param_lim'][0][1], Parameters['param_lim'][2][0], Parameters['param_lim'][2][1]), cmap='hot', origin='lower')
    kwargs23 = dict(extent=(Parameters['param_lim'][1][0], Parameters['param_lim'][1][1], Parameters['param_lim'][2][0], Parameters['param_lim'][2][1]), cmap='hot', origin='lower')


    y0 = numpy.array([ [ Parameters['prior_func'][ j1 ]( Parameters['prior_par'][ j1 ], Parameters['param_lim'][ j1 ] ) for j1 in xrange(3) ] for x in xrange( Parameters['M'] ) ])
    
    w0 = [ 1.0/len(y0) for i in y0 ]

    kde01 = gaussian_kde(y0[:,0] , weights=w0)
    y01 = kde01( sampling[0] )

    kde02 = gaussian_kde(y0[:,1], weights=w0)
    y02 = kde02( sampling[1] )

    kde03 = gaussian_kde(y0[:,2], weights=w0)
    y03 = kde03( sampling[2] )


    kde012 = gaussian_kde( numpy.array( [y0[:,0], y0[:,1] ]) )  
    y012 = kde012(( numpy.ravel(xx12), numpy.ravel(yy12)))
    zz012 = numpy.reshape( y012, xx12.shape)   
    

    kde013 = gaussian_kde( numpy.array( [y0[:,0], y0[:,2]] ) )   
    y013 = kde013(( numpy.ravel(xx13), numpy.ravel(yy13)))
    zz013 = numpy.reshape( y013, xx13.shape)    
    

    kde023 = gaussian_kde( numpy.array([ y0[:,1], y0[:,2] ]) )
    y023 = kde023(( numpy.ravel(xx23), numpy.ravel(yy23)))
    zz023 = numpy.reshape( y023, xx23.shape)

    epsilon_ev = []
    time_ev = []
    ndraws_ev = []
   
    marker = ['*','+','x','o','^','d','<']
    color = ['red', 'green','blue','brown', 'gray', 'purple']

    with PdfPages( file_output ) as pdf:


        #### Plot posteriors
        f = plt.figure()
           
            
        gs0 = gridspec.GridSpec( 2, 3, left=0.075, right=0.975, wspace=0.35, hspace=0.3 )
            
        ax1 = plt.Subplot(f, gs0[0])
        f.add_subplot( ax1 )

        ax2 = plt.Subplot( f, gs0[1] )
        f.add_subplot( ax2 )

        ax3 = plt.Subplot( f, gs0[2] )
        f.add_subplot( ax3 )

        ax4 = plt.Subplot( f, gs0[3] )
        f.add_subplot( ax4 )

        ax5 = plt.Subplot( f, gs0[4] )
        f.add_subplot( ax5 )

        ax6 = plt.Subplot( f, gs0[5] )
        f.add_subplot( ax6 )

        ax1.imshow( zz012, **kwargs12)
        ax1.set_xlabel( Parameters['param_to_fit'][0] )
        ax1.set_ylabel( Parameters['param_to_fit'][1] )    
        ax1.set_aspect('auto')     
        ax1.tick_params(axis='both', which='major', labelsize=8)
            
        ax2.set_title( 'Particle System t = 0' )
        ax2.imshow( zz013, **kwargs13)
        ax2.set_xlabel( Parameters['param_to_fit'][0] )
        ax2.set_ylabel( Parameters['param_to_fit'][2] )    
        ax2.set_aspect('auto')     
        ax2.tick_params(axis='both', which='major', labelsize=8)
   
        ax3.imshow( zz023, **kwargs23)
        ax3.set_xlabel( Parameters['param_to_fit'][1] )
        ax3.set_ylabel( Parameters['param_to_fit'][2] )    
        ax3.set_aspect('auto')     
        ax3.tick_params(axis='both', which='major', labelsize=8) 

        ax4.plot( sampling[0], y01, color='blue')
        ax4.set_xlabel( Parameters['param_to_fit'][0] )
        ax4.set_ylabel( 'density', fontsize=8 )
        ax4.tick_params(axis='both', which='major', labelsize=8)
        ax4.set_xlim( Parameters['param_lim'][0][0], Parameters['param_lim'][0][1])


        ax5.plot( sampling[1], y02, color='red')
        ax5.set_xlabel( Parameters['param_to_fit'][1] )
        ax5.set_ylabel( 'density', fontsize = 8 )
        ax5.tick_params(axis='both', which='major', labelsize=8)
        ax5.set_xlim( Parameters['param_lim'][1][0], Parameters['param_lim'][1][1]) 

        ax6.plot( sampling[2], y03, color='green')
        ax6.set_xlabel( Parameters['param_to_fit'][2] )
        ax6.set_ylabel( 'density', fontsize = 8 )
        ax6.tick_params(axis='both', which='major', labelsize=8)
        ax6.set_xlim( Parameters['param_lim'][2][0], Parameters['param_lim'][2][1]) 

        pdf.savefig()
        plt.close()

        print 'Finished plotting particle system T=0'

        for i in range( T ):
           
            #read individual particle systems
            op1 = open( Parameters['file_root'] + str( i ) + '.dat', 'r' )
            lin1 = op1.readlines()
            op1.close()


            d = [ elem.split() for elem in lin1 ]
            d1 = numpy.array([ [ float( item ) for item in line ] for line in d[1:] ])
     

            epsilon_ev.append( [ float( d[1][ d[0].index( 'dist_threshold' + str( jj + 1 ) ) ] )  for jj in xrange( len( Parameters['epsilon1'] ) ) ] )
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

            kde3 = gaussian_kde( d1[:,2], weights=w1)
            y3 = kde3( sampling[2] ) 
 

            
            kde12 = gaussian_kde( numpy.array( [d1[:,0], d1[:,1] ]) )
            xx12, yy12 = numpy.meshgrid( sampling2[0], sampling2[1])
            y12 = kde12(( numpy.ravel(xx12), numpy.ravel(yy12)))
            zz12 = numpy.reshape( y12, xx12.shape)   
            kwargs12 = dict(extent=(Parameters['param_lim'][0][0], Parameters['param_lim'][0][1], Parameters['param_lim'][1][0], Parameters['param_lim'][1][1]), cmap='hot', origin='lower') 

            kde13 = gaussian_kde( numpy.array( [d1[:,0], d1[:,2]] ) )
            xx13, yy13 = numpy.meshgrid( sampling2[0], sampling2[2])
            y13 = kde13(( numpy.ravel(xx13), numpy.ravel(yy13)))
            zz13 = numpy.reshape( y13, xx13.shape)    
            kwargs13 = dict(extent=(Parameters['param_lim'][0][0], Parameters['param_lim'][0][1], Parameters['param_lim'][2][0], Parameters['param_lim'][2][1]), cmap='hot', origin='lower')

            kde23 = gaussian_kde( numpy.array([ d1[:,1], d1[:,2] ]) )
            xx23, yy23 = numpy.meshgrid( sampling2[1], sampling2[2])
            y23 = kde23(( numpy.ravel(xx23), numpy.ravel(yy23)))
            zz23 = numpy.reshape( y23, xx23.shape)
            kwargs23 = dict(extent=(Parameters['param_lim'][1][0], Parameters['param_lim'][1][1], Parameters['param_lim'][2][0], Parameters['param_lim'][2][1]), cmap='hot', origin='lower')

            

            #### Plot posteriors
            f = plt.figure()
           
            
            gs0 = gridspec.GridSpec( 2, 3, left=0.075, right=0.975, wspace=0.35, hspace=0.3 )
            
        
            ax1 = plt.Subplot(f, gs0[0])
            f.add_subplot( ax1 )

            ax2 = plt.Subplot( f, gs0[1] )
            f.add_subplot( ax2 )

            ax3 = plt.Subplot( f, gs0[2] )
            f.add_subplot( ax3 )

            ax4 = plt.Subplot( f, gs0[3] )
            f.add_subplot( ax4 )

            ax5 = plt.Subplot( f, gs0[4] )
            f.add_subplot( ax5 )

            ax6 = plt.Subplot( f, gs0[5] )
            f.add_subplot( ax6 )

            ax1.imshow( zz12, **kwargs12)
            ax1.set_xlabel( Parameters['param_to_fit'][0] )
            ax1.set_ylabel( Parameters['param_to_fit'][1] )    
            ax1.set_aspect('auto')     
            ax1.tick_params(axis='both', which='major', labelsize=8)

            
            ax2.set_title( 'Particle System t = ' + str( i + 1) )
            ax2.imshow( zz13, **kwargs13)
            ax2.set_xlabel( Parameters['param_to_fit'][0] )
            ax2.set_ylabel( Parameters['param_to_fit'][2] )    
            ax2.set_aspect('auto')     
            ax2.tick_params(axis='both', which='major', labelsize=8)

    
            ax3.imshow( zz23, **kwargs23)
            ax3.set_xlabel( Parameters['param_to_fit'][1] )
            ax3.set_ylabel( Parameters['param_to_fit'][2] )    
            ax3.set_aspect('auto')     
            ax3.tick_params(axis='both', which='major', labelsize=8) 

            ax4.plot( sampling[0], y1, color='blue')
            ax4.set_xlabel( Parameters['param_to_fit'][0] )
            ax4.set_ylabel( 'density', fontsize=8 )
            ax4.tick_params(axis='both', which='major', labelsize=8)
            ax4.set_xlim( Parameters['param_lim'][0][0], Parameters['param_lim'][0][1])


            ax5.plot( sampling[1], y2, color='red')
            ax5.set_xlabel( Parameters['param_to_fit'][1] )
            ax5.set_ylabel( 'density', fontsize = 8 )
            ax5.tick_params(axis='both', which='major', labelsize=8)
            ax5.set_xlim( Parameters['param_lim'][1][0], Parameters['param_lim'][1][1]) 

            ax6.plot( sampling[2], y3, color='green')
            ax6.set_xlabel( Parameters['param_to_fit'][2] )
            ax6.set_ylabel( 'density', fontsize = 8 )
            ax6.tick_params(axis='both', which='major', labelsize=8)
            ax6.set_xlim( Parameters['param_lim'][2][0], Parameters['param_lim'][2][1]) 

            pdf.savefig()
            plt.close()

            print 'Finished plotting particle system T=' + str( i+1 )
    
        convergence = numpy.array( [ float(Parameters['M'])/item  for item in ndraws_ev ] )
    
        #Plot epsilon evolution
        plt.figure()
        for jj in xrange(  len( Parameters['epsilon1'] ) ):
            plt.scatter( range( 1,T + 1), numpy.array(epsilon_ev)[:, jj ]/max(numpy.array(epsilon_ev)[:,jj]), color=color[ jj ], marker=marker[ jj ], label='distance threshold' + str( jj ) )
        plt.legend()
        plt.xlabel( 'Particle System' )
        plt.ylabel( r'$\epsilon$' )
        pdf.savefig()
        plt.close()

        plt.figure()
        plt.scatter( range( 1, T+1 ), numpy.array(time_ev), color=color[len( Parameters['epsilon1'] )], marker=marker[ len( Parameters['epsilon1'] ) ], label='time' )
        plt.legend(loc='upper left')
        plt.xlabel( 'Particle System' )
        plt.ylabel( 'time' )
        pdf.savefig()
        plt.close()

        plt.figure()
        plt.scatter( range( 1, T + 1 ), convergence, color=color[len( Parameters['epsilon1'] ) + 1], marker=marker[len( Parameters['epsilon1'] ) + 1 ], label='convergence' )
        plt.legend()
        plt.xlabel( 'Particle System' )
        plt.ylabel( 'convergence criteria(M/K)' )
        pdf.savefig()
        plt.close()


def plot_4D( T, file_output, Parameters):
    """
    Make 3-dimensional plot for ABC results. 

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


    #sampling for 1D plots
    sampling = numpy.array( [ [ i for i in numpy.arange( Parameters['param_lim'][ j ][0], Parameters['param_lim'][ j ][1], (Parameters['param_lim'][ j ][1]-Parameters['param_lim'][ j ][0])/1000) ] for j in range( len( Parameters['param_to_fit'] ) )] )
   

    #sampling for 2D plots
    sampling2 = numpy.array( [ [ i for i in numpy.arange( Parameters['param_lim'][ j ][0], Parameters['param_lim'][ j ][1], (Parameters['param_lim'][ j ][1]-Parameters['param_lim'][ j ][0])/100) ] for j in range( len( Parameters['param_to_fit'] ) )] )

    #define variables for 2D plots 
    xx12, yy12 = numpy.meshgrid( sampling2[0], sampling2[1])
    xx13, yy13 = numpy.meshgrid( sampling2[0], sampling2[2])
    xx14, yy14 = numpy.meshgrid( sampling2[0], sampling2[3])
    xx23, yy23 = numpy.meshgrid( sampling2[1], sampling2[2])
    xx24, yy24 = numpy.meshgrid( sampling2[1], sampling2[3])
    xx34, yy34 = numpy.meshgrid( sampling2[2], sampling2[3])
    
    kwargs12 = dict(extent=(Parameters['param_lim'][0][0], Parameters['param_lim'][0][1], Parameters['param_lim'][1][0], Parameters['param_lim'][1][1]), cmap='hot', origin='lower') 
    kwargs13 = dict(extent=(Parameters['param_lim'][0][0], Parameters['param_lim'][0][1], Parameters['param_lim'][2][0], Parameters['param_lim'][2][1]), cmap='hot', origin='lower')
    kwargs14 = dict(extent=(Parameters['param_lim'][0][0], Parameters['param_lim'][0][1], Parameters['param_lim'][3][0], Parameters['param_lim'][3][1]), cmap='hot', origin='lower')
    kwargs23 = dict(extent=(Parameters['param_lim'][1][0], Parameters['param_lim'][1][1], Parameters['param_lim'][2][0], Parameters['param_lim'][2][1]), cmap='hot', origin='lower')
    kwargs24 = dict(extent=(Parameters['param_lim'][1][0], Parameters['param_lim'][1][1], Parameters['param_lim'][3][0], Parameters['param_lim'][3][1]), cmap='hot', origin='lower')
    kwargs34 = dict(extent=(Parameters['param_lim'][2][0], Parameters['param_lim'][2][1], Parameters['param_lim'][3][0], Parameters['param_lim'][3][1]), cmap='hot', origin='lower')


    y0 = numpy.array([ [ Parameters['prior_func'][ j1 ]( Parameters['prior_par'][ j1 ], Parameters['param_lim'][ j1 ] ) for j1 in xrange(4) ] for x in xrange( Parameters['M'] ) ])
    
    w0 = [ 1.0/len(y0) for i in y0 ]

    kde01 = gaussian_kde(y0[:,0] , weights=w0)
    y01 = kde01( sampling[0] )

    kde02 = gaussian_kde(y0[:,1], weights=w0)
    y02 = kde02( sampling[1] )

    kde03 = gaussian_kde(y0[:,2], weights=w0)
    y03 = kde03( sampling[2] )

    kde04 = gaussian_kde(y0[:,3], weights=w0)
    y04 = kde04( sampling[3] )


    kde012 = gaussian_kde( numpy.array( [y0[:,0], y0[:,1] ]) )  
    y012 = kde012(( numpy.ravel(xx12), numpy.ravel(yy12)))
    zz012 = numpy.reshape( y012, xx12.shape)   
    

    kde013 = gaussian_kde( numpy.array( [y0[:,0], y0[:,2]] ) )   
    y013 = kde013(( numpy.ravel(xx13), numpy.ravel(yy13)))
    zz013 = numpy.reshape( y013, xx13.shape)  

    kde014 = gaussian_kde( numpy.array( [y0[:,0], y0[:,3]] ) )   
    y014 = kde014(( numpy.ravel(xx14), numpy.ravel(yy14)))
    zz014 = numpy.reshape( y014, xx14.shape)   
    

    kde023 = gaussian_kde( numpy.array([ y0[:,1], y0[:,2] ]) )
    y023 = kde023(( numpy.ravel(xx23), numpy.ravel(yy23)))
    zz023 = numpy.reshape( y023, xx23.shape)

    kde024 = gaussian_kde( numpy.array([ y0[:,1], y0[:,3] ]) )
    y024 = kde024(( numpy.ravel(xx24), numpy.ravel(yy24)))
    zz024 = numpy.reshape( y024, xx24.shape)

    kde034 = gaussian_kde( numpy.array([ y0[:,2], y0[:,3] ]) )
    y034 = kde034(( numpy.ravel(xx34), numpy.ravel(yy34)))
    zz034 = numpy.reshape( y034, xx34.shape)

    epsilon_ev = []
    time_ev = []
    ndraws_ev = []
   
    marker = ['*','+','x','o','^','d','<']
    color = ['red', 'green','blue','brown', 'gray', 'purple']

    with PdfPages( file_output ) as pdf:


        #### Plot posteriors
        f = plt.figure()
           
            
        gs0 = gridspec.GridSpec( 3, 4, left=0.075, right=0.975, wspace=0.35, hspace=0.3 )
            
        ax1 = plt.Subplot(f, gs0[0])
        f.add_subplot( ax1 )

        ax2 = plt.Subplot( f, gs0[1] )
        f.add_subplot( ax2 )

        ax3 = plt.Subplot( f, gs0[2] )
        f.add_subplot( ax3 )

        ax4 = plt.Subplot( f, gs0[3] )
        f.add_subplot( ax4 )

        ax5 = plt.Subplot( f, gs0[4] )
        f.add_subplot( ax5 )

        ax6 = plt.Subplot( f, gs0[5] )
        f.add_subplot( ax6 )

        ax7 = plt.Subplot( f, gs0[6] )
        f.add_subplot( ax7 )

        ax8 = plt.Subplot( f, gs0[7] )
        f.add_subplot( ax8 )

        ax9 = plt.Subplot( f, gs0[8] )
        f.add_subplot( ax9 )

        ax10 = plt.Subplot( f, gs0[9] )
        f.add_subplot( ax10 )

        ax1.imshow( zz012, **kwargs12)
        ax1.set_xlabel( Parameters['param_to_fit'][0] )
        ax1.set_ylabel( Parameters['param_to_fit'][1] )    
        ax1.set_aspect('auto')     
        ax1.tick_params(axis='both', which='major', labelsize=8)
            
        ax2.set_title( 'Particle System t = 0' )
        ax2.imshow( zz013, **kwargs13)
        ax2.set_xlabel( Parameters['param_to_fit'][0] )
        ax2.set_ylabel( Parameters['param_to_fit'][2] )    
        ax2.set_aspect('auto')     
        ax2.tick_params(axis='both', which='major', labelsize=8)
   
        ax3.imshow( zz014, **kwargs14)
        ax3.set_xlabel( Parameters['param_to_fit'][0] )
        ax3.set_ylabel( Parameters['param_to_fit'][3] )    
        ax3.set_aspect('auto')     
        ax3.tick_params(axis='both', which='major', labelsize=8) 

        ax4.imshow( zz023, **kwargs23)
        ax4.set_xlabel( Parameters['param_to_fit'][1] )
        ax4.set_ylabel( Parameters['param_to_fit'][2] )    
        ax4.set_aspect('auto')     
        ax4.tick_params(axis='both', which='major', labelsize=8) 

        ax5.imshow( zz024, **kwargs24)
        ax5.set_xlabel( Parameters['param_to_fit'][1] )
        ax5.set_ylabel( Parameters['param_to_fit'][3] )    
        ax5.set_aspect('auto')     
        ax5.tick_params(axis='both', which='major', labelsize=8)  

        ax6.imshow( zz034, **kwargs34)
        ax6.set_xlabel( Parameters['param_to_fit'][2] )
        ax6.set_ylabel( Parameters['param_to_fit'][3] )    
        ax6.set_aspect('auto')     
        ax6.tick_params(axis='both', which='major', labelsize=8) 


        ax7.plot( sampling[0], y01, color='blue')
        ax7.set_xlabel( Parameters['param_to_fit'][0] )
        ax7.set_ylabel( 'density', fontsize=8 )
        ax7.tick_params(axis='both', which='major', labelsize=8)
        ax7.set_xlim( Parameters['param_lim'][0][0], Parameters['param_lim'][0][1])

        ax8.plot( sampling[1], y02, color='red')
        ax8.set_xlabel( Parameters['param_to_fit'][1] )
        ax8.set_ylabel( 'density', fontsize = 8 )
        ax8.tick_params(axis='both', which='major', labelsize=8)
        ax8.set_xlim( Parameters['param_lim'][1][0], Parameters['param_lim'][1][1]) 

        ax9.plot( sampling[2], y03, color='green')
        ax9.set_xlabel( Parameters['param_to_fit'][2] )
        ax9.set_ylabel( 'density', fontsize = 8 )
        ax9.tick_params(axis='both', which='major', labelsize=8)
        ax9.set_xlim( Parameters['param_lim'][2][0], Parameters['param_lim'][2][1]) 

        ax10.plot( sampling[3], y04, color='green')
        ax10.set_xlabel( Parameters['param_to_fit'][2] )
        ax10.set_ylabel( 'density', fontsize = 8 )
        ax10.tick_params(axis='both', which='major', labelsize=8)
        ax10.set_xlim( Parameters['param_lim'][3][0], Parameters['param_lim'][3][1])

        pdf.savefig()
        plt.close()

        print 'Finished plotting particle system T=0'

        for i in range( T ):
           
            #read individual particle systems
            op1 = open( Parameters['file_root'] + str( i ) + '.dat', 'r' )
            lin1 = op1.readlines()
            op1.close()


            d = [ elem.split() for elem in lin1 ]
            d1 = numpy.array([ [ float( item ) for item in line ] for line in d[1:] ])
     

            epsilon_ev.append( [ float( d[1][ d[0].index( 'dist_threshold' + str( jj + 1 ) ) ] )  for jj in xrange( len( Parameters['epsilon1'] ) ) ] )
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

            kde3 = gaussian_kde( d1[:,2], weights=w1)
            y3 = kde3( sampling[2] ) 
 
            kde4 = gaussian_kde( d1[:,3], weights=w1)
            y4 = kde3( sampling[3] ) 

            
            kde12 = gaussian_kde( numpy.array( [d1[:,0], d1[:,1] ]) )
            xx12, yy12 = numpy.meshgrid( sampling2[0], sampling2[1])
            y12 = kde12(( numpy.ravel(xx12), numpy.ravel(yy12)))
            zz12 = numpy.reshape( y12, xx12.shape)   
            kwargs12 = dict(extent=(Parameters['param_lim'][0][0], Parameters['param_lim'][0][1], Parameters['param_lim'][1][0], Parameters['param_lim'][1][1]), cmap='hot', origin='lower') 

            kde13 = gaussian_kde( numpy.array( [d1[:,0], d1[:,2]] ) )
            xx13, yy13 = numpy.meshgrid( sampling2[0], sampling2[2])
            y13 = kde13(( numpy.ravel(xx13), numpy.ravel(yy13)))
            zz13 = numpy.reshape( y13, xx13.shape)    
            kwargs13 = dict(extent=(Parameters['param_lim'][0][0], Parameters['param_lim'][0][1], Parameters['param_lim'][2][0], Parameters['param_lim'][2][1]), cmap='hot', origin='lower')

            kde14 = gaussian_kde( numpy.array( [d1[:,0], d1[:,3]] ) )
            xx14, yy14 = numpy.meshgrid( sampling2[0], sampling2[3])
            y14 = kde14(( numpy.ravel(xx14), numpy.ravel(yy14)))
            zz14 = numpy.reshape( y14, xx14.shape)    
            kwargs14 = dict(extent=(Parameters['param_lim'][0][0], Parameters['param_lim'][0][1], Parameters['param_lim'][3][0], Parameters['param_lim'][3][1]), cmap='hot', origin='lower') 

            kde23 = gaussian_kde( numpy.array([ d1[:,1], d1[:,2] ]) )
            xx23, yy23 = numpy.meshgrid( sampling2[1], sampling2[2])
            y23 = kde23(( numpy.ravel(xx23), numpy.ravel(yy23)))
            zz23 = numpy.reshape( y23, xx23.shape)
            kwargs23 = dict(extent=(Parameters['param_lim'][1][0], Parameters['param_lim'][1][1], Parameters['param_lim'][2][0], Parameters['param_lim'][2][1]), cmap='hot', origin='lower')

            kde24 = gaussian_kde( numpy.array([ d1[:,1], d1[:,3] ]) )
            xx24, yy24 = numpy.meshgrid( sampling2[1], sampling2[3])
            y24 = kde24(( numpy.ravel(xx24), numpy.ravel(yy24)))
            zz24 = numpy.reshape( y24, xx24.shape)
            kwargs24 = dict(extent=(Parameters['param_lim'][1][0], Parameters['param_lim'][1][1], Parameters['param_lim'][3][0], Parameters['param_lim'][3][1]), cmap='hot', origin='lower')

            kde34 = gaussian_kde( numpy.array([ d1[:,2], d1[:,3] ]) )
            xx34, yy34 = numpy.meshgrid( sampling2[2], sampling2[3])
            y34 = kde34(( numpy.ravel(xx34), numpy.ravel(yy34)))
            zz34 = numpy.reshape( y34, xx34.shape)
            kwargs34 = dict(extent=(Parameters['param_lim'][2][0], Parameters['param_lim'][2][1], Parameters['param_lim'][3][0], Parameters['param_lim'][3][1]), cmap='hot', origin='lower')

            

            #### Plot posteriors
            f = plt.figure()
           
            
            gs0 = gridspec.GridSpec( 3, 4, left=0.075, right=0.975, wspace=0.35, hspace=0.3 )
            
        
            ax1 = plt.Subplot(f, gs0[0])
            f.add_subplot( ax1 )

            ax2 = plt.Subplot( f, gs0[1] )
            f.add_subplot( ax2 )

            ax3 = plt.Subplot( f, gs0[2] )
            f.add_subplot( ax3 )

            ax4 = plt.Subplot( f, gs0[3] )
            f.add_subplot( ax4 )

            ax5 = plt.Subplot( f, gs0[4] )
            f.add_subplot( ax5 )

            ax6 = plt.Subplot( f, gs0[5] )
            f.add_subplot( ax6 )

            ax7 = plt.Subplot( f, gs0[6] )
            f.add_subplot( ax7 )

            ax8 = plt.Subplot( f, gs0[7] )
            f.add_subplot( ax8 )

            ax9 = plt.Subplot( f, gs0[8] )
            f.add_subplot( ax9 )

            ax10 = plt.Subplot( f, gs0[9] )
            f.add_subplot( ax10 )
   
            ax1.imshow( zz12, **kwargs12)
            ax1.set_xlabel( Parameters['param_to_fit'][0] )
            ax1.set_ylabel( Parameters['param_to_fit'][1] )    
            ax1.set_aspect('auto')     
            ax1.tick_params(axis='both', which='major', labelsize=8)

            
            ax2.set_title( 'Particle System t = ' + str( i + 1) )
            ax2.imshow( zz13, **kwargs13)
            ax2.set_xlabel( Parameters['param_to_fit'][0] )
            ax2.set_ylabel( Parameters['param_to_fit'][2] )    
            ax2.set_aspect('auto')     
            ax2.tick_params(axis='both', which='major', labelsize=8)

            ax3.set_title( 'Particle System t = ' + str( i + 1) )
            ax3.imshow( zz14, **kwargs14)
            ax3.set_xlabel( Parameters['param_to_fit'][0] )
            ax3.set_ylabel( Parameters['param_to_fit'][3] )    
            ax3.set_aspect('auto')     
            ax3.tick_params(axis='both', which='major', labelsize=8)

            ax4.imshow( zz23, **kwargs23)
            ax4.set_xlabel( Parameters['param_to_fit'][1] )
            ax4.set_ylabel( Parameters['param_to_fit'][2] )    
            ax4.set_aspect('auto')     
            ax4.tick_params(axis='both', which='major', labelsize=8) 

            ax5.imshow( zz24, **kwargs24)
            ax5.set_xlabel( Parameters['param_to_fit'][1] )
            ax5.set_ylabel( Parameters['param_to_fit'][3] )    
            ax5.set_aspect('auto')     
            ax5.tick_params(axis='both', which='major', labelsize=8) 

            ax6.imshow( zz34, **kwargs34)
            ax6.set_xlabel( Parameters['param_to_fit'][2] )
            ax6.set_ylabel( Parameters['param_to_fit'][3] )    
            ax6.set_aspect('auto')     
            ax6.tick_params(axis='both', which='major', labelsize=8)

            ax7.plot( sampling[0], y1, color='blue')
            ax7.set_xlabel( Parameters['param_to_fit'][0] )
            ax7.set_ylabel( 'density', fontsize=8 )
            ax7.tick_params(axis='both', which='major', labelsize=8)
            ax7.set_xlim( Parameters['param_lim'][0][0], Parameters['param_lim'][0][1])


            ax8.plot( sampling[1], y2, color='red')
            ax8.set_xlabel( Parameters['param_to_fit'][1] )
            ax8.set_ylabel( 'density', fontsize = 8 )
            ax8.tick_params(axis='both', which='major', labelsize=8)
            ax8.set_xlim( Parameters['param_lim'][1][0], Parameters['param_lim'][1][1]) 

            ax9.plot( sampling[2], y3, color='green')
            ax9.set_xlabel( Parameters['param_to_fit'][2] )
            ax9.set_ylabel( 'density', fontsize = 8 )
            ax9.tick_params(axis='both', which='major', labelsize=8)
            ax9.set_xlim( Parameters['param_lim'][2][0], Parameters['param_lim'][2][1]) 
   
            ax10.plot( sampling[3], y3, color='green')
            ax10.set_xlabel( Parameters['param_to_fit'][3] )
            ax10.set_ylabel( 'density', fontsize = 8 )
            ax10.tick_params(axis='both', which='major', labelsize=8)
            ax10.set_xlim( Parameters['param_lim'][3][0], Parameters['param_lim'][3][1]) 

            pdf.savefig()
            plt.close()

            print 'Finished plotting particle system T=' + str( i+1 )
    
        convergence = numpy.array( [ float(Parameters['M'])/item  for item in ndraws_ev ] )
    
        #Plot epsilon evolution
        plt.figure()
        for jj in xrange(  len( Parameters['epsilon1'] ) ):
            plt.scatter( range( 1,T + 1), numpy.array(epsilon_ev)[:, jj ]/max(numpy.array(epsilon_ev)[:,jj]), color=color[ jj ], marker=marker[ jj ], label='distance threshold' + str( jj ) )
        plt.legend()
        plt.xlabel( 'Particle System' )
        plt.ylabel( r'$\epsilon$' )
        pdf.savefig()
        plt.close()

        plt.figure()
        plt.scatter( range( 1, T+1 ), numpy.array(time_ev), color=color[len( Parameters['epsilon1'] )], marker=marker[ len( Parameters['epsilon1'] ) ], label='time' )
        plt.legend(loc='upper left')
        plt.xlabel( 'Particle System' )
        plt.ylabel( 'time' )
        pdf.savefig()
        plt.close()

        plt.figure()
        plt.scatter( range( 1, T + 1 ), convergence, color=color[len( Parameters['epsilon1'] ) + 1], marker=marker[len( Parameters['epsilon1'] ) + 1 ], label='convergence' )
        plt.legend()
        plt.xlabel( 'Particle System' )
        plt.ylabel( 'convergence criteria(M/K)' )
        pdf.savefig()
        plt.close()
