"""
Functions for plotting.
"""

import datetime
import numpy as np

from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

from weighted_gaussian_kde import gaussian_kde


def plot_1p(T, file_output, Parameters):
    """
    Make 1 free parameter plot for ABC results. 

    input:  T 		-> number of total particle systems (int)  
	    file_output -> name of output figure file (str)
            Parameters 	-> input parameter dictionary - same as used for ABC sampler (dict)

    output: summary pdf file with multiple pages for each generated particle system. 

    *******
    If you want to change the range of parameters for plotting, change values for the key 'param_lim' in the input dictionary parameter! 
    *******
    """

    par = Parameters['param_to_fit'][0]

    if file_output[-3:] != 'pdf':
        raise NameError('Name file for figure output must be a pdf!')       

    sampling = np.array([[i1 for i1 in np.arange(Parameters['prior'][element]['min'], 
                    Parameters['prior'][element]['max'], 
                   (Parameters['prior'][element]['max']-Parameters['prior'][element]['min'])/1000) ] 
                              for element in Parameters['param_to_fit']])

    par1 = Parameters['param_to_fit'][0]
    y0 = [Parameters['prior'][par1]['func'](Parameters['prior'][par1]) for x in xrange(Parameters['M'])]
    w0 = [1.0/len(y0) for i1 in y0]

    kde0 = gaussian_kde(y0 , weights=w0)
    y00 = kde0(sampling[0])

    epsilon_ev = []
    time_ev = []
    ndraws_ev = []

    marker = ['*', '+', 'x', 'o', '^', 'd', '<']
    color = ['red', 'green', 'blue','brown', 'gray', 'purple']

    with PdfPages(file_output) as pdf:

        plt.figure()
        plt.title('Particle System: 0', fontsize=15)
        plt.plot(sampling[0], y00, color='blue')
        plt.xlabel(Parameters['param_to_fit'][0])
        plt.ylabel('density', fontsize=12)
        plt.tick_params(axis='both', which='major', labelsize=12)
        plt.xlim(Parameters['prior'][par1]['min'], Parameters['prior'][par1]['max'])
        pdf.savefig()
        plt.close()

        print 'Finished plotting particle system T=0'

        for i in range(T+1):
           
            op1 = open(Parameters['file_root'] + str(i) + '.dat', 'r')
            lin1 = op1.readlines()
            op1.close()

            d = [elem.split() for elem in lin1]
            d1 = np.array([[float(item) for item in line] for line in d[1:]])
     
            epsilon_ev.append([float(d[1][d[0].index('dist_threshold' + str(jj + 1))]) 
                                               for jj in xrange(Parameters['dist_dim'])])
            time_ev.append(sum(float(line[d[0].index('time')]) for line in d1[1:]))
            ndraws_ev.append(sum(float(line[d[0].index('NDraws')]) for line in d1[1:]))

            if i > 0:
                w1 = np.loadtxt(Parameters['file_root'] + str(i) + 'weights.dat')
            else:
                w1 = np.array([1.0/len(d1) for k in range(len(d1))])

            kde1 = gaussian_kde(d1[:,0] , weights=w1)
            y1 = kde1(sampling[0])

            #### Plot posteriors
            plt.figure()
            plt.title('Particle System: ' + str(i+1), fontsize=15)
            plt.plot(sampling[0], y1, color='blue')
            plt.xlabel(Parameters['param_to_fit'][0])
            plt.ylabel('density', fontsize=12)
            plt.tick_params(axis='both', which='major', labelsize=12)
            plt.xlim(Parameters['prior'][par1]['min'], Parameters['prior'][par1]['max'])
            pdf.savefig()
            plt.close()

        convergence = np.array([float(Parameters['M'])/item  for item in ndraws_ev])
    
        #Plot epsilon evolution
        plt.figure()
        plt.title('Distance threshold evolution')
        for kk in xrange(Parameters['dist_dim']):
            plt.scatter(range(1,  T + 2), np.array(epsilon_ev)[:,kk]/max(np.array(epsilon_ev)[:,kk]), 
                        color=color[kk], marker=marker[kk], label='distance threshold' + str(kk + 1))
        plt.legend()
        plt.xlabel('Particle System')
        plt.ylabel(r'$\epsilon_{norm}$')
        pdf.savefig()
        plt.close()

        plt.figure()
        plt.title('Computational time evolution')
        plt.scatter(range(1, T + 2), np.array(time_ev), 
                   color=color[Parameters['dist_dim']], marker=marker[Parameters['dist_dim']], 
                   label='time' )
        plt.legend(loc='upper left')
        plt.xlabel('Particle System')
        plt.ylabel('time (s)')
        pdf.savefig()
        plt.close()

        plt.figure()
        plt.title('Convergence criteria evolution')
        plt.scatter(range(1, T + 2), convergence, color=color[Parameters['dist_dim'] + 1], 
                    marker=marker[Parameters['dist_dim'] + 1], label='convergence')
        plt.legend()
        plt.xlabel('Particle System')
        plt.ylabel('convergence criteria(M/K)')
        pdf.savefig()
        plt.close()

        print 'Finished plotting particle system T=' + str(i + 1)

def plot_2p(T, file_output, Parameters):
    """
    Make 2 free parameter plot for ABC results. 

    input:  T 		-> number of total particle systems (int)  
	    file_output -> name of output figure file (str)
            Parameters 	-> input parameter dictionary - same as used for ABC sampler (dict)

    output: summary pdf file with multiple pages for each generated particle system. 

    *******
    If you want to change the range of parameters for plotting, change values for the key 'param_lim' in the input dictionary parameter! 
    *******
    """

    p1 = Parameters['param_to_fit'][0]
    p2 = Parameters['param_to_fit'][1]

    if file_output[-3:] != 'pdf':
        raise NameError('Name file for figure output must be a pdf!')       

    sampling = np.array([[i for i in np.arange(Parameters['prior'][par]['min'], 
                          Parameters['prior'][par]['max'], (Parameters['prior'][par]['max'] -
                          Parameters['prior'][par]['min'])/1000)] 
                                    for par in Parameters['param_to_fit']])
   
    sampling2 = np.array([[i for i in np.arange(Parameters['prior'][par]['min'], 
                           Parameters['prior'][par]['max'], (Parameters['prior'][par]['max'] - 
                           Parameters['prior'][par]['min'])/100)] 
                                     for par in Parameters['param_to_fit']])

    y0 = np.array([[Parameters['prior'][par]['func'](Parameters['prior'][par]) for par in Parameters['param_to_fit']] 
                                 for x in xrange(Parameters['M'])])
    
    w0 = [1.0/len(y0) for i in y0]

    kde01 = gaussian_kde(y0[:,0] , weights=w0)
    y01 = kde01(sampling[0])

    kde02 = gaussian_kde(y0[:,1], weights=w0)
    y02 = kde02(sampling[1])

    d20 = np.array(y0) 
    kde30 = gaussian_kde(d20.transpose(), weights=w0)
    xx, yy = np.meshgrid(sampling2[0], sampling2[1])
    y30 = kde30((np.ravel(xx), np.ravel(yy)))
    zz0 = np.reshape(y30, xx.shape) 

    
    kwargs = dict(extent=([Parameters['prior'][par]['min'] for par in Parameters['param_to_fit']] +
                          [Parameters['prior'][par]['max'] for par in Parameters['param_to_fit']]), 
                          cmap='hot', origin='lower')

    epsilon_ev = []
    time_ev = []
    ndraws_ev = []
   
    marker = ['*', '+', 'x', 'o', '^', 'd', '<']
    color = ['red', 'green','blue','brown', 'gray', 'purple']

    with PdfPages(file_output) as pdf:

        #### Plot posteriors
        f = plt.figure()
        gs0 = gridspec.GridSpec(3, 1, left=0.1, right=0.95, wspace=0.2, hspace=0.5)
        gs1 = gridspec.GridSpecFromSubplotSpec(2,2, subplot_spec=gs0[:-1], wspace=0.0, hspace=0.0)
        
        axA = plt.Subplot(f, gs1[:,:])
        f.add_subplot(axA)

        gs2 = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=gs0[-1], wspace=0.25)
        ax1 = plt.Subplot(f, gs2[0])
        f.add_subplot(ax1)

        ax2 = plt.Subplot(f, gs2[1])
        f.add_subplot(ax2)

        axA.set_title('Particle System t = 0')
        axA.imshow(zz0, **kwargs)
        axA.set_xlabel(Parameters['param_to_fit'][0])
        axA.set_ylabel(Parameters['param_to_fit'][1])    
        axA.set_aspect('auto')     
        axA.tick_params(axis='both', which='major', labelsize=10)

        p1 = Parameters['param_to_fit']
        ax1.plot(sampling[0], y01, color='blue')
        ax1.set_xlabel(Parameters['param_to_fit'][0])
        ax1.set_ylabel('density', fontsize=8)
        ax1.tick_params(axis='both', which='major', labelsize=8)
        ax1.set_xlim(Parameters['prior'][p1]['min'], Parameters['prior'][p1]['max'])

        p2 = Parameters['param_to_fit']
        ax2.plot(sampling[1], y02, color='blue')
        ax2.set_xlabel(Parameters['param_to_fit'][1])
        ax2.set_ylabel('density', fontsize = 8)
        ax2.tick_params(axis='both', which='major', labelsize=8)
        ax2.set_xlim(Parameters['prior'][p2]['min'], Parameters['prior'][p2]['max']) 

        pdf.savefig()
        plt.close()

        print 'Finished plotting particle system T=0' 

        for i in range(T+1):
           
            op1 = open(Parameters['file_root'] + str( i ) + '.dat', 'r')
            lin1 = op1.readlines()
            op1.close()

            d = [elem.split() for elem in lin1]
            d1 = np.array([[float(item) for item in line] for line in d[1:]])
     
            epsilon_ev.append([float(d[1][d[0].index('dist_threshold' + str(jj + 1))])       
                              for jj in xrange(Parameters['dist_dim'])])

            time_ev.append(sum(float(line[d[0].index('time')]) for line in d1[1:]))
            ndraws_ev.append(sum(float(line[d[0].index('NDraws')]) for line in d1[1:]))

            if i > 0:
                w1 = np.loadtxt(Parameters['file_root'] + str(i) + 'weights.dat')
            else:
                w1 = np.array([1.0/len(d1) for k in range(len(d1))])

            kde1 = gaussian_kde(d1[:,0] , weights=w1)
            y1 = kde1(sampling[0])

            kde2 = gaussian_kde(d1[:,1], weights=w1)
            y2 = kde2(sampling[1])
 
            d2 = np.array(d1[:,:len(Parameters['param_to_fit'])]) 
            kde3 = gaussian_kde(d2.transpose(), weights=w1)
            y3 = kde3((np.ravel(xx), np.ravel(yy)))
            zz = np.reshape(y3, xx.shape)    

            #### Plot posteriors
            f = plt.figure()
            gs0 = gridspec.GridSpec(3, 1, left=0.1, right=0.95, wspace=0.2, hspace=0.5)
            gs1 = gridspec.GridSpecFromSubplotSpec(2,2, subplot_spec=gs0[:-1], wspace=0.0, hspace=0.0)
        
            axA = plt.Subplot(f, gs1[:,:])
            f.add_subplot(axA)

            gs2 = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=gs0[-1], wspace=0.25)
            ax1 = plt.Subplot(f, gs2[0])
            f.add_subplot(ax1)

            ax2 = plt.Subplot(f, gs2[1])
            f.add_subplot(ax2)

            axA.set_title('Particle System t = ' + str( i + 1))
            axA.imshow(zz, **kwargs)
            axA.set_xlabel(Parameters['param_to_fit'][0])
            axA.set_ylabel(Parameters['param_to_fit'][1])    
            axA.set_aspect('auto')     
            axA.tick_params(axis='both', which='major', labelsize=10)

            ax1.plot(sampling[0], y1, color='blue')
            ax1.set_xlabel(Parameters['param_to_fit'][0])
            ax1.set_ylabel('density', fontsize=8)
            ax1.tick_params(axis='both', which='major', labelsize=8)
            ax1.set_xlim(Parameters['prior'][p1]['min'], Parameters['prior'][p1]['max'])

            ax2.plot(sampling[1], y2, color='blue')
            ax2.set_xlabel(Parameters['param_to_fit'][1])
            ax2.set_ylabel('density', fontsize = 8)
            ax2.tick_params(axis='both', which='major', labelsize=8)
            ax2.set_xlim(Parameters['prior'][p2]['min'], Parameters['prior'][p2]['max']) 

            pdf.savefig()
            plt.close()
           
            print 'Finished plotting particle system T=' +str(i + 1)

    
        convergence = np.array([float(Parameters['M'])/item  for item in ndraws_ev])
    
        #Plot epsilon evolution
        plt.figure()
        plt.title('Distance threshold evolution')
        for jj in xrange(Parameters['dist_dim']):
            plt.scatter(range(1, T + 2), np.array(epsilon_ev)[:,jj]/max(np.array(epsilon_ev)[:,jj]), 
                        color=color[jj], marker=marker[jj], label='distance threshold' + str(jj + 1))
        plt.legend()
        plt.xlabel('Particle System')
        plt.ylabel(r'$\epsilon_{\rm norm}$')
        pdf.savefig()
        plt.close()

        plt.figure()
        plt.title('Computational time evolution')
        plt.scatter(range(1, T + 2), np.array(time_ev), 
                   color=color[Parameters['dist_dim']], 
                   marker=marker[Parameters['dist_dim']], label='time')
        plt.legend(loc='upper left')
        plt.xlabel('Particle System')
        plt.ylabel('time(s)')
        pdf.savefig()
        plt.close()

        plt.figure()
        plt.title('Convergene criteria evolution')
        plt.scatter(range(1, T + 2), convergence, 
                    color=color[Parameters['dist_dim'] + 1], 
                    marker=marker[Parameters['dist_dim'] + 1], label='convergence')
        plt.legend()
        plt.xlabel('Particle System')
        plt.ylabel('convergence criteria(M/K)')
        pdf.savefig()
        plt.close()


def plot_3p(T, file_output, Parameters):
    """
    Make 3 free parameters plot for ABC results. 

    input:  T 		-> number of total particle systems (int)  
	    file_output -> name of output figure file (str)
            Parameters 	-> input parameter dictionary - same as used for ABC sampler (dict)

    output: summary pdf file with multiple pages for each generated particle system. 

    *******
    If you want to change the range of parameters for plotting, change values for the key 'param_lim' in the input dictionary parameter! 
    *******
    """

    if file_output[-3:] != 'pdf':
        raise NameError('Name file for figure output must be a pdf!')       

    #sampling for 1D plots
    sampling = np.array([[i for i in np.arange(Parameters['prior'][par]['min'], 
                               Parameters['prior'][par]['max'], (Parameters['prior'][par]['max'] - 
                               float(Parameters['prior'][par]['min']))/1000)]  
                                     for par in Parameters['param_to_fit']])
    
    #sampling for 2D plots
    sampling2 = np.array([[i for i in np.arange(Parameters['prior'][par]['min'], 
                                Parameters['prior'][par]['max'], (Parameters['prior'][par]['max'] - 
                                Parameters['prior'][par]['min'])/100)] 
                                          for par in Parameters['param_to_fit']])

    #define variables for 2D plots 
    xx12, yy12 = np.meshgrid(sampling2[0], sampling2[1])
    xx13, yy13 = np.meshgrid(sampling2[0], sampling2[2])
    xx23, yy23 = np.meshgrid(sampling2[1], sampling2[2])
    
    p1 = Parameters['param_to_fit'][0]
    p2 = Parameters['param_to_fit'][1]
    p3 = Parameters['param_to_fit'][2]
   
    kwargs12 = dict(extent=(Parameters['prior'][p1]['min'], Parameters['prior'][p1]['max'], 
                            Parameters['prior'][p2]['min'], Parameters['prior'][p2]['max']), 
                            cmap='hot', origin='lower')
 
    kwargs13 = dict(extent=(Parameters['prior'][p1]['min'], Parameters['prior'][p1]['max'], 
                            Parameters['prior'][p3]['min'], Parameters['prior'][p3]['max']), 
                            cmap='hot', origin='lower')

    kwargs23 = dict(extent=(Parameters['prior'][p2]['min'], Parameters['prior'][p2]['max'], 
                            Parameters['prior'][p3]['min'], Parameters['prior'][p3]['max']), 
                            cmap='hot', origin='lower')

    y0 = np.array([[Parameters['prior'][par]['func'](Parameters['prior'][par]) for par in Parameters['param_to_fit']] 
                                 for x in xrange(Parameters['M'])])
    
    w0 = [1.0/len(y0) for i in y0]

    kde01 = gaussian_kde(y0[:,0] , weights=w0)
    y01 = kde01(sampling[0])

    kde02 = gaussian_kde(y0[:,1], weights=w0)
    y02 = kde02(sampling[1])

    kde03 = gaussian_kde(y0[:,2], weights=w0)
    y03 = kde03(sampling[2])

    kde012 = gaussian_kde(np.array([y0[:,0], y0[:,1]]), weights=w0)  
    y012 = kde012((np.ravel(xx12), np.ravel(yy12)))
    zz012 = np.reshape(y012, xx12.shape)   
    
    kde013 = gaussian_kde(np.array([y0[:,0], y0[:,2]]), weights=w0)   
    y013 = kde013((np.ravel(xx13), np.ravel(yy13)))
    zz013 = np.reshape(y013, xx13.shape)    
    
    kde023 = gaussian_kde(np.array([y0[:,1], y0[:,2]]), weights=w0)
    y023 = kde023((np.ravel(xx23), np.ravel(yy23)))
    zz023 = np.reshape(y023, xx23.shape)

    epsilon_ev = []
    time_ev = []
    ndraws_ev = []
   
    marker = ['*', '+', 'x', 'o', '^', 'd', '<']
    color = ['red', 'green', 'blue','brown', 'gray', 'purple']

    with PdfPages(file_output) as pdf:

        #### Plot posteriors
        f = plt.figure()
            
        gs0 = gridspec.GridSpec(2, 3, left=0.075, right=0.975, wspace=0.35, hspace=0.3)
            
        ax1 = plt.Subplot(f, gs0[0])
        f.add_subplot(ax1)

        ax2 = plt.Subplot(f, gs0[1])
        f.add_subplot(ax2)

        ax3 = plt.Subplot(f, gs0[2])
        f.add_subplot(ax3)

        ax4 = plt.Subplot(f, gs0[3])
        f.add_subplot(ax4)

        ax5 = plt.Subplot(f, gs0[4])
        f.add_subplot(ax5)

        ax6 = plt.Subplot(f, gs0[5])
        f.add_subplot(ax6)

        ax1.imshow(zz012, **kwargs12)
        ax1.set_xlabel(Parameters['param_to_fit'][0])
        ax1.set_ylabel(Parameters['param_to_fit'][1])    
        ax1.set_aspect('auto')     
        ax1.tick_params(axis='both', which='major', labelsize=8)
            
        ax2.set_title('Particle System t = 0')
        ax2.imshow(zz013, **kwargs13)
        ax2.set_xlabel(Parameters['param_to_fit'][0])
        ax2.set_ylabel(Parameters['param_to_fit'][2])    
        ax2.set_aspect('auto')     
        ax2.tick_params(axis='both', which='major', labelsize=8)
   
        ax3.imshow(zz023, **kwargs23)
        ax3.set_xlabel(Parameters['param_to_fit'][1])
        ax3.set_ylabel(Parameters['param_to_fit'][2])    
        ax3.set_aspect('auto')     
        ax3.tick_params(axis='both', which='major', labelsize=8) 

        ax4.plot(sampling[0], y01, color='blue')
        ax4.set_xlabel(Parameters['param_to_fit'][0])
        ax4.set_ylabel('density', fontsize=8)
        ax4.tick_params(axis='both', which='major', labelsize=8)
        ax4.set_xlim(Parameters['prior'][p1]['min'], Parameters['prior'][p1]['max'])

        ax5.plot(sampling[1], y02, color='red')
        ax5.set_xlabel(Parameters['param_to_fit'][1])
        ax5.set_ylabel('density', fontsize = 8)
        ax5.tick_params(axis='both', which='major', labelsize=8)
        ax5.set_xlim(Parameters['prior'][p2]['min'], Parameters['prior'][p2]['max']) 

        ax6.plot(sampling[2], y03, color='green')
        ax6.set_xlabel(Parameters['param_to_fit'][2])
        ax6.set_ylabel('density', fontsize = 8)
        ax6.tick_params(axis='both', which='major', labelsize=8)
        ax6.set_xlim(Parameters['prior'][p3]['min'], Parameters['prior'][p3]['max']) 

        pdf.savefig()
        plt.close()

        print 'Finished plotting particle system T=0'

        for i in range(T+1):
           
            #read individual particle systems
            op1 = open(Parameters['file_root'] + str(i) + '.dat', 'r')
            lin1 = op1.readlines()
            op1.close()

            d = [elem.split() for elem in lin1]
            d1 = np.array([[float(item) for item in line] for line in d[1:]])
     
            epsilon_ev.append([float(d[1][d[0].index('dist_threshold' + str(jj + 1))])  
                              for jj in xrange(Parameters['dist_dim'])])

            time_ev.append(sum(float(line[ d[0].index('time')]) for line in d1[1:]))
            ndraws_ev.append(sum(float(line[ d[0].index('NDraws')]) for line in d1[1:]))

            if i > 0:
                w1 = np.loadtxt(Parameters['file_root'] + str(i) + 'weights.dat')
            else:
                w1 = np.array([1.0/len(d1) for k in range(len(d1))])

            kde1 = gaussian_kde(d1[:,0] , weights=w1)
            y1 = kde1(sampling[0])

            kde2 = gaussian_kde(d1[:,1], weights=w1)
            y2 = kde2(sampling[1])

            kde3 = gaussian_kde(d1[:,2], weights=w1)
            y3 = kde3(sampling[2]) 
 
            kde12 = gaussian_kde(np.array([d1[:,0], d1[:,1]]), weights=w1)
            xx12, yy12 = np.meshgrid(sampling2[0], sampling2[1])
            y12 = kde12((np.ravel(xx12), np.ravel(yy12)))
            zz12 = np.reshape(y12, xx12.shape)
   
            kwargs12 = dict(extent=(Parameters['prior'][p1]['min'], Parameters['prior'][p1]['max'], 
                                    Parameters['prior'][p2]['min'], Parameters['prior'][p2]['max']), 
                                    cmap='hot', origin='lower') 

            kde13 = gaussian_kde(np.array([d1[:,0], d1[:,2]]), weights=w1)
            xx13, yy13 = np.meshgrid(sampling2[0], sampling2[2])
            y13 = kde13(( np.ravel(xx13), np.ravel(yy13)))
            zz13 = np.reshape( y13, xx13.shape)    

            kwargs13 = dict(extent=(Parameters['prior'][p1]['min'], Parameters['prior'][p1]['max'], 
                                    Parameters['prior'][p3]['min'], Parameters['prior'][p3]['max']), 
                                    cmap='hot', origin='lower')

            kde23 = gaussian_kde(np.array([d1[:,1], d1[:,2]]), weights=w1)
            xx23, yy23 = np.meshgrid(sampling2[1], sampling2[2])
            y23 = kde23((np.ravel(xx23), np.ravel(yy23)))
            zz23 = np.reshape(y23, xx23.shape)

            kwargs23 = dict(extent=(Parameters['prior'][p2]['min'], Parameters['prior'][p2]['max'], 
                                    Parameters['prior'][p3]['min'], Parameters['prior'][p3]['max']), 
                                    cmap='hot', origin='lower')

            #### Plot posteriors
            f = plt.figure()
          
            gs0 = gridspec.GridSpec(2, 3, left=0.075, right=0.975, wspace=0.35, hspace=0.3)
        
            ax1 = plt.Subplot(f, gs0[0])
            f.add_subplot(ax1)

            ax2 = plt.Subplot(f, gs0[1])
            f.add_subplot(ax2)

            ax3 = plt.Subplot(f, gs0[2])
            f.add_subplot(ax3)

            ax4 = plt.Subplot(f, gs0[3])
            f.add_subplot(ax4)

            ax5 = plt.Subplot(f, gs0[4])
            f.add_subplot(ax5)

            ax6 = plt.Subplot(f, gs0[5])
            f.add_subplot(ax6)

            ax1.imshow(zz12, **kwargs12)
            ax1.set_xlabel(Parameters['param_to_fit'][0])
            ax1.set_ylabel(Parameters['param_to_fit'][1])    
            ax1.set_aspect('auto')     
            ax1.tick_params(axis='both', which='major', labelsize=8)

            ax2.set_title('Particle System t = ' + str(i + 1))
            ax2.imshow(zz13, **kwargs13)
            ax2.set_xlabel(Parameters['param_to_fit'][0])
            ax2.set_ylabel(Parameters['param_to_fit'][2])    
            ax2.set_aspect('auto')     
            ax2.tick_params(axis='both', which='major', labelsize=8)

            ax3.imshow(zz23, **kwargs23)
            ax3.set_xlabel(Parameters['param_to_fit'][1])
            ax3.set_ylabel(Parameters['param_to_fit'][2])    
            ax3.set_aspect('auto')     
            ax3.tick_params(axis='both', which='major', labelsize=8) 

            ax4.plot(sampling[0], y1, color='blue')
            ax4.set_xlabel(Parameters['param_to_fit'][0])
            ax4.set_ylabel('density', fontsize=8)
            ax4.tick_params(axis='both', which='major', labelsize=8)
            ax4.set_xlim(Parameters['prior'][p1]['min'], Parameters['prior'][p1]['max'])

            ax5.plot(sampling[1], y2, color='red')
            ax5.set_xlabel(Parameters['param_to_fit'][1])
            ax5.set_ylabel('density', fontsize = 8)
            ax5.tick_params(axis='both', which='major', labelsize=8)
            ax5.set_xlim(Parameters['prior'][p2]['min'], Parameters['prior'][p2]['max']) 

            ax6.plot(sampling[2], y3, color='green')
            ax6.set_xlabel(Parameters['param_to_fit'][2])
            ax6.set_ylabel('density', fontsize = 8)
            ax6.tick_params(axis='both', which='major', labelsize=8)
            ax6.set_xlim(Parameters['prior'][p3]['min'], Parameters['prior'][p3]['max']) 

            pdf.savefig()
            plt.close()

            print 'Finished plotting particle system T=' + str(i + 1)
    
        convergence = np.array([float(Parameters['M'])/item  for item in ndraws_ev])
    
        #Plot epsilon evolution
        plt.figure()
        for jj in xrange(Parameters['dist_dim']):
            plt.scatter(range(1,T + 2), np.array(epsilon_ev)[:, jj]/max(np.array(epsilon_ev)[:,jj]), 
                        color=color[jj], marker=marker[jj], label='distance threshold' + str(jj + 1))
        plt.legend()
        plt.xlabel('Particle System')
        plt.ylabel(r'$\epsilon_{\rm norm}$')
        pdf.savefig()
        plt.close()

        plt.figure()
        plt.scatter(range(1, T+2 ), np.array(time_ev), 
                    color=color[Parameters['dist_dim']], 
                    marker=marker[Parameters['dist_dim']], label='time')
        plt.legend(loc='upper left')
        plt.xlabel('Particle System')
        plt.ylabel('time(s)')
        pdf.savefig()
        plt.close()

        plt.figure()
        plt.scatter(range(1, T + 2), convergence, 
                    color=color[Parameters['dist_dim'] + 1], 
                    marker=marker[Parameters['dist_dim'] + 1], label='convergence')
        plt.legend()
        plt.xlabel('Particle System')
        plt.ylabel('convergence criteria(M/K)')
        pdf.savefig()
        plt.close()


def plot_4p(T, file_output, Parameters):
    """
    Make 4 free parameter plot for ABC results. 

    input:  T 		-> number of total particle systems (int)  
	    file_output -> name of output figure file (str)
            Parameters 	-> input parameter dictionary - same as used for ABC sampler (dict)

    output: summary pdf file with multiple pages for each generated particle system. 

    *******
    If you want to change the range of parameters for plotting, change values for the key 'param_lim' in the input dictionary parameter! 
    *******
    """

    if file_output[-3:] != 'pdf':
        raise NameError('Name file for figure output must be a pdf!')       


    #sampling for 1D plots
    sampling = np.array([[i for i in np.arange( Parameters['prior'][par]['min'], 
                               Parameters['prior'][par]['max'], (Parameters['prior'][par]['max'] - 
                               Parameters['prior'][par]['min'])/1000)] for par in Parameters['param_to_fit']])
   
    #sampling for 2D plots
    sampling2 = np.array([[i for i in np.arange(Parameters['prior'][par]['min'], 
                                Parameters['prior'][par]['max'], (Parameters['prior'][par]['max'] - 
                                Parameters['prior'][par]['min'])/100)] for par in Parameters['param_to_fit']])

    #define variables for 2D plots 
    xx12, yy12 = np.meshgrid(sampling2[0], sampling2[1])
    xx13, yy13 = np.meshgrid(sampling2[0], sampling2[2])
    xx14, yy14 = np.meshgrid(sampling2[0], sampling2[3])
    xx23, yy23 = np.meshgrid(sampling2[1], sampling2[2])
    xx24, yy24 = np.meshgrid(sampling2[1], sampling2[3])
    xx34, yy34 = np.meshgrid(sampling2[2], sampling2[3])

    p1 = Parameters['param_to_fit'][0]
    p2 = Parameters['param_to_fit'][1]
    p3 = Parameters['param_to_fit'][2]
    p4 = Parameters['param_to_fit'][3]
    
    kwargs12 = dict(extent=(Parameters['prior'][p1]['min'], Parameters['prior'][p1]['max'], 
                            Parameters['prior'][p2]['min'], Parameters['prior'][p2]['max']), 
                            cmap='hot', origin='lower')
 
    kwargs13 = dict(extent=(Parameters['prior'][p1]['min'], Parameters['prior'][p1]['max'], 
                            Parameters['prior'][p3]['min'], Parameters['prior'][p3]['max']), 
                            cmap='hot', origin='lower')

    kwargs14 = dict(extent=(Parameters['prior'][p1]['min'], Parameters['prior'][p1]['max'], 
                            Parameters['prior'][p4]['min'], Parameters['prior'][p4]['max']), 
                            cmap='hot', origin='lower')

    kwargs23 = dict(extent=(Parameters['prior'][p2]['min'], Parameters['prior'][p2]['max'], 
                            Parameters['prior'][p3]['min'], Parameters['prior'][p3]['max']), 
                            cmap='hot', origin='lower')

    kwargs24 = dict(extent=(Parameters['prior'][p2]['min'], Parameters['prior'][p2]['max'], 
                            Parameters['prior'][p4]['min'], Parameters['prior'][p4]['max']), 
                            cmap='hot', origin='lower')

    kwargs34 = dict(extent=(Parameters['prior'][p3]['min'], Parameters['prior'][p3]['max'], 
                            Parameters['prior'][p4]['min'], Parameters['prior'][p4]['max']), 
                            cmap='hot', origin='lower')

    y0 = np.array([[Parameters['prior'][par]['func'](Parameters['prior'][par]) for par in Parameters['param_to_fit']] for x in xrange(Parameters['M'])])
    
    w0 = [1.0/len(y0) for i in y0]

    kde01 = gaussian_kde(y0[:,0] , weights=w0)
    y01 = kde01(sampling[0])

    kde02 = gaussian_kde(y0[:,1], weights=w0)
    y02 = kde02(sampling[1])

    kde03 = gaussian_kde(y0[:,2], weights=w0)
    y03 = kde03(sampling[2])

    kde04 = gaussian_kde(y0[:,3], weights=w0)
    y04 = kde04(sampling[3])

    kde012 = gaussian_kde(np.array([y0[:,0], y0[:,1]]), weights=w0)  
    y012 = kde012((np.ravel(xx12), np.ravel(yy12)))
    zz012 = np.reshape(y012, xx12.shape)   
    
    kde013 = gaussian_kde(np.array([y0[:,0], y0[:,2]]), weights=w0)   
    y013 = kde013((np.ravel(xx13), np.ravel(yy13)))
    zz013 = np.reshape(y013, xx13.shape)  

    kde014 = gaussian_kde(np.array([y0[:,0], y0[:,3]]), weights=w0)   
    y014 = kde014((np.ravel(xx14), np.ravel(yy14)))
    zz014 = np.reshape(y014, xx14.shape)   
    
    kde023 = gaussian_kde(np.array([y0[:,1], y0[:,2]]), weights=w0)
    y023 = kde023((np.ravel(xx23), np.ravel(yy23)))
    zz023 = np.reshape(y023, xx23.shape)

    kde024 = gaussian_kde(np.array([y0[:,1], y0[:,3]]), weights=w0)
    y024 = kde024((np.ravel(xx24), np.ravel(yy24)))
    zz024 = np.reshape(y024, xx24.shape)

    kde034 = gaussian_kde(np.array([y0[:,2], y0[:,3]]), weights=w0)
    y034 = kde034((np.ravel(xx34), np.ravel(yy34)))
    zz034 = np.reshape(y034, xx34.shape)

    epsilon_ev = []
    time_ev = []
    ndraws_ev = []
   
    marker = ['*', '+', 'x', 'o', '^', 'd', '<']
    color = ['red', 'green', 'blue', 'brown', 'gray', 'purple']

    with PdfPages(file_output) as pdf:

        #### Plot posteriors
        f = plt.figure()           
            
        gs0 = gridspec.GridSpec(3, 4, left=0.075, right=0.975, wspace=0.35, hspace=0.3)
            
        ax1 = plt.Subplot(f, gs0[0])
        f.add_subplot(ax1)

        ax2 = plt.Subplot(f, gs0[1])
        f.add_subplot(ax2)

        ax3 = plt.Subplot(f, gs0[2])
        f.add_subplot(ax3)

        ax4 = plt.Subplot(f, gs0[3])
        f.add_subplot(ax4)

        ax5 = plt.Subplot(f, gs0[4])
        f.add_subplot(ax5)

        ax6 = plt.Subplot(f, gs0[5])
        f.add_subplot(ax6)

        ax7 = plt.Subplot(f, gs0[6])
        f.add_subplot(ax7)

        ax8 = plt.Subplot(f, gs0[7])
        f.add_subplot(ax8)

        ax9 = plt.Subplot(f, gs0[8])
        f.add_subplot(ax9)

        ax10 = plt.Subplot(f, gs0[9])
        f.add_subplot(ax10)

        ax1.imshow(zz012, **kwargs12)
        ax1.set_xlabel(Parameters['param_to_fit'][0])
        ax1.set_ylabel(Parameters['param_to_fit'][1])    
        ax1.set_aspect('auto')     
        ax1.tick_params(axis='both', which='major', labelsize=8)
            
        ax2.set_title('Particle System t = 0')
        ax2.imshow(zz013, **kwargs13)
        ax2.set_xlabel(Parameters['param_to_fit'][0])
        ax2.set_ylabel(Parameters['param_to_fit'][2])    
        ax2.set_aspect('auto')     
        ax2.tick_params(axis='both', which='major', labelsize=8)
   
        ax3.imshow(zz014, **kwargs14)
        ax3.set_xlabel(Parameters['param_to_fit'][0])
        ax3.set_ylabel(Parameters['param_to_fit'][3])    
        ax3.set_aspect('auto')     
        ax3.tick_params(axis='both', which='major', labelsize=8) 

        ax4.imshow(zz023, **kwargs23)
        ax4.set_xlabel(Parameters['param_to_fit'][1])
        ax4.set_ylabel(Parameters['param_to_fit'][2])    
        ax4.set_aspect('auto')     
        ax4.tick_params(axis='both', which='major', labelsize=8) 

        ax5.imshow(zz024, **kwargs24)
        ax5.set_xlabel(Parameters['param_to_fit'][1])
        ax5.set_ylabel(Parameters['param_to_fit'][3])    
        ax5.set_aspect('auto')     
        ax5.tick_params(axis='both', which='major', labelsize=8)  

        ax6.imshow(zz034, **kwargs34)
        ax6.set_xlabel(Parameters['param_to_fit'][2])
        ax6.set_ylabel(Parameters['param_to_fit'][3])    
        ax6.set_aspect('auto')     
        ax6.tick_params(axis='both', which='major', labelsize=8) 

        ax7.plot(sampling[0], y01, color='blue')
        ax7.set_xlabel(Parameters['param_to_fit'][0])
        ax7.set_ylabel('density', fontsize=8)
        ax7.tick_params(axis='both', which='major', labelsize=8)
        ax7.set_xlim(Parameters['prior'][p1]['min'], Parameters['prior'][p1]['max'])

        ax8.plot(sampling[1], y02, color='red')
        ax8.set_xlabel(Parameters['param_to_fit'][1])
        ax8.set_ylabel('density', fontsize = 8)
        ax8.tick_params(axis='both', which='major', labelsize=8)
        ax8.set_xlim(Parameters['prior'][p2]['min'], Parameters['prior'][p2]['max']) 

        ax9.plot(sampling[2], y03, color='green')
        ax9.set_xlabel(Parameters['param_to_fit'][2])
        ax9.set_ylabel('density', fontsize = 8)
        ax9.tick_params(axis='both', which='major', labelsize=8)
        ax9.set_xlim(Parameters['prior'][p3]['min'], Parameters['prior'][p3]['max']) 

        ax10.plot(sampling[3], y04, color='purple')
        ax10.set_xlabel(Parameters['param_to_fit'][2])
        ax10.set_ylabel('density', fontsize = 8)
        ax10.tick_params(axis='both', which='major', labelsize=8)
        ax10.set_xlim(Parameters['prior'][p4]['min'], Parameters['prior'][p4]['max'])

        pdf.savefig()
        plt.close()

        print 'Finished plotting particle system T=0'

        for i in range(T+1):
           
            #read individual particle systems
            op1 = open(Parameters['file_root'] + str( i ) + '.dat', 'r')
            lin1 = op1.readlines()
            op1.close()


            d = [elem.split() for elem in lin1]
            d1 = np.array([[float(item) for item in line] for line in d[1:]])
     
            epsilon_ev.append([float(d[1][d[0].index('dist_threshold' + str(jj + 1))]) 
                              for jj in xrange(Parameters['dist_dim'])])

            time_ev.append(sum(float(line[d[0].index('time')]) for line in d1[1:]))
            ndraws_ev.append(sum(float(line[d[0].index('NDraws')]) for line in d1[1:]))

            if i > 0:
                w1 = np.loadtxt(Parameters['file_root'] + str(i) + 'weights.dat')
            else:
                w1 = np.array([1.0/len(d1) for k in range(len(d1))])

            kde1 = gaussian_kde(d1[:,0] , weights=w1)
            y1 = kde1(sampling[0])

            kde2 = gaussian_kde(d1[:,1], weights=w1)
            y2 = kde2(sampling[1])

            kde3 = gaussian_kde(d1[:,2], weights=w1)
            y3 = kde3(sampling[2]) 
 
            kde4 = gaussian_kde(d1[:,3], weights=w1)
            y4 = kde3(sampling[3]) 

            kde12 = gaussian_kde(np.array([d1[:,0], d1[:,1]]), weights=w1)
            xx12, yy12 = np.meshgrid(sampling2[0], sampling2[1])
            y12 = kde12((np.ravel(xx12), np.ravel(yy12)))
            zz12 = np.reshape(y12, xx12.shape)   

            kwargs12 = dict(extent=(Parameters['prior'][p1]['min'], Parameters['prior'][p1]['max'], 
                                    Parameters['prior'][p2]['min'], Parameters['prior'][p2]['max']), 
                                    cmap='hot', origin='lower') 

            kde13 = gaussian_kde(np.array([d1[:,0], d1[:,2]]), weights=w1)
            xx13, yy13 = np.meshgrid(sampling2[0], sampling2[2])
            y13 = kde13((np.ravel(xx13), np.ravel(yy13)))
            zz13 = np.reshape(y13, xx13.shape)
    
            kwargs13 = dict(extent=(Parameters['prior'][p1]['min'], Parameters['prior'][p1]['max'], 
                                    Parameters['prior'][p3]['min'], Parameters['prior'][p3]['max']), 
                                    cmap='hot', origin='lower')

            kde14 = gaussian_kde(np.array([d1[:,0], d1[:,3]]), weights=w1)
            xx14, yy14 = np.meshgrid(sampling2[0], sampling2[3])
            y14 = kde14((np.ravel(xx14), np.ravel(yy14)))
            zz14 = np.reshape(y14, xx14.shape)    

            kwargs14 = dict(extent=(Parameters['prior'][p1]['min'], Parameters['prior'][p1]['max'], 
                                    Parameters['prior'][p4]['min'], Parameters['prior'][p4]['max']), 
                                    cmap='hot', origin='lower') 

            kde23 = gaussian_kde(np.array([d1[:,1], d1[:,2]]), weights=w1)
            xx23, yy23 = np.meshgrid(sampling2[1], sampling2[2])
            y23 = kde23((np.ravel(xx23), np.ravel(yy23)))
            zz23 = np.reshape(y23, xx23.shape)

            kwargs23 = dict(extent=(Parameters['prior'][p2]['min'], Parameters['prior'][p2]['max'], 
                                    Parameters['prior'][p3]['min'], Parameters['prior'][p3]['max']), 
                                    cmap='hot', origin='lower')

            kde24 = gaussian_kde(np.array([d1[:,1], d1[:,3] ]), weights=w1)
            xx24, yy24 = np.meshgrid(sampling2[1], sampling2[3])
            y24 = kde24((np.ravel(xx24), np.ravel(yy24)))
            zz24 = np.reshape(y24, xx24.shape)

            kwargs24 = dict(extent=(Parameters['prior'][p2]['min'], Parameters['prior'][p2]['max'], 
                                    Parameters['prior'][p4]['min'], Parameters['prior'][p4]['max']), 
                                    cmap='hot', origin='lower')

            kde34 = gaussian_kde(np.array([d1[:,2], d1[:,3]]), weights=w1)
            xx34, yy34 = np.meshgrid(sampling2[2], sampling2[3])
            y34 = kde34((np.ravel(xx34), np.ravel(yy34)))
            zz34 = np.reshape(y34, xx34.shape)

            kwargs34 = dict(extent=(Parameters['prior'][p3]['min'], Parameters['prior'][p3]['max'], 
                                    Parameters['prior'][p4]['min'], Parameters['prior'][p4]['max']), 
                                    cmap='hot', origin='lower')

            #### Plot posteriors
            f = plt.figure()
                       
            gs0 = gridspec.GridSpec(3, 4, left=0.075, right=0.975, wspace=0.35, hspace=0.3)
        
            ax1 = plt.Subplot(f, gs0[0])
            f.add_subplot(ax1)

            ax2 = plt.Subplot(f, gs0[1])
            f.add_subplot(ax2)

            ax3 = plt.Subplot(f, gs0[2])
            f.add_subplot(ax3)

            ax4 = plt.Subplot(f, gs0[3])
            f.add_subplot(ax4)

            ax5 = plt.Subplot(f, gs0[4])
            f.add_subplot(ax5)

            ax6 = plt.Subplot(f, gs0[5])
            f.add_subplot(ax6)

            ax7 = plt.Subplot(f, gs0[6])
            f.add_subplot(ax7)

            ax8 = plt.Subplot(f, gs0[7])
            f.add_subplot(ax8)

            ax9 = plt.Subplot(f, gs0[8])
            f.add_subplot(ax9)

            ax10 = plt.Subplot(f, gs0[9])
            f.add_subplot(ax10)
   
            ax1.imshow(zz12, **kwargs12)
            ax1.set_xlabel(Parameters['param_to_fit'][0])
            ax1.set_ylabel(Parameters['param_to_fit'][1])    
            ax1.set_aspect('auto')     
            ax1.tick_params(axis='both', which='major', labelsize=8)
            
            ax2.set_title('Particle System t = ' + str(i + 1))
            ax2.imshow(zz13, **kwargs13)
            ax2.set_xlabel(Parameters['param_to_fit'][0])
            ax2.set_ylabel(Parameters['param_to_fit'][2])    
            ax2.set_aspect('auto')     
            ax2.tick_params(axis='both', which='major', labelsize=8)

            ax3.set_title('Particle System t = ' + str(i + 1))
            ax3.imshow(zz14, **kwargs14)
            ax3.set_xlabel(Parameters['param_to_fit'][0])
            ax3.set_ylabel(Parameters['param_to_fit'][3])    
            ax3.set_aspect('auto')     
            ax3.tick_params(axis='both', which='major', labelsize=8)

            ax4.imshow(zz23, **kwargs23)
            ax4.set_xlabel(Parameters['param_to_fit'][1])
            ax4.set_ylabel(Parameters['param_to_fit'][2])    
            ax4.set_aspect('auto')     
            ax4.tick_params(axis='both', which='major', labelsize=8) 

            ax5.imshow(zz24, **kwargs24)
            ax5.set_xlabel(Parameters['param_to_fit'][1])
            ax5.set_ylabel(Parameters['param_to_fit'][3])    
            ax5.set_aspect('auto')     
            ax5.tick_params(axis='both', which='major', labelsize=8) 

            ax6.imshow(zz34, **kwargs34)
            ax6.set_xlabel(Parameters['param_to_fit'][2] )
            ax6.set_ylabel(Parameters['param_to_fit'][3] )    
            ax6.set_aspect('auto')     
            ax6.tick_params(axis='both', which='major', labelsize=8)

            ax7.plot(sampling[0], y1, color='blue')
            ax7.set_xlabel(Parameters['param_to_fit'][0])
            ax7.set_ylabel('density', fontsize=8 )
            ax7.tick_params(axis='both', which='major', labelsize=8)
            ax7.set_xlim(Parameters['prior'][p1]['min'], Parameters['prior'][p1]['max'])


            ax8.plot(sampling[1], y2, color='red')
            ax8.set_xlabel(Parameters['param_to_fit'][1])
            ax8.set_ylabel('density', fontsize = 8)
            ax8.tick_params(axis='both', which='major', labelsize=8)
            ax8.set_xlim(Parameters['prior'][p2]['min'], Parameters['prior'][p2]['max']) 

            ax9.plot(sampling[2], y3, color='green')
            ax9.set_xlabel(Parameters['param_to_fit'][2])
            ax9.set_ylabel('density', fontsize = 8)
            ax9.tick_params(axis='both', which='major', labelsize=8)
            ax9.set_xlim(Parameters['prior'][p3]['min'], Parameters['prior'][p3]['max']) 
   
            ax10.plot(sampling[3], y3, color='purple')
            ax10.set_xlabel(Parameters['param_to_fit'][3] )
            ax10.set_ylabel('density', fontsize = 8)
            ax10.tick_params(axis='both', which='major', labelsize=8)
            ax10.set_xlim(Parameters['prior'][p4]['min'], Parameters['prior'][p4]['max']) 

            pdf.savefig()
            plt.close()

            print 'Finished plotting particle system T=' + str(i + 1)
    
        convergence = np.array([float(Parameters['M'])/item  for item in ndraws_ev])
    
        #Plot epsilon evolution
        plt.figure()
        for jj in xrange(Parameters['dist_dim']):
            plt.scatter(range(1, T + 2), np.array(epsilon_ev)[:,jj]/max(np.array(epsilon_ev)[:,jj]), 
                        color=color[jj], marker=marker[jj], label='distance threshold' + str(jj + 1))
        plt.legend()
        plt.xlabel('Particle System')
        plt.ylabel(r'$\epsilon_{\rm norm}$')
        pdf.savefig()
        plt.close()

        plt.figure()
        plt.scatter(range(1, T + 2), np.array(time_ev), 
                     color=color[Parameters['dist_dim']], 
                     marker=marker[Parameters['dist_dim']], label='time')
        plt.legend(loc='upper left')
        plt.xlabel('Particle System')
        plt.ylabel('time(s)')
        pdf.savefig()
        plt.close()

        plt.figure()
        plt.scatter(range(1, T + 2), convergence, color=color[Parameters['dist_dim'] + 1], 
                    marker=marker[Parameters['dist_dim'] + 1], label='convergence')
        plt.legend()
        plt.xlabel('Particle System')
        plt.ylabel('convergence criteria(M/K)')
        pdf.savefig()
        plt.close()
