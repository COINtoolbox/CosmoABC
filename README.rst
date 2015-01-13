CosmoABC - Likelihood free parameter estimation for cosmology
**************************************************************


CosmoABC is a package which enables parameter inference using an Approximate Bayesian Computation (ABC) algorithm, as described in Ishida et al., 2015 [LINK].
The code was originally designed for cosmological parameter inference from galaxy clusters number counts based on Sunyaev-Zel'dovich measurements. In this context, the cosmological simulations were performed using the NumCosmo library.

Nevertheless, the user can easily take advantadge of the ABC sampler along with his/her own simulator, as well as  test personalized summary statistics and distance functions. 


.. _examples:

Examples
========

Sample input in can be found in ~CosmoABC/examples. All example files mentioned in this section are available in that directory. 

The user input file should contain all necessary variables for simulation as well as for the ABC sampler.

.. literalinclude:: /home/emille/Dropbox/WGC/ABC/github/ABC_repository/CosmoABC/examples/example_user_input.dat




User defined simulation, distance and prior functions
-----------------------------------------------------

The most important ingredients in an ABC analysis are:

* the prior probability distributions (PDF)
* the simulator
* the distance function


CosmoABC is able to handdle user defined functions for all three elements. 
You will find example files in the corresponding directory which will help you taylor your functions for the ABC sampler. 

Built-in options for priors PDF are:

* Gaussian
* flat
* beta






NumCosmo simulations
--------------------

In order to reproduce the results of Ishida *et al.* 2015, first you need to make sure the NumCosmo library is running smoothly. 
Instructions for complete instalation and tests can be found in [LINK].






Once the simulator is installed run the complete ABC sampler + NumCosmo cluster simulations from the command line ::

    run_ABC_NumCosmo.py -i **example_NumCosmo_input.dat**

This will run the complete analysis presented in Ishida *et al.*, 2015 as well as produce
plots with the corresponding results.

** WARNING**
    
    *This might take a while! *




Documentation
=============

The complete documentation can be found in [LINK].


Requirements
============

* Python 2.7
* numpy >=1.8.2
* scipy >= 0.14.0
* statsmodels >= 0.5.0
* matplotlib >= 1.3.1     
* argparse >= 1.1
* imp
* math
* argparse


Optional
--------

* NumCosmo  [LINK]


License
=======

* GNU General Public License (GPL>=3)
