CosmoABC - Likelihood free parameter estimation for cosmology
**************************************************************


CosmoABC is a package which enables parameter inference using an Approximate Bayesian Computation (ABC) algorithm, as described in Ishida et al., 2015 [LINK].
The code was originally designed for cosmological parameter inference from galaxy clusters number counts based on Sunyaev-Zel'dovich measurements. In this context, the cosmological simulations were performed using the NumCosmo library [LINK].

Nevertheless, the user can easily use the ABC sampler along with his/her own simulation function, as well as  test personalized summary statistics and distance functions. 


Examples
========

Sample input in can be found in ~CosmoABC/examples.
To run the complete ABC sampler + NumCosmo cluster simulations from the command line ::

    run_ABC_NumCosmo.py -i example_NumCosmo_input.dat

.. warning::
    
    This might take a while!









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

Optional
--------

* NumCosmo


License
=======

* GNU General Public License (GPL>=3)
