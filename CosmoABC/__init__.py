"""
CosmoABC: Likelihood-free parameter estimation for cosmology
============================================================

.. module:: useful_1
   :platform: Unix, Mac
   :synopsis: An Approximate Bayesian Computation sampler.

.. moduleauthor:: Emille Ishida <emilleishida@gmail.com>

Documentation is available in the docstrings and online in the corresponding 
`GitHub <https://github.com/COINtoolbox/CosmoABC>`_   and 
`Read the Docs <http://cosmoabc.readthedocs.org/en/latest/>`_


Contents
--------
CosmoABC uses functionalities of a number of different libraries and
in addition also provides subpackages.


Subpackages
-----------
Using any of these subpackages requires an explicit import.  For example,
``import CosmoABC.priors``.

::

 ABC_sampler                  --- ABC acceptance/rejection algorithm
 priors                       --- prior probability distribution functions
 distances                    --- distance functions between 2 catalogs
 plots                        --- Plotting rotines
 sim_NumCosmo                 --- Simulation setup for NumCosmo library
 weighted_gaussian_kde        --- rotines for plotting

"""


