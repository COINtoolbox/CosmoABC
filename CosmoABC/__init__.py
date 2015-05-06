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

from ABC_functions import SelectParamInnerLoop, SetDistanceFromSimulation, DrawAllParams
from ABC_sampler import ABC
from distances import distance_quantiles, distance_GRBF
from priors import flat_prior, gaussian_prior, beta_prior


def __cite__():

    print '@ARTICLE{2015arXiv150406129I,'
    print 'author = {{Ishida}, E.~E.~O. and {Vitenti}, S.~D.~P. and {Penna-Lima}, M. and'
    print '          {Cisewski}, J. and {de Souza}, R.~S. and {Trindade}, A.~M.~M. and'
    print '          {Cameron}, E. and {V.~C.~Busti}},'
    print '          title = "{cosmoabc: Likelihood-free inference via Population Monte Carlo Approximate Bayesian Computation}",'
    print 'journal = {ArXiv e-prints},'
    print 'archivePrefix = "arXiv",'
    print 'eprint = {1504.06129},'
    print 'keywords = {Astrophysics - Cosmology and Nongalactic Astrophysics, Astrophysics - Instrumentation and Methods for Astrophysics},'
    print 'year = 2015,'
    print 'month = apr,'
    print 'adsurl = {http://adsabs.harvard.edu/abs/2015arXiv150406129I},'
    print 'adsnote = {Provided by the SAO/NASA Astrophysics Data System}'
    print '}'
            






