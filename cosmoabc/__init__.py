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

__author__ = ("E. E. O. Ishida, S. D. P. Vitenti, M. Penna-Lima," +
             "J. Cisewski, R. S. de Souza, A. M. M. Trindade, " + 
             "V. C. Busti, E. Cameron")
__maintainer__ = "E. E. O. Ishida"
__copyright__ = "Copyright 2015"
__version__ = "1.0.10"
__email__ = "emille@cosmostatistics-initiative.org"
__status__ = "Prototype"
__license__ = "GPL3"

import numpy as np


from .ABC_functions import SelectParamInnerLoop, SetDistanceFromSimulation, DrawAllParams
from .ABC_sampler import ABC
from .distances import distance_quantiles, distance_GRBF
from .priors import flat_prior, gaussian_prior, beta_prior
from .plots import plot_1p, plot_2p, plot_3p


def __cite__():

    print('@ARTICLE{2015A&C....13....1I,')
    print('author = {{Ishida}, E.~E.~O. and {Vitenti}, S.~D.~P. and {Penna-Lima}, M. and')
    print('         {Cisewski}, J. and {de Souza}, R.~S. and {Trindade}, A.~M.~M. and')
    print('         {Cameron}, E. and {Busti}, V.~C.},')
    print('title = "{COSMOABC: Likelihood-free inference via Population Monte Carlo Approximate Bayesian Computation}",')
    print('journal = {Astronomy and Computing},')
    print('archivePrefix = "arXiv",')
    print('eprint = {1504.06129},')
    print('keywords = {Galaxies: statistics, (cosmology:) large-scale structure of universe},')
    print('year = 2015,')
    print('month = nov,')
    print('volume = 13,')
    print('pages = {1-11},')
    print('doi = {10.1016/j.ascom.2015.09.001},')
    print('adsurl = {http://adsabs.harvard.edu/abs/2015A%26C....13....1I},')
    print('adsnote = {Provided by the SAO/NASA Astrophysics Data System}')
    print('}')
            






