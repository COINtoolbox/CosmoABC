from setuptools import setup
import cosmoabc.__init__ as cosmoabc

def readme():
    with open('README.rst') as f:
        return f.read()

setup(name='cosmoabc',
      version=cosmoabc.__version__,
      description='Python ABC sampler',
      long_description=readme(),
      long_description_content_type='text/markdown',
      url='https://github.com/COINtoolbox/cosmoabc',
      author=cosmoabc.__author__,
      author_email=cosmoabc.__email__,
      license='GNU Public License',
      packages=['cosmoabc'],
      install_requires=['numpy>=1.8.2',
                      'scipy>=1.7.0',
                      'statsmodels>=0.12.2',
                      'matplotlib>=3.4.2',
                      'distribute',
                      'datetime'               
      ],
      scripts=['cosmoabc/bin/run_ABC.py', 'cosmoabc/bin/run_ABC_NumCosmo.py',
               'cosmoabc/bin/continue_ABC.py', 'cosmoabc/bin/continue_ABC_NumCosmo.py',
               'cosmoabc/bin/plot_ABC.py', 'cosmoabc/tests/test_ABC_algorithm.py', 
               'cosmoabc/tests/test_ABC_distance.py'],
      test_suite='nose.collector',
      tests_require=['nose', 'nose-cover3'],
      include_package_data=True,
      package_dir= {'cosmoabc': 'cosmoabc', 'data': 'cosmoabc/data', 
                    'examples':'cosmoabc/examples'},
      package_data = {'cosmoabc/data':['SPT_sample.dat']},
      zip_safe=False,
      classifiers = [
        'Programming Language :: Python',
        'Natural Language :: English',
        'Environment :: X11 Applications',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU General Public License (GPL)',
        'Topic :: Software Development :: Libraries :: Python Modules',
        'Topic :: Scientific/Engineering :: Astronomy',
        ])
