from setuptools import setup
import CosmoABC.ABC_sampler as CosmoABC

def readme():
    with open('README.rst') as f:
        return f.read()

setup(name='CosmoABC',
      version=CosmoABC.__version__,
      description='Python ABC sampler',
      long_description=readme(),
      url='https://github.com/COINtoolbox/CosmoABC',
      author=CosmoABC.__author__,
      author_email=CosmoABC.__email__,
      license='GNU Public License',
      packages=['CosmoABC'],
      install_requires=[
                      'numpy>=1.8.2',
                      'scipy>=0.14.0',
                      'statsmodels>=0.5.0',
                      'matplotlib>=1.3.1',
                      'distribute',
                      'imp',
                      'math',
                      'datetime',
                      'time'                   
      ],
      scripts=['CosmoABC/bin/run_ABC.py', 'CosmoABC/bin/run_ABC_NumCosmo.py','CosmoABC/bin/continue_ABC.py', 'CosmoABC/bin/continue_ABC_NumCosmo.py','CosmoABC/bin/plot_ABC.py'],
      test_suite='nose.collector',
      tests_require=['nose', 'nose-cover3'],
      include_package_data=True,
      package_dir= {'CosmoABC': 'CosmoABC', 'data': 'CosmoABC/data', 'examples':'CosmoABC/examples'},
      package_data = {'CosmoABC/data':'SPT_sample.dat'},
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
