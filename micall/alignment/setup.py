from setuptools import setup, Extension
import numpy as np

setup (name = 'gotoh2',
       version = '0.1',
       description = "C implementation of Gotoh pairwise alignment algorithm to be wrapped in Python",
       py_modules = ['gotoh2'],
       ext_modules = [Extension('_gotoh2', sources = ['src/_gotoh2.c'])],
       data_files = [('models/', ['models/HYPHY_NUC.csv',
                                  'models/NWALIGN.csv',
                                  'models/Biopp.csv',
                                  'models/EmpHIV25.csv'])],
       include_dirs = [np.get_include()],
       zip_safe = False  #avoid headache with permissions on ~/.python-eggs
)
