from distutils.core import setup, Extension
import sys

if sys.version_info < (3,):
    print('Sorry, MiCall-Lite requires Python version 3+.')
    print('Try Micall:\n  http://github.com/cfe-lab/MiCall')
    sys.exit()

setup(
    name='micall', 
    packages=['micall', 'micall.core', 'micall.utils', 'micall.alignment'],
    scripts=['bin/micall'],
    ext_modules=[Extension('micall.alignment._gotoh2',
                           sources=['micall/alignment/src/_gotoh2.c'])],
    package_data={'micall':
                      ['projects.json',
                       'alignment/models/*.csv'],
                  }
)
