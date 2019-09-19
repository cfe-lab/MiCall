import setuptools
from distutils.core import setup, Extension
import sys

if sys.version_info < (3,):
    print('Sorry, MiCall-Lite requires Python version 3+.')
    print('Try Micall:\n  http://github.com/cfe-lab/MiCall')
    sys.exit()

with open('README.md') as fh:
    long_description = fh.read()

setup(
    name='MiCall-Lite',
    version='v0.1',
    author=['Art Poon', 'Don Kirkby', 'Eric Martin', 'Richard H. Liang'], 
    author_email='apoon42@uwo.ca',

    description='Pipeline for processing FASTQ data from an Illumina MiSeq for human RNA virus (HIV, hepatitis C virus) genotyping',
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/PoonLab/MiCall-Lite",

    packages=['micall', 'micall.core', 'micall.utils', 'micall.alignment'],
    scripts=['bin/micall'],
    ext_modules=[Extension('micall.alignment._gotoh2',
                           sources=['micall/alignment/src/_gotoh2.c'])],
    package_data={'micall':
                      ['projects.json',
                       'alignment/models/*.csv'],
                  }
)
