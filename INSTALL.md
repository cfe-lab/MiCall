# Installing MiCall-Lite

## Requirements
* [Python 3.x](https://www.python.org/)
* [bowtie2](https://github.com/BenLangmead/bowtie2/releases/tag/v2.2.8), version 2.2.8
  MiCall uses strict version control on calling external programs.  In other words, it checks whether the program that it needs to run has the same version number as what it expects.  If the version number does not match, then MiCall will fail the run.  To change this version number, you need to edit the line `BOWTIE_VERSION` at the top of the file `micall/core/prelim_map.py`.
* the Python [Levenshtein] module

## Instructions for Ubuntu (16.04) Linux
In progress.
