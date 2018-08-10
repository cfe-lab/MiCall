# Installing MiCall-Lite

## Requirements
* [Python 3.x](https://www.python.org/)
* [bowtie2](https://github.com/BenLangmead/bowtie2/releases/tag/v2.2.8), version 2.2.8

  MiCall uses strict version control on calling external programs.  In other words, it checks whether the program that it needs to run has the same version number as what it expects.  If the version number does not match, then MiCall will fail the run.  To change this version number, you need to edit the line `BOWTIE_VERSION` at the top of the file `micall/core/prelim_map.py`.
* the Python [Levenshtein](https://pypi.org/project/python-Levenshtein/) module

## Instructions for Linux (Ubuntu)
These instructions assume that you have a basic working knowledge of Linux.  If you are new to this type of operating system, you should familiarize yourself with how to navigate the file system, utilize package managers, and launch programs.  See [linuxcommand.org](http://linuxcommand.org) for a helpful free and online resource.  These particular instructions were written for Ubuntu version 16.04, which is a popular distribution of Linux.  If you are running a RHEL-based distribution, then you should be able to use the same procedure, except that you would use the `yum` package manager instead of `apt`.

These instructions also assume that you have superuser access.  If you do not, then you will either need to ask a system administrator to install MiCall-Lite for you, or adapt these instructions to perform a local installation (more difficult).

1. Obtain Python 3.x.  Many Linux distros ship with Python 2.x by default and you may need to use a package manager to install this newer major version:
   ```
   sudo apt install python3.5-dev
   ```
   
2. Install the [Levenshtein](https://en.wikipedia.org/wiki/Levenshtein_distance) module:
   ```
   sudo apt install python3-levenshtein
   ```

3. Compile and install bowtie2:
   ```
   wget https://github.com/BenLangmead/bowtie2/archive/v2.2.8.tar.gz
   gunzip v2.2.8.tar.gz
   tar -xf bowtie2-2.2.8.tar
   cd bowtie2-2.2.8
   make
   sudo cp bowtie2* /usr/local/bin/
   ```
   You might prefer keeping the executables in a directory other than `/usr/local/bin`, but you should make sure that your preferred directory is in your `PATH`.
   
4. Download this repository, either as a [release](https://github.com/PoonLab/MiCall-Lite/releases) or a clone.  To clone the repository, move to a directory in your filesystem where you want to store the source code and resources (I use `~/git`) and run the following command:
   ```
   git clone https://github.com/PoonLab/MiCall-Lite.git
   ```

5. Compile and install MiCall-Lite:
   ```
   cd MiCall-Lite
   sudo python3 setup.py install
   ```

