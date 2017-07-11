# Contributing to the MiCall Project #

If you like this project and want to make it better, help out. You could report
a bug, or pitch in with some development work.

## Bug Reports and Enhancement Requests ##

Please create issue descriptions [on GitHub][issues]. Be as specific as possible.
Which version are you using? What did you do? What did you expect to happen? Are
you planning to submit your own fix in a pull request?

[issues]: https://github.com/cfe-lab/MiCall/issues

## Development ##
This will document the installation steps to get MiCall running locally on your workstation.
The steps are for Eclipse with PyDev on Ubuntu, adapt as needed to your preferred IDE or operating system.

If you want to see what's currently being worked on, check out the [waffle board][waffle].

[waffle]: https://waffle.io/cfe-lab/micall

### Java and Python ###
1. Check that you are running a 64-bit operating system, or bowtie2 won't work.
   Check About this Computer under the gear menu.
2. If you want to edit Python code using PyDev and Eclipse, you will need to
   install Java. Check the version of Java you have installed:

        java -version
 
3. If the java version is lower than 1.7, then install JDK7:

        sudo apt-get install openjdk-7-source

4. Check that you are now using the new version. If not, configure it.

        java -version
        sudo update-alternatives --config java 
        java -version

5. Check the version of Python you have installed:

        python --version

6. If the Python version is lower than 2.7, then install it:
        
        sudo apt-get install python2.7
        
7. Install [pip the Python package manager][pip] and some packages for Python:

        sudo apt-get install python-pip
        sudo pip install testfixtures
        sudo pip install requests
        sudo pip install python-Levenshtein

[pip]: https://pip.pypa.io/en/latest/installing.html

### Eclipse ###
1. Install Eclipse, although you might prefer a more recent version from the [Eclipse web site][eclipse]:

        sudo apt-get install eclipse

2. Launch Eclipse. From the Help menu, choose either Eclipse Marketplace... or Install New Software....
3. In the marketplace, just type PyDev and search. In the install wizard, use the [PyDev update site][pydev].
4. After installing PyDev, open Window: Preferences. Navigate down to PyDev: Interpreters: Python Interpreter. 
5. Click the Quick Auto-Config button. Click OK.
6. From the File menu, choose Import.... Navigate down to Git: Projects from Git.
7. Choose Clone URI, and paste this URI: https://github.com/cfe-lab/MiCall.git
8. Take all the branches, and select master as your initial branch.
9. Select import existing projects, and finish the import.

[eclipse]: https://www.eclipse.org/downloads/
[pydev]: http://pydev.org/updates

### Python ###
Check that Python is already installed.

    python --version

We have tested with Python 3.4.

On Windows, you can install [Anaconda Python][anaconda].

[anaconda]: http://continuum.io/downloads

### Bowtie and HyPhy ###
1. Download the latest version of [bowtie2's binaries for Linux][bowtie2].
2. Right click and choose Extract Here. Change the folder owner to root, move it to /opt, and add it to the path.

        chmod g-w -R bowtie2-2.2.1
        sudo chown root:root -R bowtie2-2.2.1
        sudo mv bowtie2-2.2.1 /opt
        # Do this for the default version of bowtie2
        cd /opt/bowtie2-2.2.1
        for f in bowtie2* ; do sudo ln -s /opt/bowtie2-2.2.1/$f /usr/local/bin/$f; done
        # Do this for other versions that have to be called explicitly
        cd /opt/bowtie2-2.2.8
        for f in bowtie2* ; do sudo ln -s /opt/bowtie2-2.2.8/$f /usr/local/bin/$f-2.2.1; done
        # Now try a smoke test.
        cd ~
        bowtie2 --version

3. HyPhy is not needed by the main pipeline, only some of the helper utilities,
    so you can probably skip it. Before you can build HyPhy, you will need these
    libraries:

        sudo apt-get install build-essential python-dev libcurl4-openssl-dev libcrypto++-dev libssl-dev
        
    On CentOS 6, the newest versions of HyPhy fail to compile.  The newest version that works is
    v2.2.5.  In order to compile this version, assuming you're using the Software Collections `python27` package, 
    you only need to add two packages, which you can do using `yum`:
    
        sudo yum install libcurl-devel openssl-devel

4. There is a newer package for HyPhy in the
    [hyphy-python project][hyphy-python]. Consider testing that before the next
    installation, but so far we've just downloaded the latest source (or v2.2.5 on CentOS 6).
    
5. Download the latest [source for HyPhy][hyphy]. Right click the zip file and choose Expand Here. Then run the setup script:

        cd ~/Downloads/hyphy-master/src/lib
        sudo python setup.py install

    You can test it out if you like.

        cd Examples/Python
        python BasicHyPhy.py # Just check that there are no obvious errors.

[bowtie2]: http://sourceforge.net/projects/bowtie-bio/files/bowtie2/
[hyphy]: https://github.com/veg/hyphy
[hyphy-python]: https://github.com/veg/hyphy-python

### Cutadapt library ###
In order to support more than one version of this library installed in parallel,
install it in a Python virtual environment, then put a symbolic link to it on
the path.

    sudo virtualenv /usr/local/share/vcutadapt-1.11
    sudo /usr/local/share/vcutadapt-1.11/bin/pip install cutadapt==1.11
    sudo ln -s /usr/local/share/vcutadapt-1.11/bin/cutadapt /usr/local/bin/cutadapt-1.11

### Gotoh library ###
MiCall uses an implementation of a modified Gotoh algorithm for pairwise sequence alignment.
This is written in the C++ source file `gotoh.cpp`, so you will need the
Python 3 development tools.  To compile this into a shared library
that can be accessed from Python, go to `micall/alignment` and enter the following:
```
sudo python setup.py install
```
This assumes that you have superuser permissions on your system.  We have tested this
installation on OS-X and Ubuntu. If you're installing on Windows, you will need
to [install Visual C++ for Python][vcpp].

[vcpp]: http://stackoverflow.com/a/26127562/4794


### Matplotlib ###
To install it on Ubuntu, use pip:

    sudo pip install matplotlib

### BaseSpace ###
Set up the [native apps virtual machine][bsvm], and configure a shared folder
called MiCall that points to the source code. Make sure you have a developer
account on illumina.com.

Here's a little bash script to build the Docker image from the local version
of the source code, and launch it for testing:

    #!/bin/sh
    sudo docker build -t docker.illumina.com/cfelab/micall /media/sf_MiCall/ && \
    sudo docker push docker.illumina.com/cfelab/micall && \
    sudo spacedock -a <agent id from Form Builder> -m https://mission.basespace.illumina.com

That doesn't put a version number tag on the Docker image. You'll need to put
your own agent id in place, get it from the Form Builder page.

Here's a related script that pulls master from GitHub, builds the
Docker image, tags it with the version number, and launches it for testing:

    #!/bin/sh
    if [ -z $1 ]; then
      echo Missing version tag.
      exit 1
    fi
    sudo docker build -t docker.illumina.com/cfe_lab/micall:$1 https://github.com/cfe-lab/MiCall.git && \
    sudo docker push docker.illumina.com/cfe_lab/micall:$1 && \
    sudo spacedock -a <agent id from Form Builder> -m https://mission.basespace.illumina.com

[bsvm]: https://developer.basespace.illumina.com/docs/content/documentation/native-apps/setup-dev-environment

### PyInstaller ###
If you want to distribute a stand-alone Windows executable version, you need
to use [PyInstaller][pyinstaller].

1. If you haven't done all the previous steps, install Python and PyWin32.
2. If you want to rebuild the bowtie2 binaries, install it as described above.
    You might also want ActivePerl if you plan to play with bowtie2's Perl
    wrapper scripts. Micall doesn't require Perl, and you can just build with
    the bowtie2 binaries included in the git repository.
3. Install [git for Windows][wingit], and clone the MiCall repository.
4. Follow the instructions above to install the Gotoh package from the
    `micall/alignment` folder.
5. Copy `settings_default.py` to `settings.py` and edit the settings. Point
    bowtie2 at the copy in the bin folder.
6. Try running micall.py and processing the FASTQ files in
    `micall/tests/microtest`.
6. Use pip to install pyinstaller.

    pip install pyinstaller

6. Run pyinstaller.

    cd git\micall
    pyinstaller micall.spec

The application is created as `dist\micall.exe`.


[pyinstaller]: https://github.com/pyinstaller/pyinstaller
[wingit]: https://git-scm.com/download/win

### Running the code ###
1. Copy settings_default.py to settings.py, and open it for editing.
3. Change `counting_processes` to match the number of processors on your
   computer, and set `mapping_processes` to be that number divided by four.
5. Copy hostfile_default to hostfile, and open it for editing.
6. You probably just want to uncomment the localhost line.
7. Try the launch configurations. They are saved in the `micall/tests/working`
    directory, but you should see them if you open the Run menu and choose Run
    configurations.... If you want to run all steps at once, skip to the next
    step, otherwise go through the numbered launch configurations in order. If
    you are not running under Eclipse, just run each command to display the
    list of command-line parameters.
8. Copy or link all the files from the microtest folder to the working folder.
9. Run the sample_pipeline or run_processor launch configurations. They will
    process all the sample files in the working folder. If you are not running
    under Eclipse, both commands take the run folder as a command-line parameter.
12. Run the unit tests. Either run them from Eclipse, or run them from the
    command line like this:

        cd ~/git/MiCall
        python -m unittest discover -p '*_test.py'
    
### Test data ###
If you want to run MISEQ_MONITOR.py, you have to set up data folders for raw
data and for the working folders. You'll also need to set up the QAI project
and the MiseqQCReport so you can download QC data and upload results.

1. Create a data folder somewhere on your workstation, like ~/data. Create
   subdirectories called miseq and RAW_DATA. Add folders RAW_DATA/MiSeq/runs.
2. Connect to the shared drive [using CIFS][cifs] and mount 
   smb://192.168.68.144/RAW_DATA as /media/RAW_DATA.
3. Navigate down to /media/RAW_DATA/MiSeq/runs, pick a recent folder, and make
   sure it has a file named needsprocessing.
4. Copy SampleSheet.csv to a sample run folder under your local
   RAW_DATA/MiSeq/runs folder.
5. Navigate down to Data\Intensities\BaseCalls, and copy a few of the .fastq.gz
   files to your sample run folder under Data/Intensities/BaseCalls.
5. Copy the Interop folder and the files RunInfo.xml and runParameters.xml.
6. Open settings.py for editing.
7. Point `home` at your local data/miseq folder.
8. Point `rawdata_mount` at your local RAW_DATA folder.
9. Set the Oracle connection information to a test database where you can upload
   sequence data.
10. Run the Ruby console for QAI and `LabMiseqRun.import('01-Jan-2000')` for the
    date of your sample run.
11. Run the MiSeqQCReport script to upload the QC data from the sample run
    folder.
12. Run MISEQ_MONITOR.py, it doesn't take any arguments.

[cifs]: https://wiki.ubuntu.com/MountWindowsSharesPermanently

### Looking at SAM files ###
When you don't understand the pipeline's output, it can be helpful to look at
the raw reads in a sequence viewer like [Tablet][tablet]. Run the micall_basespace
script on a run with a single sample, like this:

    python micall_basespace.py --debug_remap --all_projects --link_run /path/to/run /working/path

The options tell it to write the debug files, use all projects, link to the run
with the sample you're interested in, and put all the working files in the 
given folder. Look through the scratch folders under the working path to find
the one for the sample you're interested in. The remap step writes the mapping
results as `debug_remapX_debug.sam` and `debug_remapX_debug_ref.fasta`, where
`X` is the remapping iteration number. You should be able to open an assembly
in Tablet using those two files. If the SAM file contains multiple regions,
you'll probably have to sort it with the `micall/utils/sort_sam.py` script.

[tablet]: http://ics.hutton.ac.uk/tablet/

### Releases ###
This section assumes you already have a working server up and running, and you
just want to publish a new release. If you're setting up a new server, follow
similar steps to setting up a development workstation. Follow these steps:

1. Check that all the issues in the current milestone are closed, and make sure
    the code works in your development environment. Run all the unit
    tests as described above, process the microtest data set in your local copy
    of Kive, and process all the samples from test_samples.csv using the
    `release_test_*.py` scripts to compare the results of the new release with
    the previous version. Get the comparison signed off to begin the release
    process.
2. Check if the kiveapi package needs a new release by looking for new commits.
    Make sure you tested with the latest version.
3. Determine what version number should be used next. Update the version number
    in `settings_default.py` if it hasn't been updated already, commit, and push.
4. Copy the previous pipeline on QAI/lab_miseq_pipelines to make a new version.
    Use the `projects_dump.py` script and compare `projects.json` to check that
    the projects match.
5. Check the history of the `micall.alignment` folder. If it has changed since
    the last release, then update the version number in `setup.py`.
6. [Create a release][release] on Github. Use "vX.Y" as the tag, where X.Y
    matches the version you used in `settings_default.py`. If you have to redo
    a release, you can create additional releases with tags vX.Y.1, vX.Y.2, and
    so on. Mark the release as pre-release until you finish deploying it.
7. Upgrade the scripts in Kive, and record the id of the new pipeline. You might
    find the Kive project's `dump_pipeline.py` and `upload_pipeline.py` scripts
    helpful. They are in the `utils` folder. First upgrade them on the test
    server, run them on a few samples, then upgrade them on the production
    server and run them on a few samples.
8. Stop the `MISEQ_MONITOR.py` process after you check that it's not processing
    any important runs.

        ssh user@server
        tail /data/miseq/micall.log
        ps aux|grep MISEQ_MONITOR.py
        sudo kill -int <process id from grep output>

9. Get the code from Github into the server's environment.

        ssh user@server
        cd /usr/local/share/MiCall
        git fetch
        git checkout tags/vX.Y

10. Check if you need to set any new settings by running
    `diff micall/settings_default.py micall/settings.py`. You will probably need
    to modify the version number and pipeline id, at least. Make sure that
    `production = True`.
11. Check if the gotoh package is up to date. If not, install it.

        cd /usr/local/share/MiCall/micall/alignment
        pip3 show gotoh
        cat setup.py  # compare version numbers
        sudo python3 setup.py install

12. Check that the kiveapi package is the same version you tested with. If not,
    do a Kive release first.

        cd /usr/local/share/Kive
        pip3 show kiveapi
        cat api/setup.py

13. Start the monitor, and tail the log to see that it begins processing all the
    runs with the new version of the pipeline. Before you launch, change all
    the working folders to be owned by the pipeline group.

        cd /usr/local/share/MiCall
        sudo chgrp -R micall /data/miseq
        ./MISEQ_MONITOR.py &
        tail -f /data/miseq/micall.log

14. Launch the basespace virtual machine, and build a new Docker image
    from GitHub. Tag it with the release number. See the bash scripts above for
    an easy way to do this.

    sudo docker build -t docker.illumina.com/cfe_lab/micall:vX.Y https://github.com/cfe-lab/MiCall.git

15. Push the new image to the repository. You might have to log in to docker
    before running this.
    
        sudo docker push docker.illumina.com/cfe_lab/micall:vX.Y

16. Edit the `callbacks.js` in the form builder, and add the `:vX.Y` tag to the
    `containerImageId` field.
17. Activate the new revisions in the form builder and the report builder.
18. Send an e-mail to users describing the major changes in the release.
19. Close the milestone for this release, create one for the next release, and
    decide which issues you will include in that milestone.

[release]: https://help.github.com/categories/85/articles
