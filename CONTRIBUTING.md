Setting up a developer workstation
==================================

This will document the installation steps to get the miseq pipeline running locally on your workstation.
The steps are for Eclipse with PyDev on Ubuntu, adapt as needed to your preferred IDE or operating system.

Java, Python, and Oracle
------------------------
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
        
7. Install pip the Python package manager and the testfixtures package for Python:

        sudo apt-get install python-pip
        sudo pip install testfixtures

8.  Install [Oracle Instant Client][oracle]. Use the basic lite version, and 
    test that sqlplus works. You will probably have to follow the steps to set 
    up the libraries, and you may have to run sqlplus64 instead of sqlplus.

        sqlplus USER@\"//192.168.67.9:1521/CFE9ir2\"

    If you want to have history and tab expansion in sqlplus, install rlwrap:

        sudo apt-get install rlwrap
        alias sqlplus="rlwrap sqlplus"

    You also need to set the `ORACLE_HOME` environment variable and make it
    visible to sudo.

        sudo vi /etc/profile.d/oracle.sh # Add the next line:
        export ORACLE_HOME=/usr/lib/oracle/12.1/client64

    [oracle]: https://help.ubuntu.com/community/Oracle%20Instant%20Client

9. To use Oracle from Python, you will need the cx_Oracle package.

        sudo env ORACLE_HOME=$ORACLE_HOME pip install cx_Oracle

Eclipse
-------
1. Install Eclipse, although you might prefer a more recent version from the [Eclipse web site][eclipse]:

        sudo apt-get install eclipse

2. Launch Eclipse. From the Help menu, choose either Eclipse Marketplace... or Install New Software....
3. In the marketplace, just type PyDev and search. In the install wizard, use the [PyDev update site][pydev].
4. After installing PyDev, open Window: Preferences. Navigate down to PyDev: Interpreters: Python Interpreter. 
5. Click the Quick Auto-Config button. Click OK.
6. From the File menu, choose Import.... Navigate down to Git: Projects from Git.
7. Choose Clone URI, and paste this URI: https://github.com/ArtPoon/MiseqPipeline.git
8. Take all the branches, and select dev as your initial branch.
9. Select import existing projects, and finish the import.

Bowtie, Samtools, and Hyphy
---------------------------
1. Download the latest version of [bowtie2's binaries for Linux][bowtie2].
2. Right click and choose Extract Here. Change the folder owner to root, move it to /opt, and add it to the path.

        chmod g-w -R bowtie2-2.2.1
        sudo chown root:root -R bowtie2-2.2.1
        sudo mv bowtie2-2.2.1 /opt
        sudo vi /etc/environment # add :/opt/bowtie2-2.2.1 and logout/login
        bowtie2 --version # smoke test

3. Before you can build samtools, you will need these libraries:

        sudo apt-get install zlib1g-dev libncurses5-dev

4. Download the latest version of the [source for samtools][samtools].
5. Extract the files, and follow the instructions in the INSTALL file. Copy the samtools executable to /usr/bin.
6. Before you can build HyPhy, you will need these libraries:

        sudo apt-get install build-essential python-dev libcurl4-openssl-dev libcrypto++-dev libssl-dev

7. Download the latest [source for HyPhy][hyphy]. Right click the zip file and choose Expand Here. Then run the setup script:

        cd ~/Downloads/hyphy-master/src/lib
        sudo python setup.py install

    You can test it out if you like.

        cd Examples/Python
        python BasicHyPhy.py # Just check that there are no obvious errors.

Ruby
----
1. Install Ruby for the fasta2g2p step. Check what version you have:

        ruby -v

2. If you don't have version 1.8.6, install Ruby Version Manager, and Ruby 1.8.6.
        sudo apt-get install curl
        curl -sSL https://get.rvm.io | bash -s stable
        # exit, then start a new shell so rvm will work
        sudo ls
        rvm requirements
        rvm install 1.8.6
        gem install bio
        gem install fastercsv

3. It sounds like executable-hooks is installed with Ruby 1.8.6, even though
    it is incompatible. It causes the warning: `parenthesize argument(s) for 
    future version`. Follow the removal instructions from
    [the bug report][rvm-bug]. The first command complains a lot, but it seems
    to work.

        executable-hooks-uninstaller
        rvm @global do gem uninstall -ax rubygems-bundler executable-hooks bundler-unload

    [rvm-bug]: https://github.com/wayneeseguin/rvm/issues/2325        
4. Build the alignment library.

        cd ~/git/MiseqPipeline
        ./build_alignment.sh

5. From the Help menu in Eclipse, choose Eclipse Marketplace…
6. Search for Ruby, and install Ruby (DLTK).
7. From the Window menu, select Preferences. Navigate down to Ruby: Interpreters.
8. Click Add... and browse for the Interpreter executable. Look under your
    `.rvm` folder for a path like this:
    
        ~/.rvm/rubies/ruby-1.8.6-p420/bin/ruby

9. For Interpreter arguments, type `-rubygems`.

R
-
1. Install R. The last two commands are run in the R console, and you should
    check the [StatET installation page][statet] to see exactly which version
    of the rj package is compatible with the version of StatET you are going to
    install.

        sudo apt-get install r-base r-base-dev
        sudo R
        install.packages(c("rj", "rj.gd"), repos="http://download.walware.de/rj-2.0")
        q()

    [statet]: http://www.walware.de/it/statet/installation.mframe

2. Launch Eclipse. For some reason, you can't currently install StatET from the
    Eclipse Marketplace, so from the Help menu, choose Install New Software....
3. Go to the [StatET installation page][statet], and find the update site for
    your version of Eclipse. Paste that address in the install wizard, and 
    select the StatET for R component. Finish the installation.
4. From the Window menu, choose Preferences. Navigate down to StatET: 
    Run/Debug: R Environments.
5. Click the Add... button.
6. Next to the Location (R_HOME) field, press the + button, and choose Try
    find automatically. It should find the R you just installed.
7. Click the Detect Default Properties/Settings button. Click OK. Click OK.
8. From the Window menu, choose Preferences. Navigate down to 
    StatET: R Code Formatting. Change the policy to use spaces. Click OK.
9. If you want an R console, open the Run menu, and choose 
    Run Configurations.... Select R Console and click the add button. Click Run.
10. To run an R script with command-line arguments, modify the R console 
    configuration by setting the working directory and adding this to the 
    Options/Arguments field with whatever CSV file name was created by the
    previous step:
    
        --args /path/to/amino_frequencies.csv /path/to/coverage_maps
    
    Then you can use `source("coverage_plots.R")` in the console to launch it.

MPI
---
1. Install Open MPI.

    sudo apt-get install openmpi-bin
    
2. Install mpi4py.

    sudo apt-get install python-mpi4py

Running the code
----------------
1. Copy settings_default.py to settings.py, and open it for editing.
2. Set `base_path` to '../', and comment out the next line with the
   production/development extensions.
3. Change `counting_processes` to match the number of processors on your
   computer, and set `mapping_processes` to be that number divided by four.
4. If you want to reduce the combinations that run, remove all but the first 
    value in g2p_fpr_cutoffs, v3_mincounts, conseq_mixture_cutoffs. Remove all 
    but 15 from sam2csf_q_cutoffs.
5. Copy hostfile_default to hostfile, and open it for editing.
6. You probably just want to uncomment the localhost line.
7. Try the launch configurations. They are saved in the `working` directory,
    but you should see them if you open the Run menu and choose Run
    configurations.... If you want to run all steps at once, skip to the next
    step, otherwise go through the numbered launch configurations in order.
8. Copy all the files from the microtest folder to the working folder.
9. Run the sample_pipeline or run_processor launch configurations. They will
    process all the sample files in the working folder. 
12. Run the unit tests. Either run them from Eclipse, or run them from the
    command line like this:

        cd ~/git/MiseqPipeline
        python -m unittest discover -p '*_test.py'
    
Test data
---------
If you want to run MISEQ_MONITOR.py, you have to set up data folders for raw
data and for the working folders.

1. Create a data folder somewhere on your workstation, like ~/data. Create
   subdirectories called miseq and RAW_DATA. Add folders RAW_DATA/MiSeq/runs.
2. Connect to the shared drive [using CIFS][cifs] and mount 
   smb://192.168.68.144/RAW_DATA as /media/RAW_DATA.
3. Navigate down to /media/RAW_DATA/MiSeq/runs, pick a recent folder, and make
   sure it has a file named needsprocessing.
4. Copy SampleSheet.csv to a sample run folder under your local
   RAW_DATA/MiSeq/runs folder.
5. Navigate down to Data\Intensities\BaseCalls, and copy a few of the .fastq.gz
   files to your sample run folder.
6. Open settings.py for editing.
7. Point `home` at your local data/miseq folder.
8. Point `rawdata_mount` at your local RAW_DATA folder.
9. Set the Oracle connection information to a test database where you can upload
   sequence data.
10. Run MISEQ_MONITOR.py, it doesn't take any arguments.

Releases
--------
This section assumes you already have a working server up and running, and you
just want to publish a new release. If you're setting up a new server, follow
similar steps to setting up a development workstation. Follow these steps:
1. Make sure the code works in your development environment. Run all the unit
    tests as described above, process the microtest data set, and process a full
    run using MISEQ_MONITOR.py. Check the logs for errors. Also check that all
    the issues in the current milestone are closed.
2. Determine what version number should be used next. Update the version number
    in `settings_default.py` if it hasn't been updated already, commit, and push.
3. [Create a release][release] on Github. Use "vX.Y" as the tag, where X.Y
    matches the version you used in `settings_default.py`. If you have to redo
    a release, you can create additional releases with tags vX.Y.1, vX.Y.2, and
    so on.
4. Get the code from Github into the server's development environment.
        ssh user@server
        cd /usr/local/share/miseq/development/
        git fetch github
        git checkout tags/vX.Y
5. Check if you need to set any new settings by running
    `diff settings_default.py settings.py`. You will probably need to modify
    the version number, at least. Make sure that `production = False`, and the
    process counts are half the production values. Do the same comparison of
    `hostfile`.
6. Process one full run of data.
        cd /usr/local/share/miseq/development/
        ./run_processor.py /data/miseq/YYMMDD*
7. Stop the `MISEQ_MONITOR.py` process after you check that it's not processing
    any runs.
        ssh user@server
        tail /data/miseq/MISEQ_MONITOR_OUTPUT.log
        ps aux|grep MISEQ_MONITOR.py
        sudo kill -9 <process id from grep output>
8. Get the code from Github into the server's production environment.
        ssh user@server
        cd /usr/local/share/miseq/production/
        git fetch
        git checkout tags/vX.Y
9. Review the settings and host file just as you did in the development
    environment, but make sure that `production = True`.
10. Start the monitor, and tail the log to see that it begins processing all the
    runs with the new version of the pipeline.
        cd /usr/local/share/miseq/production/
        python MISEQ_MONITOR.py &>/dev/null &
        tail -f /data/miseq/MISEQ_MONITOR_OUTPUT.log
11. Close the milestone for this release, create one for the next release, and
    decide which issues you will include in that milestone.

[eclipse]: https://www.eclipse.org/downloads/
[pydev]: http://pydev.org/updates
[bowtie2]: http://sourceforge.net/projects/bowtie-bio/files/bowtie2/
[samtools]: http://sourceforge.net/projects/samtools/files/
[hyphy]: https://github.com/veg/hyphy
[cifs]: https://wiki.ubuntu.com/MountWindowsSharesPermanently
[release]: https://help.github.com/categories/85/articles