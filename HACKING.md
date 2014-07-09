Setting up a developer workstation
==================================

This will document the installation steps to get the miseq pipeline running locally on your workstation.
The steps are for Eclipse with PyDev on Ubuntu, adapt as needed to your preferred IDE or operating system.

1. Check that you are running a 64-bit operating system, or bowtie2 won't work. Check About this Computer under the gear menu.
1. Check the version of Java you have installed:

        java -version
 
2. If the java version is lower than 1.7, then install JDK7:

        sudo apt-get install openjdk-7-source

3. Check that you are now using the new version. If not, configure it.

        java -version
        sudo update-alternatives --config java 
        java -version

3. Check the version of Python you have installed:

        python --version

4. If the Python version is lower than 2.7, then install it:
        
        sudo apt-get install python2.7
        
4. Install pip the Python package manager and the testfixtures package for Python:

        sudo apt-get install python-pip
        sudo pip install testfixtures

5. Install Eclipse, although you might prefer a more recent version from the [Eclipse web site][eclipse]:

        sudo apt-get install eclipse

6. Launch Eclipse. From the Help menu, choose either Eclipse Marketplace... or Install New Software....
7. In the marketplace, just type PyDev and search. In the install wizard, use the [PyDev update site][pydev].
7. After installing PyDev, open Window: Preferences. Navigate down to PyDev: Interpreters: Python Interpreter. 
7. Click the Quick Auto-Config button. Click OK.
8. From the File menu, choose Import.... Navigate down to Git: Projects from Git.
9. Choose Clone URI, and paste this URI: https://github.com/emartin-cfe/fifo_scheduler.git
10. Use the defaults, select import existing projects, and finish the import.
11. From the File menu, choose Import.... Navigate down to Git: Projects from Git.
12. Choose Clone URI, and paste this URI: https://github.com/ArtPoon/MiseqPipeline.git
14. Take all the branches, and select dev as your initial branch.
15. Select import existing projects, and finish the import.
15. Download the latest version of [bowtie2's binaries for Linux][bowtie2].
15. Right click and choose Extract Here. Change the folder owner to root, move it to /opt, and add it to the path.

        chmod g-w -R bowtie2-2.2.1
        sudo chown root:root -R bowtie2-2.2.1
        sudo mv bowtie2-2.2.1 /opt
        sudo vi /etc/environment # add :/opt/bowtie2-2.2.1 and logout/login
        bowtie2 --version # smoke test

16. Before you can build samtools, you will need these libraries:

        sudo apt-get install zlib1g-dev libncurses5-dev

16. Download the latest version of the [source for samtools][samtools].
16. Extract the files, and follow the instructions in the INSTALL file. Copy the samtools executable to /usr/bin.
16. Before you can build HyPhy, you will need these libraries:

        sudo apt-get install build-essential python-dev libcurl4-openssl-dev libcrypto++-dev libssl-dev

16. Download the latest [source for HyPhy][hyphy]. Right click the zip file and choose Expand Here. Then run the setup script:

        cd ~/Downloads/hyphy-master/src/lib
        sudo python setup.py install

    You can test it out if you like.

        cd Examples/Python
        python BasicHyPhy.py # Just check that there are no obvious errors.

16. Install Ruby for the fasta2g2p step. Check what version you have:

        ruby -v

17. If you don't have version 1.8.6, install Ruby Version Manager, and Ruby 1.8.6.
        sudo apt-get install curl
        curl -sSL https://get.rvm.io | bash -s stable
        # exit, then start a new shell so rvm will work
        sudo ls
        rvm requirements
        rvm install 1.8.6
        gem install bio
        
17. Build the alignment library.

        cd ~/git/MiseqPipeline
        ./build_alignment.sh

17. From the Help menu in Eclipse, choose Eclipse Marketplace…
17. Search for Ruby, and install Ruby (DLTK).
17. From the Window menu, select Preferences. Navigate down to Ruby: Interpreters.
17. Click Add... and browse for the Interpreter executable. Look under your
    `.rvm` folder for a path like this:
    
        ~/.rvm/rubies/ruby-1.8.6-p420/bin/ruby

17. For Interpreter arguments, type `-rubygems`.
17. Install R. The last two commands are run in the R console, and you should
    check the [StatET installation page][statet] to see exactly which version
    of the rj package is compatible with the version of StatET you are going to
    install.

        sudo apt-get install r-base r-base-dev
        sudo R
        install.packages(c("rj", "rj.gd"), repos="http://download.walware.de/rj-2.0")
        q()

    [statet]: http://www.walware.de/it/statet/installation.mframe

17. Launch Eclipse. For some reason, you can't currently install StatET from the
    Eclipse Marketplace, so from the Help menu, choose Install New Software....
17. Go to the [StatET installation page][statet], and find the update site for
    your version of Eclipse. Paste that address in the install wizard, and 
    select the StatET for R component. Finish the installation.
18. From the Window menu, choose Preferences. Navigate down to StatET: 
    Run/Debug: R Environments.
19. Click the Add... button.
20. Next to the Location (R_HOME) field, press the + button, and choose Try
    find automatically. It should find the R you just installed.
21. Click the Detect Default Properties/Settings button. Click OK. Click OK.
22. If you want an R console, open the Run menu, and choose 
    Run Configurations.... Select R Console and click the add button. Click Run.
22. To run an R script with command-line arguments, modify the R console 
    configuration by setting the working directory and adding this to the 
    Options/Arguments field with whatever CSV file name was created by the
    previous step:
    
        --args /path/to/amino_frequencies.csv /path/to/coverage_maps
    
    Then you can use `source("coverage_plots.R")` in the console to launch it.
22. Create a data folder somewhere on your workstation, like ~/data. Create subdirectories called miseq and RAW_DATA.
22. Connect to the shared drive [using CIFS][cifs] and mount smb://192.168.68.144/RAW_DATA as /media/RAW_DATA.
22. Navigate down to /media/RAW_DATA/MiSeq/runs, pick a recent folder, and make sure it has a file named needsprocessing.
22. Copy SampleSheet.csv to a sample run folder under your local RAW_DATA folder.
22. Navigate down to Data\Intensities\BaseCalls, and copy a few of the .fastq.gz files to your sample run folder.
22. Select all the .fastq.gz files you copied, right click, and choose Extract Here.
22. Delete the compressed versions of the files.
22. Copy settings_default.py to settings.py, and open it for editing.
23. Point macdatafile_mount at your local RAW_DATA folder, and set mapping_ref_path to "reference_sequences/cfe".
24. Set both mapping_factory_resources and single_thread_resources to [("", 2)]
25. If you want to reduce the combinations that run, remove all but the first 
    value in g2p_fpr_cutoffs, v3_mincounts, conseq_mixture_cutoffs. Remove all 
    but 10 from sam2csf_q_cutoffs.
25. Copy or create a symbolic link for two fastq files in the working directory.
    The files should be a matched pair: forward and reverse. Call them 
    `read1.fastq` and `read2.fastq`.
25. Try the launch configurations. They are saved in the `working` directory,
    but you should see them if you open the Run menu and choose Run
    configurations.... If you want to run all steps at once, skip down to the 
    MISEQ_PIPELINE.py file, otherwise try launching each step in the following
    order:
    * `prelim_map.py`
    * `remap.py`
    * `sam2csf.py`
    * `csf2counts.py`
    * `fasta_to_g2p.rb`
    * All together: `MISEQ_PIPELINE.py`.
... to be continued ...

[eclipse]: https://www.eclipse.org/downloads/
[pydev]: http://pydev.org/updates
[bowtie2]: http://sourceforge.net/projects/bowtie-bio/files/bowtie2/
[samtools]: http://sourceforge.net/projects/samtools/files/
[hyphy]: https://github.com/veg/hyphy
[cifs]: https://wiki.ubuntu.com/MountWindowsSharesPermanently
