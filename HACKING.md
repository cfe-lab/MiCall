Setting up a developer workstation
==================================

This will document the installation steps to get the miseq pipeline running locally on your workstation.
The steps are for Eclipse with PyDev on Ubuntu, adapt as needed to your preferred IDE or operating system.

1. Check that you are running a 64-bit operating system, or bowtie2 won't work. Check About this Computer under the gear menu.
2. Check the version of Java you have installed:

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

7. Install Eclipse, although you might prefer a more recent version from the [Eclipse web site][eclipse]:

        sudo apt-get install eclipse

8. Launch Eclipse. From the Help menu, choose either Eclipse Marketplace... or Install New Software....
9. In the marketplace, just type PyDev and search. In the install wizard, use the [PyDev update site][pydev].
10. After installing PyDev, open Window: Preferences. Navigate down to PyDev: Interpreters: Python Interpreter. 
11. Click the Quick Auto-Config button.
12. From the File menu, choose Import.... Navigate down to Git: Projects from Git.
13. Choose Clone URI, and paste this URI: https://github.com/emartin-cfe/fifo_scheduler.git
14. Use the defaults, and import existing projects.
15. Change the folder to point at the new fifo_scheduler folder created by git, and finish the import.
16. From the File menu, choose Import.... Navigate down to Git: Projects from Git.
17. Choose Clone URI, and paste this URI: https://github.com/ArtPoon/MiseqPipeline.git
18. Take all the branches, and select dev as your initial branch.
19. Select import existing projects, and finish the import.
20. Download the latest version of [bowtie2's binaries for Linux][bowtie2].
21. Right click and choose Expand Here. Change the folder owner to root, move it to /opt, and add it to the path.

        chmod g-w -R bowtie2-2.2.1
        sudo chown root:root -R bowtie2-2.2.1
        sudo mv bowtie2-2.2.1 /opt
        sudo vi /etc/environment # add :/opt/bowtie2-2.2.1 and logout/login
        bowtie2 --version # smoke test

22. Before you can build samtools, you will need these libraries:

        sudo apt-get install zlib1g-dev libncurses5-dev

23. Download the latest version of the [source for samtools][samtools].
24. Extract the files, and follow the instructions in the INSTALL file. Copy the samtools executable to /usr/bin.
25. Before you can build HyPhy, you will need these libraries:

        sudo apt-get install build-essential python-dev libcurl4-openssl-dev libcrypto++-dev libssl-dev

26. Download the latest [source for HyPhy][hyphy]. Right click the zip file and choose Expand Here. Then run the setup script:

        cd ~/Downloads/hyphy-master/src/lib
        sudo python setup.py install

27. Create a data folder somewhere on your workstation, like ~/data. Create subdirectories called miseq and RAW_DATA.
28. Connect to the shared drive [using CIFS][cifs] and mount smb://192.168.68.144/RAW_DATA as /media/RAW_DATA.
29. Navigate down to /media/RAW_DATA/MiSeq/runs, pick a recent folder, and make sure it has a file named needsprocessing.
30. Copy SampleSheet.csv to a sample run folder under your local RAW_DATA folder.
31. Navigate down to Data\Intensities\BaseCalls, and copy a few of the .fastq.gz files to your sample run folder.
32. Select all the .fastq.gz files you copied, right click, and choose Extract Here.
33. Delete the compressed versions of the files.
34. Copy settings_default.py to settings.py, and open it for editing.
35. Point `rawdata_mount` at your local RAW_DATA folder.
36. Set `base_path` to './', and comment out the next line with the production/development extensions.
37. Set both mapping_factory_resources and single_thread_resources to [("", 2)]
38. If you want to reduce the combinations that run, remove all but the first value in g2p_fpr_cutoffs, v3_mincounts, conseq_mixture_cutoffs. Remove all but 10 from sam2csf_q_cutoffs.
39. Open the MISEQ_PIPELINE.py file, and type Ctrl-F11 to try and run it. It will fail.
40. From the Run menu, choose Run Configurations....
41. Go to the Arguments tab, and type the path to the sample run folder you created above.
42. Click the Run button.
... to be continued ...

[eclipse]: https://www.eclipse.org/downloads/
[pydev]: http://pydev.org/updates
[bowtie2]: http://sourceforge.net/projects/bowtie-bio/files/bowtie2/
[samtools]: http://sourceforge.net/projects/samtools/files/
[hyphy]: https://github.com/veg/hyphy
[cifs]: https://wiki.ubuntu.com/MountWindowsSharesPermanently
