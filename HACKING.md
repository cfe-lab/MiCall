Setting up a developer workstation
==================================

This will document the installation steps for Eclipse with PyDev on Windows, adapt as needed to your preferred IDE or operating system.

1. Check the version of Java you have installed:

    java -version
 
2. If the java version is lower than 1.7, then download and install JDK7 or higher from the [Oracle web site][oracle].
3. Download and install the latest version of Python 2 from the [Python web site][python].
4. Download and unzip the Eclipse standard package from the [Eclipse web site][eclipse].
5. Move the eclipse folder under your C:\Program Files folder.
6. Drag a shortcut to eclipse.exe onto your Start menu.
7. Launch Eclipse, and open Window: Preferences. Navigate down to PyDev: Interpreters: Python Interpreter. Click the Quick Auto-Config button.
8. From the File menu, choose Import.... Navigate down to Git: Projects from Git.
9. Choose Clone URI, and paste this URI: 192.168.69.159:/usr/local/git/miseqpipeline.git
10. Get the user and password from your supervisor.
11. Take all the branches, and select dev as your initial branch.
12. Run the new project wizard, choose PyDev, and change the folder to point at the new folder created by git.
13. Create a data folder somewhere on your workstation, like C:\data. Create subdirectories called miseq and RAW_DATA.
14. Start Windows Explorer by typing Windows-E, right click on your network, and choose Map network drive....
15. Choose a drive letter, like "M", and map it to \\192.168.68.144\RAW_DATA.
16. Navigate down to M:\MiSeq\runs, pick a recent folder, and make sure it has a file named needsprocessing.
17. Copy SampleSheet.csv to a sample run folder under your local RAW_DATA folder.
18. Navigate down to Data\Intensities\BaseCalls, and copy a few of the .fastq.gz files to your sample run folder.
19. Copy settings_default.py to settings.py, and open it for editing.
20. Point macdatafile_mount at your local RAW_DATA folder, replace all the bpsh strings with empty strings, and set mapping_ref_path to "reference_sequences/cfe".
21. Open the MISEQ_PIPELINE.py file, and type Ctrl-F11 to try and run it. It will fail.
... to be continued ...

[oracle]: http://www.oracle.com/technetwork/java/javase/downloads/index.html
[python]: https://www.python.org/download/
[eclipse]: https://www.eclipse.org/downloads/ 