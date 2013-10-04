"""
Persistent monitoring of macdatafile for runs that need processing

Within run folders, a 'needsprocessing' flag triggers this script to:
	1) Copy + unzip the file from macdatafile to a local disk
	2) Call the pipeline (0_MPI_wrapper.py)
	3) Upload results back to macdatafile
	4) Replace the 'needsprocessing' flag with 'processingcomplete'

Processing will be done in serial on a 'first come first serve' basis:
no asynchronous processing of multiple concurrent runs
"""

import os, sys
from datetime import datetime
from glob import glob
from seqUtils import timestamp
from time import sleep

home='/data/miseq/'
delay = 3600

while 1:
	runs = glob('/media/macdatafile/MiSeq/runs/*/needsprocessing')
	if len(runs) == 0:
		timestamp('No runs need processing')
		sleep(delay)
		continue
	
	# If any exist, sort list by runID (run_name)
	runs.sort()
	root = runs[0].replace('needsprocessing', '')
	timestamp ('Processing {}'.format(root))
	run_name = root.split('/')[-2]
	if not os.path.exists(home+run_name): os.system('mkdir {}{}'.format(home, run_name))
	
	# Extract assay/description/chemistry from SampleSheet.csv
	infile = open(runs[0].replace('needsprocessing', 'SampleSheet.csv'), 'rU')
	assay, mode, chemistry = '', '', ''
	for line in infile:
		if line.startswith('Assay,'):
			assay = line.strip('\n').split(',')[1]
		elif line.startswith('Description,'):
			mode = line.strip('\n').split(',')[1]
		elif line.startswith('Chemistry,'):
			chemistry = line.strip('\n').split(',')[1]
			break
	infile.close()
	
	# Mode must be Nextera or Amplicon: if not, mark as an error and proceed
	if mode not in ['Nextera', 'Amplicon']:
		timestamp('Error - \'{}\' is not a recognized mode: skipping'.format(mode))
		os.remove(runs[0])
		flag = open(runs[0].replace('needsprocessing', 'processed'), 'w')
		flag.close()
		continue
	
	# If run is valid, transfer fasta.gz files to cluster and unzip (Except undetermined reads)
	files = glob(root+'Data/Intensities/BaseCalls/*.fastq.gz')
	for file in files:
		filename = file.split('/')[-1]
		if filename.startswith('Undetermined'): continue
		local_file = home + run_name + '/' + filename
		timestamp('cp and gunzip {}'.format(filename))
		os.system('cp {} {}'.format(file, local_file))
		os.system('gunzip -f {}'.format(local_file))

	load_mpi = "module load openmpi/gnu"
	script_path = "/usr/local/share/miseq/development/miseqpipeline/0_MPI_wrapper.py"
	qCutoff = 20
	command = "{}; mpirun -machinefile mfile python -u {} {} {} {}".format(load_mpi, script_path, home+run_name, mode, qCutoff)
	timestamp(command)
	os.system(command)

	# Replace the 'needsprocesing' flag with a 'processed' flag
	os.remove(runs[0])
	flag = open(runs[0].replace('needsprocessing', 'processed'), 'w')
	flag.close()
	result_path = runs[0].replace('needsprocessing', 'Results')

	# Post files to macdatafile
	if not os.path.exists(result_path): os.mkdir(result_path)
	results_files = []
	results_files += glob(home + run_name + '/*.HXB2.sam')
	results_files += glob(home + run_name + '/HXB2.nuc_poly.summary')
	results_files += glob(home + run_name + '/HXB2.amino_poly.summary')

	for file in results_files:
		filename = file.split('/')[-1]
		command = 'cp {} {}/{}'.format(file, result_path, filename)
		timestamp(command)
		os.system(command)
