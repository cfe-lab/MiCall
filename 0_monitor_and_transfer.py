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

import os, sys, time
from datetime import datetime
from glob import glob
from seqUtils import timestamp
from time import sleep

if sys.version_info[:2] != (2, 7):
	print "Monitor requires python 2.7"
	sys.exit()

## settings
pipeline_version = "1.0"
home='/data/miseq/'			# Where we will write the data
delay = 3600				# Delay between checking Miseq runs
load_mpi = "module load openmpi/gnu"
qCutoff = 20

log_file_name = "pipeline_output.log"

## main loop
while 1:
	# For now, ignore standard error (Otherwise timestamp gives duplicate output)
	sys.stderr = open(os.devnull, 'w')

	# Check if any MiSeq folders have been flagged for processing
	runs = glob('/media/macdatafile/MiSeq/runs/*/needsprocessing')
	if len(runs) == 0:
		timestamp('No runs need processing')
		sleep(delay)
		continue

	# If any exist, sort list by runID (run_name) and do the first one
	runs.sort()
	root = runs[0].replace('needsprocessing', '')
	timestamp ('Processing {}'.format(root))
	run_name = root.split('/')[-2]
	if not os.path.exists(home+run_name): os.system('mkdir {}{}'.format(home, run_name))

	# Assign standard error to a log file
	log_file = open(home + run_name + "/" + log_file_name, "w")
	sys.stderr = log_file

	# Extract description (mode) from SampleSheet.csv
	# Assay is ALWAYS "Nextera" and chemistry is ALWAYS "Amplicon"
	infile = open(runs[0].replace('needsprocessing', 'SampleSheet.csv'), 'rU')
	mode = ''
	for line in infile:
		if line.startswith('Description,'):
			mode = line.strip('\n').split(',')[1]
			break
	infile.close()

	# Mode must be Nextera or Amplicon: if not, mark as an error and proceed
	if mode not in ['Nextera', 'Amplicon']:
		timestamp('Error - \'{}\' is not a recognized mode: skipping'.format(mode))
		os.remove(runs[0])
		flag = open(runs[0].replace('needsprocessing', 'processed'), 'w')
		flag.close()
		continue
	
	# If run is valid, transfer fasta.gz files to cluster and unzip
	files = glob(root+'Data/Intensities/BaseCalls/*.fastq.gz')
	for file in files:
		filename = file.split('/')[-1]

		# Skip reads that failed to demultiplex
		if filename.startswith('Undetermined'):
			continue
		
		local_file = home + run_name + '/' + filename
		timestamp('rsync + gunzip {}'.format(filename))
		os.system('rsync -a {} {}'.format(file, local_file))
		os.system('gunzip -f {}'.format(local_file))

	# Call MPI wrapper
	command = "{}; mpirun -machinefile mfile python -u {} {} {} {}".format(load_mpi, "0_MPI_wrapper.py", home+run_name, mode, qCutoff)
	timestamp(command)
	os.system(command)

	# Replace the 'needsprocessing' flag with a 'processed' flag
	os.remove(runs[0])
	flag = open(runs[0].replace('needsprocessing', 'processed'), 'w')
	flag.close()

	result_path = runs[0].replace('needsprocessing', 'Results')
	result_path_final = result_path + '/version_' + pipeline_version

	if not os.path.exists(result_path): os.mkdir(result_path)		# Outer results folder
	if not os.path.exists(result_path_final): os.mkdir(result_path_final)	# Inner version folder

	results_files = []

	if mode == 'Amplicon':
		results_files += glob(home + run_name + '/*.HIV1B-env.remap.sam')
		results_files += glob(home + run_name + '/*.v3prot')
		results_files += glob(home + run_name + '/v3prot.summary')

	elif mode == 'Nextera':
		results_files += glob(home + run_name + '/*.HXB2.sam')
		results_files += glob(home + run_name + '/HXB2.nuc_poly.summary')
		results_files += glob(home + run_name + '/HXB2.amino_poly.summary')
		results_files += glob(home + run_name + '/HXB2.nuc_poly.summary.conseq')
		results_files += glob(home + run_name + '/HXB2.amino_poly.summary.conseq')

	timestamp("Posting {} run to macdatafile".format(mode))

	for file in results_files:
		filename = file.split('/')[-1]
		command = 'rsync -a {} {}/{}'.format(file, result_path_final, filename)
		timestamp(command)
		os.system(command)


	# The run is complete - close the log, post it, re-assign sys.stderr to null
	log_file.close()

	# Merge pipeline_output.MPI_rank_*.log files (+ the monitor log) into one
	# IE, sort the list "logs" containing (timestamp, line) tuples
	log_files = glob(home + run_name + '/' + log_file_name)
	log_files += glob(home + run_name + '/pipeline_output.MPI_rank_*.log')
	logs = []
	for log_file in log_files:
		f = open(log_file, 'r')
		for line in f.readlines():
			if line.rstrip("\n") == "": continue

			try:
				fields = line.rstrip("\n").split("\t")
				myDateTime = datetime.strptime(fields[0], "%Y-%m-%d %H:%M:%S.%f")
				tuple = (myDateTime, fields[0] + "\t" + fields[1])
				logs.append(tuple)
			except:
				pass
		f.close()
		os.remove(log_file)
	logs.sort()

	# Create and post the final log file
	f = open(home + run_name + '/pipeline_output.txt', 'w')
	for tuple in logs: f.write("{}\n".format(tuple[1]))
	f.close()

	# Decouple stderr
	sys.stderr = open(os.devnull, 'w')
	command = 'rsync -a {} {}/{}'.format(home + run_name + '/pipeline_output.txt', result_path_final, 'pipeline_output.txt')
	timestamp(command)
	os.system(command)
	timestamp("========== {} successfully completed! ==========\n".format(run_name))
