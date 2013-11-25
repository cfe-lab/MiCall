"""
Persistent monitoring of macdatafile for runs needing processing.
Within run folders, a 'needsprocessing' flag triggers this script:
	1) Copy + unzip fastqs to the local disk
	2) Call the pipeline (0_MPI_wrapper.py)
	3) Upload results back to macdatafile
	4) Replace the 'needsprocessing' flag with 'processed'

Processing is done in serial on most recent runs first.
"""

import os, subprocess, sys, time
from datetime import datetime
from glob import glob
from miseqUtils import timestamp, sampleSheetParser
from time import sleep

if sys.version_info[:2] != (2, 7):
	timestamp("Monitor requires python 2.7")
	sys.exit()

## settings
home='/data/miseq/'		# Local path for writing data
delay = 3600			# Delay for polling macdatafile for unprocessed runs
pipeline_version = "4.1"
macdatafile_mount = '/media/macdatafile/'

def set_run_as_processed(runpath):
	os.remove(runpath)
	flag = open(runpath.replace('needsprocessing', 'processed'), 'w')
	flag.close()

# Monitor continuously checks for MiSeq folders flagged for processing
while 1:
	runs = glob(macdatafile_mount + 'MiSeq/runs/*/needsprocessing')
	if len(runs) == 0:
		timestamp('No runs need processing')
		sleep(delay)
		continue

	# Process recent runs first
	runs.sort()
	runs.reverse()
	root = runs[0].replace('needsprocessing', '')
	run_name = root.split('/')[-2]

	# Create folder locally on cluster as temp space
	if not os.path.exists(home+run_name):
		command = 'mkdir {}{}'.format(home, run_name)
		subprocess.call(command, shell=True, stdin=subprocess.PIPE, stdout = subprocess.PIPE)

	# Start a processing log
	log_file = home + run_name + '/pipeline_output.txt'
	pipeline_output = open(log_file, "wb")
	pipeline_output.write(timestamp('Processing {}'.format(root)))

	# Copy SampleSheet.csv from macdatafile to cluster
	remote_file = runs[0].replace('needsprocessing', 'SampleSheet.csv')
	local_file = home + run_name + '/SampleSheet.csv'
	command = 'rsync -a {} {}'.format(remote_file, local_file)
	pipeline_output.write(timestamp(command))
	p = subprocess.call(command, shell=True, stdin=subprocess.PIPE, stdout = subprocess.PIPE)

	try:
		infile = open(local_file, 'rU')
		run_info = sampleSheetParser(infile)
		mode = run_info['Description']
		infile.close()
	except:
		message = "!!!!!!!!!! Error parsing sample sheet - SKIPPING RUN !!!!!!!!!!"
		pipeline_output.write(timestamp(message))
		set_run_as_processed(runs[0])
		continue

	if mode not in ['Nextera', 'Amplicon']:
		message = "Error - {} is not a recognized mode: skipping'.format(mode)"
		pipeline_output.write(timestamp(message))
		set_run_as_processed(runs[0])
		continue

	# Run appears valid - copy fastq.gz files to cluster and unzip
	files = glob(root+'Data/Intensities/BaseCalls/*.fastq.gz')
	for file in files:
		filename = file.split('/')[-1]

		# Report number of reads failing to demultiplex in the log
		if filename.startswith('Undetermined'):
			output = subprocess.check_output(['wc', '-l', file])
			failed_demultiplexing = output.split(" ")[0]
			message= "{} reads failed to demultiplex in {} for run {}".format(failed_demultiplexing, filename, run_name)
			pipeline_output.write(timestamp(message))
			continue

		local_file = home + run_name + '/' + filename
		command = 'rsync -a {} {}; gunzip -f {};'.format(file, local_file, local_file)
		pipeline_output.write(timestamp(command))
		subprocess.call(command, shell=True, stdin=subprocess.PIPE, stdout = subprocess.PIPE)

	# Open file handle to the end of the current log
	read_output = open(log_file, "rb").read()
	
	# 0_MPI_wrapper will write standard output / standard error to the log (-tag-output displays MPI rank)
	command = "module load openmpi/gnu; mpirun -tag-output -machinefile mfile python -u {} {}".format("0_MPI_wrapper.py", home+run_name)
	p = subprocess.Popen(command, shell=True, stdin=subprocess.PIPE, stdout = pipeline_output)
	pipeline_output.write(timestamp(command))

	# Poll log every second for additional output and display it to the console
	while True:
		sleep(1)
		if p.poll() == 0: break
		sys.stdout.write(read_output.read())
	read_output.close()

	# Determine path to upload results (Denote the pipeline version in this path)
	result_path = runs[0].replace('needsprocessing', 'Results')
	result_path_final = result_path + '/version_' + pipeline_version
	if not os.path.exists(result_path): os.mkdir(result_path)
	if not os.path.exists(result_path_final): os.mkdir(result_path_final)

	# Determine which results files to post to macdatafile
	results_files = []
	if mode == 'Amplicon':
		results_files += glob(home + run_name + '/*.v3prot')
		results_files += glob(home + run_name + '/v3prot.summary')
	results_files += glob(home + run_name + '/*.counts')
	results_files += glob(home + run_name + '/*.csv')
	results_files += [f for f in glob(home + run_name + '/*.conseq') if 'pileup' not in f]
	timestamp("Posting {} run to macdatafile".format(mode))

	# Copy the chosen results
	for file in results_files:
		filename = file.split('/')[-1]
		command = 'rsync -a {} {}/{}'.format(file, result_path_final, filename)
		pipeline_output.write(timestamp(command))
		subprocess.call(command, shell=True, stdin=subprocess.PIPE, stdout = subprocess.PIPE)

	# Finally, copy the log file
	command = 'rsync -a {} {}/{}'.format(log_file, result_path_final, os.path.basename(log_file))
	p = subprocess.call(command, shell=True, stdin=subprocess.PIPE, stdout = subprocess.PIPE)
 	set_run_as_processed(runs[0])
	pipeline_output.write(timestamp(command))
	pipeline_output.write(timestamp("========== {} successfully completed! ==========\n".format(run_name)))
	pipeline_output.close()
