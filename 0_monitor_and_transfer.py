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

import os, subprocess, sys, time
from datetime import datetime
from glob import glob
from miseqUtils import timestamp, sampleSheetParser
from time import sleep

if sys.version_info[:2] != (2, 7):
	print "Monitor requires python 2.7"
	sys.exit()

## settings
pipeline_version = "4"
home='/data/miseq/'				# Local path for writing data
macdatafile_mount = '/media/macdatafile/'
delay = 3600					# Delay between polling macdatafile
						# for runs that need processing

def set_run_as_processed(runpath):
	os.remove(runpath)
	flag = open(runpath.replace('needsprocessing', 'processed'), 'w')
	flag.close()

while 1:

	# Check if any MiSeq folders have been flagged for processing
	runs = glob(macdatafile_mount + 'MiSeq/runs/*/needsprocessing')

	if len(runs) == 0:
		print 'No runs need processing'
		sleep(delay)
		continue

	# Process most recent runs first
	runs.sort()
	runs.reverse()

	root = runs[0].replace('needsprocessing', '')
	timestamp('Processing {}'.format(root))
	run_name = root.split('/')[-2]

	# Start logging this run
	log_file = home + run_name + '/pipeline_output.txt'
	pipeline_output = open(log_file, "wb")

	if not os.path.exists(home+run_name):
		subprocess.call('mkdir {} {}'.format(home, run_name), shell=True, stdin=subprocess.PIPE, stdout = subprocess.PIPE)

	# Copy SampleSheet.csv from macdatafile to local
	remote_file = runs[0].replace('needsprocessing', 'SampleSheet.csv')
	local_file = home + run_name + '/SampleSheet.csv'
	command = 'rsync -a {} {}'.format(remote_file, local_file)
	message = timestamp(command)
	pipeline_output.write(message)
	p = subprocess.call(command, shell=True, stdin=subprocess.PIPE, stdout = subprocess.PIPE)

	try:
		infile = open(local_file, 'rU')
		run_info = sampleSheetParser(infile)
		mode = run_info['Description']
		infile.close()
	except:
		set_run_as_processed(runs[0])
		message = "Error parsing sample sheet - SKIPPING RUN"
		message = timestamp(message)
		pipeline_output.write("{}".format(message))
		continue

	if mode not in ['Nextera', 'Amplicon']:
		message = "Error - {} is not a recognized mode: skipping'.format(mode)"
		message = timestamp(message)
		pipeline_output.write("{}".format(message))
		set_run_as_processed(runs[0])
		continue

	# If run is valid, transfer fasta.gz files to cluster and unzip
	files = glob(root+'Data/Intensities/BaseCalls/*.fastq.gz')
	for file in files:
		filename = file.split('/')[-1]

		# Skip reads that failed to demultiplex
		if filename.startswith('Undetermined'):
			output = subprocess.check_output(['wc', '-l', file])
			failed_demultiplexing = output.split(" ")[0]
			message = timestamp("{} reads failed to demultiplex in {} for run {}".format(failed_demultiplexing, filename, run_name))
			pipeline_output.write("{}".format(message))
			continue

		local_file = home + run_name + '/' + filename
		command = 'rsync -a {} {}; gunzip -f {};'.format(file, local_file, local_file)
		message = timestamp(command)
		pipeline_output.write("{}".format(message))
		subprocess.call(command, shell=True, stdin=subprocess.PIPE, stdout = subprocess.PIPE)

	# -tag-output: Tag each line of output to stdout, stderr, and stddiag with [jobid, rank]<stdxxx>
	command = "module load openmpi/gnu; mpirun -tag-output -machinefile mfile python -u {} {}".format("0_MPI_wrapper.py", home+run_name)

	# Save the console output to pipeline_output.txt - but also display it to the console
	p = subprocess.Popen(command, shell=True, stdin=subprocess.PIPE, stdout = pipeline_output)
	message = timestamp(command)
	pipeline_output.write("{}".format(message))

	read_output = open(log_file, "rb")
	while True:
		sleep(1)
		if p.poll() == 0: break
		sys.stdout.write(read_output.read())
	read_output.close()

	# Determine path to upload results (a pipeline version-specific path)
	result_path = runs[0].replace('needsprocessing', 'Results')
	result_path_final = result_path + '/version_' + pipeline_version
	if not os.path.exists(result_path): os.mkdir(result_path)
	if not os.path.exists(result_path_final): os.mkdir(result_path_final)

	# Determine which files to post to macdatafile
	results_files = []
	if mode == 'Amplicon':
		results_files += glob(home + run_name + '/*.v3prot')
		results_files += glob(home + run_name + '/v3prot.summary')
	results_files += glob(home + run_name + '/*.counts')
	results_files += glob(home + run_name + '/*.csv')
	results_files += [f for f in glob(home + run_name + '/*.conseq') if 'pileup' not in f]
	timestamp("Posting {} run to macdatafile".format(mode))

	# Copy the results files
	for file in results_files:
		filename = file.split('/')[-1]
		command = 'rsync -a {} {}/{}'.format(file, result_path_final, filename)
		message = timestamp(command)
		pipeline_output.write("{}".format(message))
		subprocess.call(command, shell=True, stdin=subprocess.PIPE, stdout = subprocess.PIPE)

	# Copy the log file
	command = 'rsync -a {} {}/{}'.format(log_file, result_path_final, os.basename(log_file))
	message = timestamp(command)
	pipeline_output.write("{}".format(message))
	pipeline_output.close()
	p = subprocess.call(command, shell=True, stdin=subprocess.PIPE, stdout = subprocess.PIPE)
 	set_run_as_processed(runs[0])
	timestamp("========== {} successfully completed! ==========\n".format(run_name))
