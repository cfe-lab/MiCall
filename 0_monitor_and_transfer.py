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

import logging
from miseq_logging import init_logging

if sys.version_info[:2] != (2, 7):
	timestamp("Monitor requires python 2.7")
	sys.exit()

## Settings
home='/data/miseq/'		# Local path on cluster for writing data
delay = 3600			# Delay for polling macdatafile for unprocessed runs
pipeline_version = "4.3"
macdatafile_mount = '/media/macdatafile/'

def mark_run_as_processed(runpath):
	os.rename(runpath, runpath.replace('needsprocessing', 'processed'))

# Continuously look for MiSeq folders flagged as needing processing
while 1:
	runs = glob(macdatafile_mount + 'MiSeq/runs/*/needsprocessing')

	# FIXME DEBUG
	runs = glob(macdatafile_mount + 'MiSeq/runs/0_testing_amplicon/needsprocessing')
	# DEBUG FIXME

	if len(runs) == 0:
		timestamp('No runs need processing')
		sleep(delay)
		continue

	# Process most recently generated run and work backwards
	runs.sort()

	# FIXME DEBUG
	#runs.reverse()
	# DEBUG FIXME

	curr_run = runs[0]
	root = curr_run.replace('needsprocessing', '')
	run_name = root.split('/')[-2]

	# Make folder on the cluster for intermediary files
	if not os.path.exists(home+run_name):
		command = 'mkdir {}{}'.format(home, run_name)
		subprocess.call(command, shell=True, stdin=subprocess.PIPE, stdout = subprocess.PIPE)

	# Record all standard input / output to the pipeline log
	log_file = home + run_name + '/monitor_output.txt'
	try:
		logger = init_logging(log_file, file_log_level=logging.DEBUG, console_log_level=logging.INFO)
		logger.info('===== Processing {} with pipeline version {} ====='.format(root, pipeline_version))
	except:
		print "!!! CANNOT INITIATE LOGGING - PIPELINE HALTING RIGHT NOW !!!"
		sys.exit()

	# Copy SampleSheet.csv from macdatafile to the cluster
	remote_file = curr_run.replace('needsprocessing', 'SampleSheet.csv')
	local_file = home + run_name + '/SampleSheet.csv'
	command = 'rsync -a {} {}'.format(remote_file, local_file)
	logging.info(command)

	stdout, stderr = subprocess.Popen(command, shell=True, stdout = subprocess.PIPE, stderr = subprocess.PIPE).communicate()
	if stdout != "": logging.info(stdout)
	if stderr != "": logging.warn(stderr)

	try:
		with open(local_file, 'rU') as sample_sheet:
			run_info = sampleSheetParser(sample_sheet)
			mode = run_info['Description']
	except:
		message = "ERROR: Couldn't parse sample sheet - SKIPPING RUN"
		logging.error(message)
		mark_run_as_processed(curr_run)
		continue

	if mode not in ['Nextera', 'Amplicon']:
		message = "Error - {} is not a recognized mode: skipping'.format(mode)".format(mode)
		logging.error(message)
		mark_run_as_processed(curr_run)
		continue

	# Run appears valid - copy fastq.gz files to cluster and unzip
	files = glob(root+'Data/Intensities/BaseCalls/*.fastq.gz')
	for file in files:
		filename = file.split('/')[-1]

		# Report number of reads failing to demultiplex to the log
		if filename.startswith('Undetermined'):
			output = subprocess.check_output(['wc', '-l', file])
			failed_demultiplexing = output.split(" ")[0]
			message= "{} reads failed to demultiplex in {} for run {}".format(failed_demultiplexing, filename, run_name)
			logging.info(message)
			continue

		local_file = home + run_name + '/' + filename
		command = 'rsync -a {} {}; gunzip -f {};'.format(file, local_file, local_file)
		logging.info(command)
		stdout, stderr = subprocess.Popen(command, shell=True, stdout = subprocess.PIPE, stderr = subprocess.PIPE).communicate()
		if stdout != "": logging.info(stdout) 
		if stderr != "": logging.warn(stderr)

	# Create log where 0_MPI_wrapper will write it's output
	MPI_log_path = home + run_name + '/0_MPI_wrapper_output.log'
	MPI_wrapper_log = open(MPI_log_path, "w")

	# stdout/stderr of 0_MPI_wrapper concatenates to the log file (It includes stderr due to -tag-output)
	command = "module load openmpi/gnu; mpirun -tag-output -machinefile mfile python -u {} {}".format("0_MPI_wrapper.py", home+run_name)
	p = subprocess.Popen(command, shell=True, stdout = MPI_wrapper_log, stderr = MPI_wrapper_log)
	logging.info(command)

	# Poll the log (stdout/stderr of 0_MPI_wrapper) and display additional output to console
	with open(MPI_log_path, 'r') as f:
		MPI_wrapper_log.flush()
		sleep(1)
		sys.stdout.write(f.read())
		while True:
			sleep(1)
			if p.poll() == 0: break
			sys.stdout.write(f.read())
	MPI_wrapper_log.close()

	# Determine where to upload results (Denote version number in path)
	result_path = curr_run.replace('needsprocessing', 'Results')
	result_path_final = '{}/version_{}'.format(result_path, pipeline_version)
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

	# Copy the selected files over
	logging.info("Posting {} run to macdatafile".format(mode))
	for file in results_files:
		filename = file.split('/')[-1]
		command = 'rsync -a {} {}/{}'.format(file, result_path_final, filename)
		logging.info(command)
		stdout, stderr = subprocess.Popen(command, shell=True, stdout = subprocess.PIPE, stderr = subprocess.PIPE).communicate()
		if stdout != "": logging.info(stdout)
		if stderr != "": logging.warn(stderr)

	# Lastly, close the log and copy it over: this run is done
	command = 'rsync -a {} {}/{}'.format(log_file, result_path_final, os.path.basename(log_file))
	logging.info(command)
	logging.info("===== {} successfully completed! (Closing log) =====\n".format(run_name))
	subprocess.call(command, shell=True)
	mark_run_as_processed(curr_run)
