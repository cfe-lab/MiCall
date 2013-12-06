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
pipeline_version = "4.4-HCV-EXPERIMENTAL"
delay = 3600					# Delay for polling macdatafile for unprocessed runs
home='/data/miseq/'				# Local path on cluster for writing data
macdatafile_mount = '/media/macdatafile/'

def post_files(files, destination):
	for file in files:
		filename = os.path.basename(file)
		command = 'rsync -a {} {}/{}'.format(file, destination, os.path.basename(file))
		logging.debug(command)
		stdout, stderr = subprocess.Popen(command, shell=True, stdout = subprocess.PIPE, stderr = subprocess.PIPE).communicate()
		if stdout != "": logging.info(stdout)
		if stderr != "": logging.warn(stderr)

# Continuously look for MiSeq folders flagged as needing processing
while 1:

	# From runs marked needsprocessing, process those unprocessed by the current pipeline version
	#runs = glob(macdatafile_mount + 'MiSeq/runs/*/needsprocessing')
	runs = glob(macdatafile_mount + 'MiSeq/runs/0_testing_amplicon/needsprocessing')

	runs_needing_processing = []
	for run in runs:
		result_path = '{}/version_{}'.format(run.replace('needsprocessing', 'Results'), pipeline_version)
		if not os.path.exists(result_path):
			runs_needing_processing.append(run)

	if len(runs_needing_processing) == 0:
		timestamp('No runs need processing')
		sleep(delay)
		continue

	# Process most recently generated run and work backwards
	runs_needing_processing.sort(reverse=True)
	curr_run = runs_needing_processing[0]
	root = curr_run.replace('needsprocessing', '')
	run_name = root.split('/')[-2]

	# Make folder on the cluster for intermediary files
	if not os.path.exists(home+run_name):
		os.mkdir(home+run_name)

	# Record all standard input / output to the pipeline log
	log_file = home + run_name + '/monitor_output.txt'
	try:
		logger = init_logging(log_file, file_log_level=logging.DEBUG, console_log_level=logging.INFO)
		logger.info('===== Processing {} with pipeline version {} ====='.format(root, pipeline_version))
	except:
		print "!!! CANNOT INITIATE LOGGING - PIPELINE HALTING NOW !!!"
		sys.exit()

	# Copy SampleSheet.csv from macdatafile to the cluster
	remote_file = curr_run.replace('needsprocessing', 'SampleSheet.csv')
	local_file = home + run_name + '/SampleSheet.csv'
	command = 'rsync -a {} {}'.format(remote_file, local_file)
	logging.debug(command)
	stdout, stderr = subprocess.Popen(command, shell=True, stdout = subprocess.PIPE, stderr = subprocess.PIPE).communicate()
	if stdout != "": logging.info(stdout)
	if stderr != "": logging.warn(stderr)

	try:
		with open(local_file, 'rU') as sample_sheet:

			# FIXME: Denote each sample as amplicon or nextera
			run_info = sampleSheetParser(sample_sheet)
			mode = run_info['Description']
	except:
		message = "ERROR: Couldn't parse sample sheet - SKIPPING RUN"
		logging.error(message)
		continue

	if mode not in ['Nextera', 'Amplicon']:
		message = "Error - {} is not a recognized mode: skipping'.format(mode)".format(mode)
		logging.error(message)
		continue

	# Run appears valid - copy fastq.gz files to cluster and unzip
	files = glob(root+'Data/Intensities/BaseCalls/*.fastq.gz')
	for file in files:
		filename = os.path.basename(file)

		# Report number of reads failing to demultiplex to the log
		if filename.startswith('Undetermined'):
			output = subprocess.check_output(['wc', '-l', file])
			failed_demultiplexing = output.split(" ")[0]
			message = "{} reads failed to demultiplex in {} for run {}".format(failed_demultiplexing, filename, run_name)
			logging.info(message)
			continue

		local_file = home + run_name + '/' + filename
		command = 'rsync -a {} {}; gunzip -f {};'.format(file, local_file, local_file)
		logging.info(command)
		stdout, stderr = subprocess.Popen(command, shell=True, stdout = subprocess.PIPE, stderr = subprocess.PIPE).communicate()
		if stdout != "": logging.info(stdout) 
		if stderr != "": logging.warn(stderr)

	# Store the output of 1_MPI_wrapper in a log, and poll it so we can continuously display output to console
	MPI_log_path = home + run_name + '/MPI_wrapper_output.log'
	with open(MPI_log_path, "w") as MPI_wrapper_log:

		# stdout/stderr of 0_MPI_wrapper concatenates to the log file (It includes stderr due to -tag-output)
		command = "module load openmpi/gnu; mpirun -tag-output -machinefile mfile python -u {} {}".format("1_MPI_wrapper.py", home+run_name)
		p = subprocess.Popen(command, shell=True, stdout = MPI_wrapper_log, stderr = MPI_wrapper_log)
		logging.info(command)

		# Poll the log (stdout/stderr of 0_MPI_wrapper) and display additional output to console
		with open(MPI_log_path, 'r') as f:
			MPI_wrapper_log.flush()
			sleep(1)
			sys.stdout.write(f.read())
			while True:
				sleep(1)
				sys.stdout.write(f.read())
				if p.poll() == 0: break

	# Determine output paths
	result_path = curr_run.replace('needsprocessing', 'Results')
	result_path_final = '{}/version_{}'.format(result_path, pipeline_version)
	log_path = '{}/logs'.format(result_path_final)
	conseq_path = '{}/consensus_sequences'.format(result_path_final)
	frequencies_path = '{}/frequencies'.format(result_path_final)

	# Create folders if needed
	if not os.path.exists(result_path): os.mkdir(result_path)
	if not os.path.exists(result_path_final): os.mkdir(result_path_final)
	if not os.path.exists(log_path): os.mkdir(log_path)
	if not os.path.exists(conseq_path): os.mkdir(conseq_path)
	if not os.path.exists(frequencies_path): os.mkdir(frequencies_path)

	# Post files to appropriate sub-folders
	logging.info("Posting results to {}".format(result_path))
	if mode == 'Amplicon':
		v3_path = '{}/v3_tropism'.format(result_path_final)
		if not os.path.exists(v3_path): os.mkdir(v3_path)
		post_files(glob(home + run_name + '/*.v3prot'), v3_path)
		post_files(glob(home + run_name + '/v3_tropism_summary.txt'), v3_path)
	post_files(glob(home + run_name + '/*.counts'), log_path)
	post_files(glob(home + run_name + '/*.log'), log_path)
	post_files([f for f in glob(home + run_name + '/*.conseq') if 'pileup' not in f], conseq_path)
	post_files(glob(home + run_name + '/*.csv'), frequencies_path)

	# Close the log and copy it over: this run is now done
	command = 'rsync -a {} {}/{}'.format(log_file, log_path, os.path.basename(log_file))
	logging.debug(command)
	logging.info("===== {} successfully completed! (Closing log) =====\n".format(run_name))
	subprocess.call(command, shell=True)
