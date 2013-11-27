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

## Settings
home='/data/miseq/'		# Local path on cluster for writing data
delay = 3600			# Delay for polling macdatafile for unprocessed runs
pipeline_version = "4.1"
macdatafile_mount = '/media/macdatafile/'

def mark_run_as_processed(runpath):
	os.rename(runpath, runpath.replace('needsprocessing', 'processed'))

# Continuously look for MiSeq folders flagged as needing processing
while 1:
	runs = glob(macdatafile_mount + 'MiSeq/runs/*/needsprocessing')
	if len(runs) == 0:
		timestamp('No runs need processing')
		sleep(delay)
		continue

	# Process most recently generated run and work backwards
	runs.sort()
	runs.reverse()
	curr_run = runs[0]
	root = curr_run.replace('needsprocessing', '')
	run_name = root.split('/')[-2]

	# Make folder on the cluster for intermediary files
	if not os.path.exists(home+run_name):
		command = 'mkdir {}{}'.format(home, run_name)
		subprocess.call(command, shell=True, stdin=subprocess.PIPE, stdout = subprocess.PIPE)

	# Record all standard input / output to the pipeline log
	log_file = home + run_name + '/pipeline_output.txt'
	pipeline_output = open(log_file, "wb")
	pipeline_output.write(timestamp('===== Processing {} with pipeline version {} ====='.format(root, pipeline_version)))

	# Copy SampleSheet.csv from macdatafile to the cluster
	remote_file = curr_run.replace('needsprocessing', 'SampleSheet.csv')
	local_file = home + run_name + '/SampleSheet.csv'
	command = 'rsync -a {} {}'.format(remote_file, local_file)
	pipeline_output.write(timestamp(command))
	subprocess.call(command, shell=True, stdin=subprocess.PIPE, stdout = subprocess.PIPE)

	try:
		infile = open(local_file, 'rU')
		run_info = sampleSheetParser(infile)
		mode = run_info['Description']
		infile.close()
	except:
		message = "ERROR: Couldn't parse sample sheet - SKIPPING RUN"
		pipeline_output.write(timestamp(message))
		mark_run_as_processed(curr_run)
		continue

	if mode not in ['Nextera', 'Amplicon']:
		message = "Error - {} is not a recognized mode: skipping'.format(mode)".format(mode)
		pipeline_output.write(timestamp(message))
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
			pipeline_output.write(timestamp(message))
			continue

		local_file = home + run_name + '/' + filename
		command = 'rsync -a {} {}; gunzip -f {};'.format(file, local_file, local_file)
		pipeline_output.write(timestamp(command))
		subprocess.call(command, shell=True, stdin=subprocess.PIPE, stdout = subprocess.PIPE)

	# Direct a file handle to the end of the current log
	read_output = open(log_file, "rb")
	read_output.read()
	
	# stdout/stderr of 0_MPI_wrapper concatenates to the end of the log (It includes stderr due to -tag-output)
	command = "module load openmpi/gnu; mpirun -tag-output -machinefile mfile python -u {} {}".format("0_MPI_wrapper.py", home+run_name)
	p = subprocess.Popen(command, shell=True, stdin=subprocess.PIPE, stdout = pipeline_output)
	pipeline_output.write(timestamp(command))

	# Continuously poll the log (stdout/stderr of 0_MPI_wrapper) and display any new output to console
	while True:
		sleep(1)
		if p.poll() == 0: break
		sys.stdout.write(read_output.read())
	read_output.close()

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
	pipeline_output.write(timestamp("Posting {} run to macdatafile".format(mode)))
	for file in results_files:
		filename = file.split('/')[-1]
		command = 'rsync -a {} {}/{}'.format(file, result_path_final, filename)
		pipeline_output.write(timestamp(command))
		subprocess.call(command, shell=True, stdin=subprocess.PIPE, stdout = subprocess.PIPE)

	# Lastly, close the log and copy it over: this run is done
	command = 'rsync -a {} {}/{}'.format(log_file, result_path_final, os.path.basename(log_file))
	pipeline_output.write(timestamp(command))
	pipeline_output.write(timestamp("===== {} successfully completed! (Closing log) =====\n".format(run_name)))
	pipeline_output.close()
	subprocess.call(command, shell=True, stdin=subprocess.PIPE, stdout = subprocess.PIPE)
 	mark_run_as_processed(curr_run)
