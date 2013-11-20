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


# run = a path to the run
def set_run_as_processed(run):
	os.remove(run)
	flag = open(run.replace('needsprocessing', 'processed'), 'w')
	flag.close()

if sys.version_info[:2] != (2, 7):
	print "Monitor requires python 2.7"
	sys.exit()


## settings
pipeline_version = "4"
home='/data/miseq/'				# Local space for writing data
macdatafile_mount = '/media/macdatafile/'
delay = 3600					# Delay between checking Miseq runs
load_mpi = "module load openmpi/gnu"
log_file_name = "pipeline_output.log"



## main loop
while 1:
	# For now, ignore standard error (Otherwise timestamp gives duplicate output)
	#sys.stderr = open(os.devnull, 'w')

	# Check if any MiSeq folders have been flagged for processing
	runs = glob(macdatafile_mount + 'MiSeq/runs/*/needsprocessing')

	if len(runs) == 0:
		timestamp('No runs need processing')
		sleep(delay)
		continue

	# Process the most recent runs first
	runs.sort()
	runs.reverse()
	root = runs[0].replace('needsprocessing', '')
	timestamp ('Processing {}'.format(root))
	run_name = root.split('/')[-2]
	if not os.path.exists(home+run_name): os.system('mkdir {}{}'.format(home, run_name))


	# Assign standard error to a log file
	log_file = open(home + run_name + "/" + log_file_name, "w")
	#sys.stderr = log_file


	# copy SampleSheet.csv from NAS to local filesystem
	remote_file = runs[0].replace('needsprocessing', 'SampleSheet.csv')
	local_file = home + run_name + '/SampleSheet.csv'
	os.system('rsync -a {} {}'.format(remote_file, local_file))

	# Extract description (mode) from SampleSheet.csv
	# Assay is ALWAYS "Nextera" and chemistry is ALWAYS "Amplicon"
	try:
		infile = open(local_file, 'rU')
		run_info = sampleSheetParser(infile)
		mode = run_info['Description']
		infile.close()
	except:
		set_run_as_processed(runs[0])
		continue

	# Mode must be Nextera or Amplicon: if not, mark as an error and proceed
	if mode not in ['Nextera', 'Amplicon']:
		timestamp('Error - \'{}\' is not a recognized mode: skipping'.format(mode))
		set_run_as_processed(runs[0])
		continue

	# If run is valid, transfer fasta.gz files to cluster and unzip
	files = glob(root+'Data/Intensities/BaseCalls/*.fastq.gz')
	for file in files:
		filename = file.split('/')[-1]

		# Skip reads that failed to demultiplex
		# FIXME: Have this save somewhere specific in a file
		if filename.startswith('Undetermined'):
			p = subprocess.Popen(['wc', '-l', path], stdout = subprocess.PIPE)
			stdout, stderr = p.communicate()
			failed_demultiplexing = (int(stdout.split()[0])/4
			timestamp("Run {} had {} reads failing to demultiplex".format(run_name, failed_demultiplexing))
			continue

		local_file = home + run_name + '/' + filename
		timestamp('rsync + gunzip {}'.format(filename))
		os.system('rsync -a {} {}'.format(file, local_file))
		os.system('gunzip -f {}'.format(local_file))
		time.sleep(1)

	command = "{}; mpirun -machinefile mfile python -u {} {}".format(load_mpi, "0_MPI_wrapper.py", home+run_name)
	timestamp(command)
	os.system(command)

	# Replace 'needsprocessing' flag with 'processed'
	set_run_as_processed(runs[0])
	result_path = runs[0].replace('needsprocessing', 'Results')
	result_path_final = result_path + '/version_' + pipeline_version

	if not os.path.exists(result_path): os.mkdir(result_path)		# Outer results folder
	if not os.path.exists(result_path_final): os.mkdir(result_path_final)	# Inner version folder

        # Post files to macdatafile
	results_files = []
	if mode == 'Amplicon':
		results_files += glob(home + run_name + '/*.v3prot')
		results_files += glob(home + run_name + '/v3prot.summary')

	results_files += glob(home + run_name + '/*.counts')
	results_files += glob(home + run_name + '/*.csv')
	results_files += [f for f in glob(home + run_name + '/*.conseq') if 'pileup' not in f]
	timestamp("Posting {} run to macdatafile".format(mode))

	for file in results_files:
		filename = file.split('/')[-1]
		command = 'rsync -a {} {}/{}'.format(file, result_path_final, filename)
		timestamp(command)
		os.system(command)

	# The run is complete - close the log, post it, re-assign sys.stderr to null
	log_file.close()

	# Merge pipeline_output.MPI_rank_*.log files (+ the monitor log) into one
	# Sort the list "logs" containing (timestamp, line) tuples
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
	with open(home + run_name + '/pipeline_output.txt', 'w') as f:
		for tuple in logs: f.write("{}\n".format(tuple[1]))

	# Decouple stderr
	#sys.stderr = open(os.devnull, 'w')
	command = 'rsync -a {} {}/{}'.format(home + run_name + '/pipeline_output.txt', result_path_final, 'pipeline_output.txt')
	timestamp(command)
	os.system(command)
	timestamp("========== {} successfully completed! ==========\n".format(run_name))
