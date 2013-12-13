"""
MISEQ_MONITOR.py
1) For runs flagged 'needsprocessing' + copy/unzip fastqs to local disk
2) Call the pipeline (MISEQ_PIPELINE.py)
3) Upload results back to macdatafile
"""

import logging, miseq_logging, miseqUtils, os, subprocess, sys, time
from glob import glob

if sys.version_info[:2] != (2, 7):
	raise Exception("Python 2.7 not detected")

## Settings
pipeline_version = "4.5b-scheduler"
delay = 3600					# Delay for polling macdatafile for unprocessed runs
home='/data/miseq/'				# Local path on cluster for writing data
macdatafile_mount = '/media/macdatafile/'

def execute_command(command):
	logger.info(" ".join(command))
	stdout, stderr = subprocess.Popen(command, stdout = subprocess.PIPE, stderr = subprocess.PIPE).communicate()
	if stderr != "": logging.warn(stderr)
	return stdout

def post_files(files, destination):
	for file in files:
		execute_command(['rsync', '-a', file, '{}/{}'.format(destination, os.path.basename(file))])

# Process runs marked as 'needsprocessing' not already processed by the current version of the pipeline
while 1:
	runs = glob(macdatafile_mount + 'MiSeq/runs/*/needsprocessing')
	runs_needing_processing = []
	for run in runs:
		result_path = '{}/version_{}'.format(run.replace('needsprocessing', 'Results'), pipeline_version)
		if not os.path.exists(result_path):
			runs_needing_processing.append(run)

	if not runs_needing_processing:
		logging.info('No runs need processing')
		time.sleep(delay)
		continue

	# Process most recently generated run and work backwards
	runs_needing_processing.sort(reverse=True)
	curr_run = runs_needing_processing[0]
	root = curr_run.replace('needsprocessing', '')
	run_name = root.split('/')[-2]

	# Make folder on the cluster for intermediary files
	if not os.path.exists(home+run_name):
		os.mkdir(home+run_name)

	# Record standard input / output of monitor
	log_file = home + run_name + '/MISEQ_MONITOR_OUTPUT.log'
	try:
		logger = miseq_logging.init_logging(log_file, file_log_level=logging.DEBUG, console_log_level=logging.INFO)
		logger.info('===== Processing {} with pipeline version {} ====='.format(root, pipeline_version))
	except:
		raise Exception("Could not create logging file - halting")

	# SampleSheet.csv will be needed to determine the mode and T-primer
	remote_file = curr_run.replace('needsprocessing', 'SampleSheet.csv')
	local_file = home + run_name + '/SampleSheet.csv'
	execute_command(['rsync', '-a', remote_file, local_file])

	try:
		with open(local_file, 'rU') as sample_sheet:
			run_info = miseqUtils.sampleSheetParser(sample_sheet)
			mode = run_info['Description']
	except:
		logging.error("Couldn't parse sample sheet: skipping run")
		continue

	if mode not in ['Nextera', 'Amplicon']:
		logging.error("{} not a valid mode: skipping run".format(mode))
		continue

	# Copy fastq.gz files to the cluster and unzip them
	for gz_file in glob(root+'Data/Intensities/BaseCalls/*.fastq.gz'):
		filename = os.path.basename(gz_file)
		local_file = home + run_name + '/' + filename
		execute_command(['rsync', '-a', gz_file, local_file])
		execute_command(['gunzip', '-f', local_file])

		# Report number of reads failing to demultiplex to the log
		if filename.startswith('Undetermined'):
			output = execute_command(['wc', '-l', local_file.replace('.gz', '')])
			failed_demultiplexing = output.split(" ")[0]
			logging.info("{} reads failed to demultiplex in {} (removing file)".format(failed_demultiplexing, filename))
			os.remove(local_file.replace('.gz', ''))
			continue

	# Store output of MISEQ_PIPELINE.py in a log + poll continuously to display output to console
	pipeline_log_path = home + run_name + '/MISEQ_PIPELINE_OUTPUT.log'
	with open(pipeline_log_path, "wb") as PIPELINE_log:

		# Standard out/error concatenates to the log
		command = ['python', '-u', 'MISEQ_PIPELINE.py', home+run_name]
		p = subprocess.Popen(command, stdout = PIPELINE_log, stderr = PIPELINE_log)
		logging.info(" ".join(command))

		# Poll the log of MISEQ_PIPELINE.py and display output to console as it appears
		with open(pipeline_log_path, 'rb') as cursor:
			while True:
				if p.poll() == 0:
					PIPELINE_log.flush()
					sys.stdout.write(cursor.read())
					break
				time.sleep(1)
				PIPELINE_log.flush()
				sys.stdout.write(cursor.read())

	# Determine output paths
	result_path = curr_run.replace('needsprocessing', 'Results')
	result_path_final = '{}/version_{}'.format(result_path, pipeline_version)
	log_path = '{}/logs'.format(result_path_final)
	conseq_path = '{}/consensus_sequences'.format(result_path_final)
	frequencies_path = '{}/frequencies'.format(result_path_final)

	# Create sub-folders if needed
	for path in [result_path, result_path_final, log_path, conseq_path, frequencies_path]:
		if not os.path.exists(path): os.mkdir(path)

	# Post files to each appropriate sub-folder
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

	# Close the log and copy it to macdatafile
	logging.info("===== {} successfully completed! =====\n".format(run_name))
	logging.shutdown()
	execute_command(['rsync', '-a', log_file, '{}/{}'.format(log_path, os.path.basename(log_file))])
