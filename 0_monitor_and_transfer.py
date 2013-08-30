"""
Persistent monitoring of MiSeq runs folder mounted at /media/macdatafile
for 'needsprocessing' flag, which will trigger a pipeline on the cluster.

The role of this script is to detect this flag, transfer the *.fastq.gz
files in Data/Intensities/BaseCalls/, uncompress the files, and call
the pipeline.  When the pipeline is complete, it will upload the end
results back to macdatafile, create an empty file named 'processingcomplete'
and remove 'needsprocessing'.

Processing will be done one a "first come first serve" basis, i.e., no
asynchronous processing where another MiSeq run becomes available for 
processing while the pipeline is being executed.
"""

import os
import sys
from glob import glob
from datetime import datetime
from time import sleep


#home = '/usr/local/share/miseq/data/'
home='/data/miseq/'
delay = 3600
"""
try:
	np = int(sys.argv[1])
except:
	print '\nUsage:\npython 0_monitor_and_transfer.py [number of MPI processes]\n'
	raise
"""

def timestamp(msg):
	# Display the date/time along with a message
	print '[%s] %s' % (datetime.now().strftime('%Y-%m-%d %H:%M:%S'), msg)


while 1:

	# Check for any runs with a needsprocessing flag
	runs = glob('/media/macdatafile/MiSeq/runs/*/needsprocessing')
	if len(runs) == 0:
		timestamp('no runs need processing')
		sleep(delay)
		continue
	
	# If any exist, sort list by runID (run_name)
	runs.sort()
	root = runs[0].replace('needsprocessing', '')
	timestamp ('processing %s' % root)
	run_name = root.split('/')[-2]
	if not os.path.exists(home+run_name):
		os.system('mkdir %s%s' % (home, run_name))
	
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
	
	# mode must be Nextera or Amplicon: if not, mark as processed and proceed
	if mode not in ['Nextera', 'Amplicon']:
		print 'ERROR: Unrecognized mode "', mode, '"... skipping'

		# Eliminate 'needsprocessing' flag, add 'processed' flag
		os.remove(runs[0])
		flag = open(runs[0].replace('needsprocessing', 'processed'), 'w')
		flag.close()

		# Start over: this run will no longer be grepped
		continue
	
	
	# If run is valid, transfer fasta.gz files to cluster and unzip
	files = glob(root+'Data/Intensities/BaseCalls/*.fastq.gz')
	nfiles = 0

	for file in files:
		filename = file.split('/')[-1]
		if filename.startswith('Undetermined'):
			# Do not process undetermined read files
			continue
		
		local_file = home + run_name + '/' + filename
		timestamp('cp and gunzip %s' % filename)
		os.system('cp %s %s' % (file, local_file))
		os.system('gunzip -f %s' % local_file)
		nfiles += 1

	# Why do we divide nFiles by 2?
	# paired-end reads, each sample has two FASTQ files
	nfiles /= 2
		
	"""
	# prepare a machine file 
	outfile = open('mfile', 'w')
	for node in ['n0', 'Bulbasaur']:
		outfile.write('%s slots=%d max_slots=24\n' % (node, np))
	outfile.close()
	"""
	
	timestamp('spawning MPI processes...')

	# module load openmpi/gnu: sets up the correct alias for mpirun

	# mpirun: reads mfile and loads k processes of 0_MPI_wrapper.py
	#         (The mfile specifies the cluster, and the # CPUs to use)

	os.system('module load openmpi/gnu; mpirun -machinefile mfile python /usr/local/share/miseq/scripts/0_MPI_wrapper.py %s %s' % (home+run_name, mode))

	#os.system('module load openmpi/gnu; mpirun -machinefile mfile python /usr/local/share/miseq/scripts/0_MPI_wrapper_TEST.py %s %s' % (home+run_name, mode))
	

	# Replace the 'needsprocesing' flag with a 'processed' flag
	os.remove(runs[0])
	flag = open(runs[0].replace('needsprocessing', 'processed'), 'w')
	flag.close()

	# move results to macdatafile
	result_path = runs[0].replace('needsprocessing', 'Results')
	if not os.path.exists(result_path):
		os.mkdir(result_path)
	
	if mode == 'Amplicon':
		files = glob(home + run_name + '/*.fasta')
		files += glob(home + run_name + '/_g2p.csv')
	else:
		files = glob(home + run_name + '/*.csf')
		
	files += glob(home + run_name + '/*.amino.csv')
	files += glob(home + run_name + '/*.v3prot')
	
	for file in files:
		filename = file.split('/')[-1]
		os.system('cp %s %s/%s' % (file, result_path, filename))

