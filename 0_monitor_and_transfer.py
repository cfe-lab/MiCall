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
from glob import glob
from datetime import datetime


home = '/usr/local/share/miseq/data/'
delay = 3600

def timestamp(msg):
	print '[%s] %s' % (datetime.now().strftime('%Y-%m-%d %H:%M:%S'), msg)


while 1:
	runs = glob('/media/macdatafile/MiSeq/runs/*/needsprocessing')
	if len(runs) == 0:
		timestamp('no runs need processing')
		sleep(delay)
		continue
	
	# sort with respect to run 
	runs.sort()
	root = runs[0].replace('needsprocessing', '')
	timestamp ('processing %s' % root)
	run_name = root.split('/')[-2]
	if not os.path.exists(home+run_name):
		os.system('mkdir %s%s' % (home, run_name))
	
	# transfer *.fasta.gz files
	files = glob(root+'Data/Intensities/BaseCalls/*.fastq.gz')
	nfiles = 0
	for file in files:
		filename = file.split('/')[-1]
		if filename.startswith('Undetermined'):
			# skip undetermined read files
			continue
		
		local_file = home + run_name + '/' + filename
		timestamp('cp and gunzip %s' % filename)
		os.system('cp %s %s' % (file, local_file))
		os.system('gunzip -f %s' % local_file)
		nfiles += 1
	
	# prepare a machine file 
	outfile = open('mfile', 'w')
	if nfiles <= 6:
		outfile.write('n0 slots=%d max_slots=24\n' % (nfiles, ))
	elif nfiles <= 12:
		outfile.write('n0 slots=6 max_slots=24\n')
		outfile.write('Bulbasaur slots=%d max_slots=24\n' % ((nfiles-6), ))
	else:
		outfile.write('n0 slots=6 max_slots=24\n')
		outfile.write('Bulbasaur slots=6 max_slots=24\n')
		
	outfile.close()
	
	timestamp('spawning MPI processes...')
	os.system('/opt/scyld/openmpi/1.6.4/gnu/bin/mpirun -machinefile mfile python 0_MPI_wrapper.py %s' % home+run_name)
	
	# at this point, erase the needsprocessing file
	break




