"""
Manage the parallel execution of the first step of the pipeline.
"""

import logging, miseq_logging, miseq_modules, miseqUtils, os, subprocess, sys, time
sys.path.append('/usr/local/share/fifo_scheduler')
from fifo_scheduler import Factory
from glob import glob

## Arguments
#root = sys.argv[1]
root = "/data/miseq/130711_M01841_0010_000000000-A3TCY"

## Reference sequences
mapping_ref_path = "/usr/local/share/miseq/refs/cfe"
consensus_q_cutoff = 20
mode = "Amplicon"
REMAP_THRESHOLD = 0.95
MAX_REMAPS = 3 
bowtie_threads = 8

# Assign two workers to node -1, two workers to node 0, etc
assigned_resources = [("bpsh -1", 1), ("bpsh 0", 1), ("bpsh 1", 1), ("bpsh 2", 1)]
mapping_factory = Factory(assigned_resources)
for rule in assigned_resources:
	for i in range(rule[1]):
		mapping_factory.hire_worker(rule[0])

fastq_files = glob(root + '/*R1*.fastq')
fastq_files = [f for f in fastq_files if not f.endswith('.Tcontaminants.fastq')]

for i, fastq in enumerate(fastq_files):
	fastq_filename = os.path.basename(fastq)
	sample_name = fastq_filename.split('_')[0]
	is_t_primer = False				# Manual hack for proof of concept

	command = "python STEP_ONE_MODULE_CALL.py {} {} {} {} {} {} {} {}".format(mapping_ref_path,
			fastq, consensus_q_cutoff, mode, is_t_primer, REMAP_THRESHOLD, MAX_REMAPS, bowtie_threads)

	queue_request = mapping_factory.queue_work(command)
	if queue_request:
		p, command = queue_request
		print "STARTING: {} (pID {})".format(command, p.pid)

while True:
	processes_spawned = mapping_factory.supervise()
	if processes_spawned:
		for tuple in processes_spawned:
			p, command = tuple
			print "STARTING: {} (pID {})".format(command, p.pid)
	if mapping_factory.completely_idle(): break
	time.sleep(1)

# First barrier
