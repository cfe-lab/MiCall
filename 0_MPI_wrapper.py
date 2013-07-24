"""
Distribute pipeline processes across the cluster.
"""

import sys
import os
from glob import glob
from mpi4py import MPI

my_rank = MPI.COMM_WORLD.Get_rank()
nprocs = MPI.COMM_WORLD.Get_size()

refpath = '/usr/local/share/miseq/refs/cfe'

# iterative mapping of FASTQs to references
root = sys.argv[1]
files = glob(root + '/*R1*.fastq')

for i in range(len(files)):
	if i % nprocs != my_rank:
		continue
	
	filename = files[i].split('/')[-1]
	print '... process %d of %d starting task 1_mapping.py on %s' % (my_rank, nprocs, filename)
	os.system('python 1_mapping.py %s %s' % (refpath, files[i]))

MPI.COMM_WORLD.Barrier()

# collect all remapped SAM files
files = glob(root + '/*.remap.sam')

for i in range(len(files)):
	if i % nprocs != my_rank:
		continue
	filename = files[i].split('/')[-1]
	print '... process %d of %d starting task 2_sam2fasta on %s' % (my_rank, nprocs, filename)
	for qcut in [0, 10, 15, 20]:
		os.system('python 2_sam2fasta_with_base_censoring.py %s %d' % (files[i], qcut))


# collect all FASTA files
files = glob(root + '/*.fasta')

for i, file in enumerate(files):
	if i % nprocs != my_rank:
		continue
	filename = files[i].split('/')[-1]
        print '... process %d of %d starting task 3_g2p_scoring on %s' % (my_rank, nprocs, filename)
	

MPI.COMM_WORLD.Barrier()
MPI.Finalize()

