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

root = sys.argv[1]
files = glob(root + '/*R1*.fastq')
print len(files), 'files'

for i in range(len(files)):
	if i % nprocs != my_rank:
		continue
	
	filename = files[i].split('/')[-1]
	print '... process %d of %d starting task 1_mapping.py on %s' % (my_rank, nprocs, filename)
	os.system('python 1_mapping.py %s %s' % (refpath, files[i]))

MPI.COMM_WORLD.Barrier()
MPI.Finalize()

