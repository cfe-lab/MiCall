"""
Distribute pipeline processes across the cluster.
"""

import sys
from glob import glob
from mpi4py import MPI

my_rank = MPI.COMM_WORLD.Get_rank()
nprocs = MPI.COMM_WORLD.Get_size()

refpath = '/usr/local/share/miseq/refs/cfe'

root = sys.argv[1]
files = glob(root + '/*R1*.fastq')

for i in range(len(files)):
	if i % nprocs != my_rank:
		continue
	
	filename = files[i].split('/')[-1]
	sys.stdout('process %d of %d starting task 1_mapping.py on %s' % filename)
	os.system('python 1_mapping.py %s %s' % (refpath, files[i])
	sys.stdout('process %d of %d completed task 1_mapping.py on %s' % filename)


MPI.COMM_WORLD.Barrier()
MPI.Finalize()

