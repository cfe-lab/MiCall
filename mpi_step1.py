import os
import sys
from glob import glob
from datetime import datetime
from mpi4py import MPI

my_rank = MPI.COMM_WORLD.Get_rank()
nprocs = MPI.COMM_WORLD.Get_size()

def timestamp(msg):
        print '[%s] (%d/%d) %s' % (datetime.now().strftime('%Y-%m-%d %H:%M:%S'), 
                        my_rank,
                        nprocs,
                        msg)

refpath = '/usr/local/share/miseq/refs/cfe'
root = sys.argv[1]

# From fasta files
files = glob(root + '/*R1*.fastq')

for i, file in enumerate(files):

        if i % nprocs != my_rank:
                continue

        filename = files[i].split('/')[-1]
        timestamp('1_mapping %s' % filename)
        os.system('python 1_mapping.py %s %s 20' % (refpath, file))

# Wait for all threads to complete
MPI.COMM_WORLD.Barrier()
MPI.Finalize()

