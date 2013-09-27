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

root = sys.argv[1]
mode = sys.argv[2]

# From fasta files
files = glob(root + '/*.fasta' if mode == 'Amplicon' else root+'/*.csf')

for i, file in enumerate(files):

        if i % nprocs != my_rank:
                continue

        filename = files[i].split('/')[-1]
        timestamp('5_amino_freqs %s' % filename)
        os.system('python 5_amino_freqs.py %s %s' % (file, mode))

# Wait for all threads to complete
MPI.COMM_WORLD.Barrier()
MPI.Finalize()
