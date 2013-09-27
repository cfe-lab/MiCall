"""
Because of missing Barrier() in 0_MPI step, some HIV env fastas were generated
before they could be catalogued by a glob statement.  As a result, they were
never converted into *.v3prot format.
"""

from glob import glob
import os
from mpi4py import MPI
from datetime import datetime

my_rank = MPI.COMM_WORLD.Get_rank()
nprocs = MPI.COMM_WORLD.Get_size()

def timestamp(msg):
	print '[%s] (%d/%d) %s' % (datetime.now().strftime('%Y-%m-%d %H:%M:%S'), 
			my_rank,
			nprocs,
			msg)

# get a list of all HIV1B-env fastas
paths = glob('/data/miseq/*/*.HIV1B-env.*.fasta')
todo = [path for path in paths if not os.path.exists(path.replace('.fasta', '.v3prot'))]


for i, path in enumerate(todo):
	if i % nprocs != my_rank:
		continue
	filename = path.split('/')[-1]
	timestamp(filename)
	os.system('python /usr/local/share/miseq/scripts/3_g2p_scoring.py %s' % path)

MPI.COMM_WORLD.Barrier()
MPI.Finalize()

