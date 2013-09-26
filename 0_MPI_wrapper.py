"""
Distribute pipeline processes across the cluster.
"""

import os,sys
from glob import glob
from mpi4py import MPI
from seqUtils import convert_fasta, ambig_dict
from datetime import datetime

my_rank = MPI.COMM_WORLD.Get_rank()
nprocs = MPI.COMM_WORLD.Get_size()

def timestamp(msg): print '[{}] (rank={}/{}) {}'.format(datetime.now().strftime('%Y-%m-%d %H:%M:%S'),my_rank,nprocs,msg)

refpath = '/usr/local/share/miseq/refs/cfe'	# Consensus B refs
root = sys.argv[1]				# Path to folder containing fastqs
mode = sys.argv[2] 				# Nextera or Amplicon - from SampleSheet.csv
qCutoff = int(sys.argv[3])			# INCORRECT - FIXME

# For each fastq file, map sequence to reference sequences: 1_mapping(refseq,files[i])
files = glob(root + '/*R1*.fastq')
for i in range(len(files)):
	if i % nprocs != my_rank: continue
	filename = files[i].split('/')[-1]
	timestamp('1_mapping {}'.format(filename))
	os.system("python 1_mapping.py {} {} {}".format(refpath, files[i], qCutoff))
MPI.COMM_WORLD.Barrier()

if my_rank == 0: timestamp('MPI.COMM_WORLD.Barrier() #1 passed\n')

# For each remapped SAM, generate CSFs (done), generate pileup(Need to be done)
files = glob(root + '/*.remap.sam')
for i in range(len(files)):
	if i % nprocs != my_rank: continue
	filename = files[i].split('/')[-1]

	# For each q-cutoff, generate fasta for amplicon, OR a csf file for Nextera
	for qcut in [0, 10, 15, 20]:
		timestamp('2_sam2fasta {} {} {}'.format(filename, qcut, mode))
		os.system('python 2_sam2fasta_with_base_censoring.py {} {} {}'.format(files[i], qcut, mode))
MPI.COMM_WORLD.Barrier()

if my_rank == 0: timestamp('MPI.COMM_WORLD.Barrier() #2 passed\n')

# Fpr amplicon sequencing runs, compute g2p scores for env-mapped FASTAs
if mode == 'Amplicon':
	files = glob(root + '/*env*.fasta')
	for i, file in enumerate(files):
		if i % nprocs != my_rank: continue
		filename = files[i].split('/')[-1]
		timestamp('3_g2p_scoring {}'.format(filename))
		os.system('python 3_g2p_scoring.py {}'.format(file))
MPI.COMM_WORLD.Barrier()

if my_rank == 0: timestamp('MPI.COMM_WORLD.Barrier() #3 passed\n')

MPI.Finalize()
