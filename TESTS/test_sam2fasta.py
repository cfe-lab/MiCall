import sys
sys.path.append("/usr/local/share/miseq/development/miseqpipeline")

from miseq_modules import sam2fasta_with_base_censoring
from glob import glob

read_mapping_cutoff = 10
sam2fasta_q_cutoffs = [0,10,15]
max_prop_N = 0.5
mode = "Amplicon"

root = "/data/miseq/131112_M01841_0040_000000000-A5F9H"
files = glob(root + '/*.remap.sam')

for i in range(len(files)):
	filename = files[i].split('/')[-1]
	for qcut in sam2fasta_q_cutoffs:
		print "sam2fasta_with_base_censoring({}, {}, {}, {}, {})".format(files[i], qcut, read_mapping_cutoff, mode, max_prop_N)
		sam2fasta_with_base_censoring(files[i], qcut, read_mapping_cutoff, mode, max_prop_N)
