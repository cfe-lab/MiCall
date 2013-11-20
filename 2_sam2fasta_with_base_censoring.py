"""
Input: <sample>.<region>.remap.sam files
Output: <sample>.<region>.remap.sam.<qScore>.(fasta|csf)
"""

import sys
import os
from miseqUtils import *

samfile = sys.argv[1]
qCutoff = int(sys.argv[2])
HXB2_mapping_cutoff = int(sys.argv[3])
mode = sys.argv[4]

filename = samfile.split('/')[-1]
prefix, region = filename.split('.')[:2]


# convert SAM to FASTA by parsing CIGAR strings, censoring bases and merging reads
infile = open(samfile, 'rU')

# Note: sam2fasta has a 4th argument: the proportion of N's needed before data are dropped (Default = 0.5)
fasta = sam2fasta(infile, qCutoff, HXB2_mapping_cutoff)
infile.close()

if fasta is None:
	sys.exit()


if mode == 'Amplicon':
	# collate identical sequences
	d = {}
	for h, s in fasta:
		if d.has_key(s):
			d[s] += 1
		else:
			d.update({s: 1})
	
	# sort sequences by count
	intermed = [(count, s) for s, count in d.iteritems()]
	intermed.sort(reverse=True)
	
	ofname = '.'.join(map(str, [samfile.replace('.remap.sam', ''), qCutoff, 'fasta']))
	outfile = open(ofname, 'w')
	for i, (count, seq) in enumerate(intermed):
		outfile.write('>%s_variant_%d_count_%d\n%s\n' % (prefix, i, count, seq))
	outfile.close()
	
elif mode == 'Nextera':
	# sort reads by gap prefix
	intermed = [(len_gap_prefix(s), h, s) for h, s in fasta]
	intermed.sort()
	
	# output in compact comma-separated FASTA
	ofname = '.'.join(map(str, [samfile.replace('.remap.sam', ''), qCutoff, 'csf']))
	outfile = open(ofname, 'w')
	for (gp, h, seq) in intermed:
		outfile.write('%s,%d,%s\n' % (h, gp, seq.strip('-')))
	
	outfile.close()
