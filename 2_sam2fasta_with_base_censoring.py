"""
Input: <sample>.<region>.remap.sam files
Output: <sample>.<region>.remap.sam.<qScore>.fasta
"""

import sys
import os

from StringIO import StringIO
from sam2fasta import *
from seqUtils import *

samfile = sys.argv[1]
qCutoff = int(sys.argv[2])

filename = samfile.split('/')[-1]
prefix, region = filename.split('.')[:2]


# convert SAM to FASTA by parsing CIGAR strings, censoring bases and merging reads
infile = open(samfile, 'rU')
fasta_str = sam2fasta(infile, qCutoff)
infile.close()

if fasta_str is None:
	sys.exit()

handle = StringIO(fasta_str)
fasta = convert_fasta(handle.readlines())


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

ofname = '.'.join(map(str, [samfile, qCutoff, 'fasta']))

outfile = open(ofname, 'w')

for i, (count, seq) in enumerate(intermed):
	outfile.write('>%s_variant_%d_count_%d\n%s\n' % (prefix, i, count, seq))

outfile.close()

