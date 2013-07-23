"""
Collapse identical reads and track of read count.

Fasta headers are renamed to only reference the
sample, as determined by the prefix in the filename
of the fasta file itself (As opposed to the header)

Input: <sample>.<region>.remap.sam.<qScore>.fasta
Output: <sample>.<region>.remap.sam.<qScore>.fasta.<minCount>.seq
"""

import sys
from seqUtils import convert_fasta


f = sys.argv[1]
minCount = int(sys.argv[2])

filename = f.split('/')[-1]
prefix = filename.split('.')[0]


infile = open(f, 'rU')
fasta = convert_fasta(infile.readlines())
infile.close()


# Collect identical sequences
d = {}
for h, s in fasta:
	if d.has_key(s):
		d[s] += 1
	else:
		d.update({s: 1})


# sort sequences by count
intermed = [(count, s) for s, count in d.iteritems()]
intermed.sort(reverse=True)

outfilename = "{}.{}.seq".format(f,minCount)
outfile = open(outfilename, 'w')

for i, (count, seq) in enumerate(intermed):
	if count >= minCount:
		outfile.write('>%s_variant_%d_count_%d\n%s\n' % (prefix, i, count, seq))

outfile.close()
