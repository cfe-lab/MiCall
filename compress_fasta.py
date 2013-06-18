"""
Collapse all identical reads and keep track of read count.
Output new FASTA with extension *.seq
"""

import sys
from seqUtils import convert_fasta

files = sys.argv[1:]


if len(files) == 0:
	print 'usage: python compress_fasta.py (files)'
	print 'you can use $(ls PATH)'
	sys.exit()

for f in files:
	filename = f.split('/')[-1]
	prefix = filename.split('.')[0]
	
	infile = open(f, 'rU')
	try:
		fasta = convert_fasta(infile.readlines())
	except:
		print 'failed to convert', f
		continue
	
	infile.close()
	
	# collect identical sequences
	d = {}
	for h, s in fasta:
		if d.has_key(s):
			d[s] += 1
		else:
			d.update({s: 1})
	
	intermed = [(count, s) for s, count in d.iteritems()]
	intermed.sort(reverse=True)
	
	outfile = open(f+'.seq', 'w')
	
	for i, (count, seq) in enumerate(intermed):
		outfile.write('>%s_variant_%d_count_%d\n%s\n' % (prefix, i, count, seq))
	
	outfile.close()
