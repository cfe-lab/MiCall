"""
Resume pipeline3.py from the sam2fasta stage

Starting with *.remap.sam, generate *.qScore.fasta

Dependencies: sam2fasta
"""

import sys
import os
from glob import glob
from sam2fasta import *

if len(sys.argv) != 3:
	print 'Usage: python pipeline_3_resume_from_sam2fasta.py /path/to/remap_samfiles/ 20'
	sys.exit()

remapPath = sys.argv[1]
cutoff = int(sys.argv[2])
globPath = remapPath + '*.remap.sam'

print "Globbing: {} ...\n".format(globPath)
files = glob(globPath)

for f in files:
	print "sam2fasta({},{})".format(os.path.basename(f), cutoff)
	infile = open(f, 'rU')
	fasta = sam2fasta(infile, cutoff)
	infile.close()
	
	if fasta:
		outfilename = "{}.{}.fasta".format(f, cutoff)
		outfile = open(outfilename, 'w')
		outfile.write(fasta)
		outfile.close()
