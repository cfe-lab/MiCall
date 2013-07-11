"""
2_pipeline_resume_from_sam2fasta.py

Input: <sample>.<region>.remap.sam files
Output: <sample>.<region>.remap.sam.<qScore>.fasta

Followup step: 3_compress_fasta.py to generate <sample>.<region>.remap.<qScore>.fasta.<minCount>.seq

Dependencies: sam2fasta
"""

import sys
import os
from glob import glob
from sam2fasta import *

if len(sys.argv) != 3:
	print 'Usage: python 2_pipeline_resume_from_sam2fasta.py /path/to/remap_samfiles/ 20'
	sys.exit()

remapPath = sys.argv[1]
qCutoff = int(sys.argv[2])
globPath = remapPath + '*.remap.sam'

print "Globbing: {} ...\n".format(globPath)
files = glob(globPath)

for f in files:
	print "sam2fasta({},{})".format(os.path.basename(f), qCutoff)
	infile = open(f, 'rU')
	fasta = sam2fasta(infile, qCutoff)
	infile.close()
	
	if fasta:
		outfilename = "{}.{}.fasta".format(f, qCutoff)
		outfile = open(outfilename, 'w')
		outfile.write(fasta)
		outfile.close()
