"""
Input: <sample>.<region>.remap.sam files
Output: <sample>.<region>.remap.sam.<qScore>.fasta
"""

import sys
import os
from glob import glob
from sam2fasta import *

if len(sys.argv) != 3:
	print 'Usage: python 2_sam2fasta_with_censoring.py 130621_M01841_0006_000000000-A3RWC/Data/Intensities/BaseCalls/remap_sams/ 20'
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
