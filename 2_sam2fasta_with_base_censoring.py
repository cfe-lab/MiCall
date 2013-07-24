"""
Input: <sample>.<region>.remap.sam files
Output: <sample>.<region>.remap.sam.<qScore>.fasta
"""

import sys
import os

from sam2fasta import *


samfile = sys.argv[1]
qCutoff = int(sys.argv[2])


infile = open(samfile, 'rU')
fasta = sam2fasta(infile, qCutoff)
infile.close()

if fasta is None:
	sys.exit()

# for now, do not reject sequences based on the proportion that are 'N'
outfilename = "{}.{}.fasta".format(samfile, qCutoff, max_prop_N=1.0)
outfile = open(outfilename, 'w')
outfile.write(fasta)
outfile.close()
