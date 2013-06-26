"""
determine_proportion_X4.py

For each .v3 file, enumerate through fasta contents, and determine
proportion X4 for that patient/QC rule, as determined by sequence count
and G2PFPR cutoff.

Input: HIV1B-env.remap.sam.<qCutoff>.fasta.<minCount>.seq.V3
Output: summary.X4 (Contains name of file and percent X4)
"""

import os
import sys
import re
from glob import glob
from seqUtils import convert_fasta

helpOutput = """Usage: python determine_proportion_X4.py <G2PFPR cutoff> <folderContainingV3Files>
Example (X4 if G2P <= 3.5): python determine_proportion_X4.py 3.5 ../path/to/files/"""

if len(sys.argv) != 3:
	print helpOutput
	sys.exit()

G2P_FPR_cutoff = float(sys.argv[1])
globPath = sys.argv[2] + '*.HIV1B-env.remap.sam.*.fasta.*.seq.v3'

print "filename,sample,qCutoff,minCount,G2PFPRcutoff,total_X4_seqs,total_seqs,proportion_X4"


files = glob(globPath)

# For each v3 file containing G2P scores, determine proportion X4
for f in files:

	# Parse the patient information out of the filename
	# F00142.HIV1B-env.remap.sam.0.fasta.50.seq.v3
	filename = os.path.basename(f)
	filenameFields = re.split('[.]', filename)
	sample =filenameFields[0]
	qCutoff = filenameFields[4]
	minCount = filenameFields[6]

	infile = open(f, 'rU')
	try:
		fasta = convert_fasta(infile.readlines())
	except:
		print 'failed to convert', f
		continue
	infile.close()

	total_seqs = 0
	total_X4_seqs = 0

	# For each (header, sequence), extract V3
	for header, v3Seq in fasta:

		# Drop fields for which GP2_FPR is nill (Could not be aligned in Conan's G2P)
		try:
			G2P_FPR = float(re.split('G2PFPR_([0-9]+[.][0-9])+', header)[1])
			count = int(re.split('count_([0-9]+)',header)[1])

			if G2P_FPR < G2P_FPR_cutoff:
				total_X4_seqs += count

			total_seqs += count

		except:
			continue

	proportion_X4 = float(total_X4_seqs) / float(total_seqs)
	print filename+","+sample+","+qCutoff+","+minCount+","+str(G2P_FPR_cutoff)+"," +str(total_X4_seqs)+","+str(total_seqs)+","+str(proportion_X4)
