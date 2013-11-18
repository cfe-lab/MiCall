"""
Slice 4_csf2counts.py output (Count CSVs + conseq files) into sub-regions.
"""

import os, sys
from glob import glob
from miseqUtils import convert_fasta

# Define slices with (Label, region to slice, start, end)
# Coordinates are in nucleotide space: start/end are inclusive
region_slices = [	("PROTEASE", "HIV1B-pol", 1, 297),
			("PRRT", "HIV1B-pol", 1, 1617),
			("INTEGRASE", "HIV1B-pol", 1977, 2843)]

if len(sys.argv) != 2:
	print 'Usage: python {} /path/to/outputs'.format(sys.argv[0])
	sys.exit()

root = sys.argv[1]

for rule in region_slices:
	slice, region, start, end = rule

	# STEP 1: CONVERT NUC/AMINO COUNTS
	files = glob(root + '/*.{}.*nuc.csv'.format(region))
	files += glob(root + '/*.{}.*amino.csv'.format(region))

	for path in files:
		fileName = os.path.basename(path)
		sample,old_region = fileName.split(".")[:2]

		with open(path, 'rU') as f:
			lines = f.readlines()

		newFileName = fileName.replace(region,slice)

		dirName = os.path.dirname(path)
		f = open("{}/{}".format(dirName, newFileName), 'w')
		for i,line in enumerate(lines):
			line = line.rstrip("\n")

			# Retain CSV header
			if (i == 0):
				f.write("{}\n".format(line))
				continue

			query_pos, hxb2_pos = map(int, line.split(",")[:2])

			# Change coordinates, only displaying regions of interest
			if "nuc.csv" in path:
				if hxb2_pos < start or hxb2_pos > end: continue
				region_pos = hxb2_pos - start + 1

			elif "amino.csv" in path:
				if hxb2_pos < start/3 or hxb2_pos > end/3: continue
				region_pos = hxb2_pos - start/3

			counts = line.split(",")[2:]
			f.write("{},{},{}\n".format(query_pos,region_pos,",".join(counts)))
		f.close()

	# STEP 2: SLICE CONSEQ FILES
	files = glob(root + '/*.{}.*conseq'.format(region))
	for path in files:
		with open(path, 'rU') as f:
			fasta = convert_fasta(f.readlines())

		fileName = os.path.basename(path)
		dirName = os.path.dirname(path)
		newFileName = fileName.replace(region,slice)
		f = open("{}/{}".format(dirName, newFileName), 'w')
		for j, (h,s) in enumerate(fasta):
			f.write(">{}\n{}\n".format(h, s[start:end+1]))
		f.close()
