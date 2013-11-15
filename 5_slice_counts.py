"""
Define the bounds we want to extract from a particular region

Label, region to slice, start of slice, end of slice
"""

import os, sys
from glob import glob

#root = sys.argv[1]
root = "/usr/local/share/miseq/development/miseqpipeline/testing"

region_slices = [	("PRRT", "HIV1B-pol", 1, 199),
			("INTEGRASE", "HIV1B-pol", 200, 300)]

# For each rule, get *.<region>.*nuc.csv
for rule in region_slices:
	label, region, start, end = rule
	files = glob(root + '/*.{}.*nuc.csv'.format(region))

	# Slice each file according to the rules!
	for path in files:
		fileName = os.path.basename(path)
		sample,old_region = fileName.split(".")[:2]

		# Read the file in
		f = open(path, 'rU')
		lines = f.readlines()
		f.close()

		# Change the coordinate space
		for line in lines:
			print line.rstrip("\n")

	sys.exit()
