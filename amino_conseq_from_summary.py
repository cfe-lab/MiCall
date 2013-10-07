"""
Takes an amino file and outputs a consensus

Assumptions:
	The first 3 columns are Sample, region, and coord respectively
"""

import os,sys
from glob import glob

#summary_file = "/data/miseq/0_testing/HXB2.amino_poly.summary"
summary_file = sys.argv[1]
min_count = int(sys.argv[2])

infile = open(summary_file, 'rU')
lines = infile.readlines()
infile.close()

def set_conseq(sample,region,coord,amino):
	if sample not in amino_conseq: amino_conseq[sample] = {}
	sample_dict = amino_conseq[sample]
	if region not in sample_dict: sample_dict[region] = {}
	sample_region_dict = sample_dict[region]
	if coord in sample_region_dict:
		print "Error! (sample,region,coord) already exists!"
		sys.exit()
	else: sample_region_dict[coord] = amino

# Key the conseq based on sample,region,coord
amino_conseq = {}
for i,line in enumerate(lines):
	fields = line.rstrip("\n").split(',')
	sample, region, coord = fields[:3]

	if i == 0:
		alphabet = fields[3:]
		continue

	freqs = {}
	for j,char in enumerate(alphabet): freqs[char] = int(fields[j+3])
	dominant_amino = max(freqs, key=lambda n: freqs[n])
	max_count = freqs[dominant_amino]

	if max_count < min_count:
		set_conseq(sample,region,coord,'?')
		continue

	# Get all chars that have this count
	matches = filter(lambda n: freqs[n] == max_count, freqs)
	if len(matches) == 1:
		set_conseq(sample,region,coord,matches[0])
		continue

	# If there is a tie, do not report at this point
	else:
		set_conseq(sample,region,coord,'?')
		continue

# Traverse amino_conseq and output the results
outfile = open(summary_file + '.conseq', 'w')
for sample in sorted(amino_conseq.keys()):
        sample_dict = amino_conseq[sample]
	for region in sorted(sample_dict.keys()):
		sample_region_dict = sample_dict[region]
		conseq = ""
		for coord in sorted(sample_region_dict.keys(), key=int):
			conseq += sample_region_dict[coord]
		outfile.write(">{}_{}\n{}\n".format(sample,region,conseq))
outfile.close()
