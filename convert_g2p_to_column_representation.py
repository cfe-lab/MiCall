import sys
import os

g2p_file = sys.argv[1]
infile = open(g2p_file, 'rU')

propX4s = {}
QC_filters = set()
samples = set()

for line in infile:
	#sample,qCut,fprCut,minCount,x4Count,totalCount,prop.x4
	line = line.replace('\n','')
	sample, qCut, fprCut, minCount, x4Count, totalCount, propX4 = line.split(',')

	if sample == 'sample':
		continue

	QC_filter = "{}-{}-{}".format(qCut, minCount, fprCut[:4])

	# Add the propX4 for this (enum,filter)
	if sample not in propX4s.keys():
		propX4s[sample] = {}
		propX4s[sample][QC_filter]=propX4
	else:
		propX4s[sample][QC_filter]=propX4

	# Track all samples and filter types
	samples.add(sample)
	QC_filters.add(QC_filter)

infile.close()

samples = sorted(samples)
QC_filters = sorted(QC_filters)

# Give a legend
print "Filter legend: qCut-minCount-fprCut"

# Print the header
sys.stdout.write('sample,')
for QC_filter in QC_filters:
	sys.stdout.write(QC_filter + ",")
print ""

# For each sample, retrieve propX4 for each filter
for sample in samples:
	sys.stdout.write(sample + ",")
	for QC_filter in QC_filters:
		if sample in propX4s.keys() and QC_filter in propX4s[sample].keys():
			sys.stdout.write(propX4s[sample][QC_filter] + ",")
		else:
			sys.stdout.write('NA,')
	print ""
