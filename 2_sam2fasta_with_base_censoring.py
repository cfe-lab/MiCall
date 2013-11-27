import os, sys
from miseqUtils import len_gap_prefix, sam2fasta, timestamp

## Arguments
samfile = sys.argv[1]
qCutoff = int(sys.argv[2])
HXB2_mapping_cutoff = int(sys.argv[3])
mode = sys.argv[4]
max_prop_N = float(sys.argv[5])

filename = samfile.split('/')[-1]
prefix, region = filename.split('.')[:2]

# Convert SAM to fasta-structured variable
with open(samfile, 'rU') as infile:
	fasta = sam2fasta(infile, qCutoff, HXB2_mapping_cutoff, max_prop_N)

# Send warning to standard out if sam2fasta didn't return anything
if fasta == None:
	timestamp("WARNING (sam2fasta): {} likely empty or invalid - halting".format(samfile))
	sys.exit()

# For Amplicon runs, write a fasta file to disk
if mode == 'Amplicon':

	# Compress identical sequences
	d = {}
	for h, s in fasta:
		if d.has_key(s):
			d[s] += 1
		else:
			d.update({s: 1})

	# Write fasta to disk, sorted by the prevalence of each unique sequence
	intermed = [(count, s) for s, count in d.iteritems()]
	intermed.sort(reverse=True)
	fasta_filename = '.'.join(map(str,[samfile.replace('.remap.sam', ''), qCutoff, 'fasta']))
	with open(fasta_filename, 'w') as outfile:
		for i, (count, seq) in enumerate(intermed):
			outfile.write('>%s_variant_%d_count_%d\n%s\n' % (prefix, i, count, seq))

# For Nextera runs, write csf file (Our proprietary format) to disk
elif mode == 'Nextera':

	# Sort csf by left-gap prefix: the offset of the read relative to the ref seq
	intermed = [(len_gap_prefix(s), h, s) for h, s in fasta]
	intermed.sort()
	csf_filename = '.'.join(map(str,[samfile.replace('.remap.sam', ''), qCutoff, 'csf']))
	with open(csf_filename, 'w') as outfile:
		for (gp, h, seq) in intermed:
			outfile.write('%s,%d,%s\n' % (h, gp, seq.strip('-')))
