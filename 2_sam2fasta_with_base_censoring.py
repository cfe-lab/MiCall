import os,sys
from miseqUtils import sam2fasta

#python 2_sam2fasta_with_base_censoring.py /data/miseq/131119_M01841_0041_000000000-A5EPY/F00113-V3LOOP_S86.HIV1B-env.remap.sam 0 10 Amplicon
samfile = sys.argv[1]
qCutoff = int(sys.argv[2])
HXB2_mapping_cutoff = int(sys.argv[3])
mode = sys.argv[4]

filename = samfile.split('/')[-1]
prefix, region = filename.split('.')[:2]

infile = open(samfile, 'rU')
fasta = sam2fasta(infile, qCutoff, HXB2_mapping_cutoff, 0.5)
infile.close()

# If Amplicon, generate FASTA
if mode == 'Amplicon':

	# Compress identical sequences
	d = {}
	for h, s in fasta:
		if d.has_key(s):
			d[s] += 1
		else:
			d.update({s: 1})

	# Sort sequences by count and save to a FASTA
	intermed = [(count, s) for s, count in d.iteritems()]
	intermed.sort(reverse=True)
	fasta_filename = '.'.join(map(str,[samfile.replace('.remap.sam', ''), qCutoff, 'fasta']))

	with open(fasta_filename, 'w') as outfile:
		for i, (count, seq) in enumerate(intermed):
			outfile.write('>%s_variant_%d_count_%d\n%s\n' % (prefix, i, count, seq))

# If Nextera, generate CSF	
elif mode == 'Nextera':

	# Sort reads by gap prefix
	intermed = [(len_gap_prefix(s), h, s) for h, s in fasta]
	intermed.sort()
	
	csf_filename = '.'.join(map(str,[samfile.replace('.remap.sam', ''), qCutoff, 'csf']))

	with open(csf_filename, 'w') as outfile:
		for (gp, h, seq) in intermed:
			outfile.write('%s,%d,%s\n' % (h, gp, seq.strip('-')))
