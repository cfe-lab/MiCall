"""
For each .v3 file, enumerate through fasta contents, and determine
proportion X4 for that patient/QC rule, as determined by sequence count
and G2PFPR cutoff.

Input: HIV1B-env.remap.sam.<qCutoff>.fasta.<minCount>.seq.V3
Output: summary.X4 (Contains name of file and percent X4)
"""
from seqUtils import convert_fasta

def prop_x4 (f, cutoff, minCount):
	infile = open(f, 'rU')
	try:
		fasta = convert_fasta(infile.readlines())
	except:
		print 'failed to convert', f
		sys.exit()
	infile.close()

	total_count = 0
	total_x4_count = 0

	for h, s in fasta:
		# >F00309-IL_variant_0_count_27_fpr_4.0
		tokens = h.split('_')
		try:
			variant = int(tokens[tokens.index('variant')+1])
			count = int(tokens[tokens.index('count')+1])
			fpr = float(tokens[tokens.index('fpr')+1])
		except:
			continue
	
		if count < minCount: continue
		if fpr < cutoff: total_x4_count += count
		total_count += count

	return (total_x4_count, total_count)
