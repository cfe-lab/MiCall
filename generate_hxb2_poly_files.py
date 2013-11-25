import os, sys
from sam2fasta import apply_cigar
from seqUtils import translate_nuc, timestamp

def generate_nuc_counts(HXB2_sam, mapping_cutoff=10, debug=0):
	"""
	Generates a .HXB2.nuc_poly file from a .HXB2.sam
	"""
	nuc_alphabet = 'GATCN-'
	
	# parse the filename
	prefix, region, qCut = HXB2_sam.split('.')[:3]
	sample = prefix.split("/")[-1]
	infile = open(HXB2_sam, 'rU')
	lines = infile.readlines()
	infile.close()

	# Determine where the SAM header ends - then for each read in the SAM, tally nucleotide counts
	nucs = {}
	for start, line in enumerate(lines):
		if not line.startswith('@'): break
	j = start

	while j < len(lines):
		qname, flag, refname, pos, mapq, cigar, rnext, pnext, tlen, seq, qual = lines[j].strip('\n').split('\t')[:11]
		shift, seq_aligned, qual_aligned = apply_cigar(cigar, seq, qual)

		# Check the bitwise flag and filter out unmapped sequences
		if int(flag) & 4 == 4 or int(mapq) < mapping_cutoff:
			j += 1
			continue

		# Iterate through the nucleotide seq and tally (coordinate,nuc) counts
		for coord, nuc in enumerate(seq_aligned):
			if not nucs.has_key(coord): nucs.update({coord:{}})
			if not nucs[coord].has_key(nuc): nucs[coord].update({nuc: 0})
			nucs[coord][nuc] += 1
		j += 1

	# For each sample, write nucleotide frequency data to a .nuc_poly file
	nuc_poly = open(HXB2_sam.replace('.HXB2.sam','.HXB2.nuc_poly'), 'w')
	keys = nucs.keys()
	keys.sort()
	nuc_poly.write('Sample,region,nuc.pos,' + ','.join(nuc_alphabet) + '\n')
	for k in keys:
		nuc_poly.write("{},{},{}".format(sample,region,str(k+1)))
		for nuc in nuc_alphabet: nuc_poly.write(",{}".format(nucs[k].get(nuc,0)))
		nuc_poly.write('\n')
	nuc_poly.close()



def generate_amino_counts(HXB2_sam,mapping_cutoff=10,debug=0):
	"""
	Generates a .HXB2.amino_poly file from a .HXB2.sam
	"""
	amino_alphabet = 'ACDEFGHIKLMNPQRSTVWY*-'
	prefix, region, qCut = HXB2_sam.split('.')[:3]
	sample = prefix.split("/")[-1]
	infile = open(HXB2_sam, 'rU')
	lines = infile.readlines()
	infile.close()

	# Skip the SAM header, and for each SAM line, generate amino seq + tally aa counts
	aminos = {}
	for start, line in enumerate(lines):
		if not line.startswith('@'): break

	j = start

	while j < len(lines):

		qname, flag, refname, pos, mapq, cigar, rnext, pnext, tlen, seq, qual = lines[j].strip('\n').split('\t')[:11]
		shift, seq_aligned, qual_aligned = apply_cigar(cigar, seq, qual)

		# Check the bitwise flag and filter out unmapped sequences
		if int(flag) & 4 == 4 or int(mapq) < mapping_cutoff:
			j += 1
			continue

		# Determine ORF offset, translate on the correct ORF, add a left offset to the amino seq
		ORF = (int(pos)-1) % 3
		p = translate_nuc(seq_aligned, ORF)
		amino_offset = (int(pos)-1)/3
		p = "{}{}".format("?"*amino_offset,p)

                # Iterate through the amino seq and tally the (coordinate,aminos) counts
		for coord, aa in enumerate(p):
			if aa in '?': continue
			if not aminos.has_key(coord): aminos.update({coord:{}})
			if not aminos[coord].has_key(aa): aminos[coord].update({aa: 0})
			aminos[coord][aa] += 1
		j += 1

	# For each sample, write amino frequency data to a .amino_poly file
	HXB2_poly_filename = HXB2_sam.replace('.HXB2.sam','.HXB2.amino_poly')
	amino_poly = open(HXB2_poly_filename, 'w')

	timestamp("Writing amino_poly {}".format(HXB2_poly_filename))

	keys = aminos.keys()
	keys.sort()
	amino_poly.write('Sample,region,AA.pos,' + ','.join(amino_alphabet) + '\n')
	for k in keys:
		if debug == 1:sys.stdout.write("{},{},{}".format(sample,region,str(k+1)))
		amino_poly.write("{},{},{}".format(sample,region,str(k+1)))
		for aa in amino_alphabet:
			if debug == 1:sys.stdout.write(",{}".format(aminos[k].get(aa,0)))
			amino_poly.write(",{}".format(aminos[k].get(aa,0)))
		if debug == 1:sys.stdout.write('\n')
		amino_poly.write('\n')
	amino_poly.close()
