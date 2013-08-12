"""
Quickly convert FASTQ from MiSeq to FASTA using only 
pairwise alignment to map the amino acid sequence to
the reference. 

*** THIS ASSUMES THAT THE READ IS ALREADY IN FRAME ***

"""

import os
import sys
from Bio import SeqIO
from seqUtils import *
from hyphyAlign import *
from g2pScorer import *



hyphy = HyPhy._THyPhy (os.getcwd(), 1)
change_settings(hyphy) # default settings are for protein alignment
refseq = translate_nuc(refSeqs['V3, clinical'], 0)


try:
	infile = sys.argv[1]
	if not infile.endswith('.fastq'):
		print 'file extension must be .fastq'
		sys.exit()
	QCUTOFF = int(sys.argv[2])
except:
	print 'Usage: python fastq2fasta.py [FASTQ] [QCUTOFF]'
	sys.exit()	


handle = open(infile, 'rU')
fasta = []
for record in SeqIO.parse(handle, 'fastq'):
	header = record.id
	seq = list(record.seq)
	quals = record.letter_annotations["phred_quality"]
	# censor bases with quality < 25
	for i, q in enumerate(quals):
		if q < QCUTOFF:
			seq[i] = 'N'
	seq = ''.join(seq)
	fasta.append([header, seq])



# convert FASTQ to FASTA
prefix = infile.split('/')[-1].split('_')[0] # sample name
print prefix

outfile = infile.replace('.fastq', '.q%d.fasta' % QCUTOFF)
"""
SeqIO.convert(infile, 'fastq', outfile, 'fasta')
fasta = convert_fasta(open(outfile, 'rU').readlines())
"""


# collect identical sequences
d = {}
for h, s in fasta:
	if d.has_key(s):
		d[s] += 1
	else:
		d.update({s: 1})

intermed = [(count, s) for s, count in d.iteritems()]
intermed.sort(reverse=True)

# output nuc variants, translate, align against V3 reference and clip
seqfile = open(outfile.replace('.fasta', '.seq'), 'w')
badfile = open(outfile.replace('.fasta', '.bad'), 'w')
v3nucfile = open(outfile.replace('.fasta', '.v3nuc'), 'w')

v3 = {}
for i, (count, seq) in enumerate(intermed):
	seqfile.write('>%s_variant_%d_count_%d\n%s\n' % (prefix, i, count, seq))
	
	aaseq = translate_nuc(seq, 0)
	aquery, aref, ascore = pair_align(hyphy, refseq, aaseq)
	left, right = get_boundaries(aref)
	
	v3prot = aquery[left:right]
	
	# apply alignment to nucleotide sequence
	v3nuc = apply2nuc(seq[(3*left):], v3prot, aref[left:right], keepIns=True, keepDel=False)
	if 'N' in v3nuc: # we could also screen on '?' in v3prot
		continue
	
	v3nucfile.write('>%s_variant_%d_count_%d\n%s\n' % (prefix, i, count, v3nuc))
	
	# screen for bad V3 sequences, provide reason(s)
	if not v3prot.startswith('C') or not v3prot.endswith('C') or '*' in v3prot or ascore < 50:
		badfile.write('>%s_variant_%d_count_%d_reason_%s\n%s\n' % (prefix, i, count,
			'|'.join(['stopcodon' if '*' in v3prot else '',
					'lowscore' if ascore < 50 else '',
					'cystines' if not v3prot.startswith('C') or not v3prot.endswith('C') else '']),
			seq))
		continue
	
	if v3.has_key(v3prot):
		v3[v3prot] += count
	else:
		v3.update({v3prot: count})
	

v3nucfile.close()
badfile.close()
seqfile.close()


# collect identical V3 amino acid sequences and output
intermed = [(count, s) for s, count in v3.iteritems()]
intermed.sort(reverse=True)

# emulate g2p settings
change_settings(hyphy, 	alphabet=gonnetAlphabet, 
						scoreMatrix=scoreMatrixGonnet,
						gapOpen=10,
						gapOpen2=5,
						gapExtend=1,
						gapExtend2=1,
						noTerminalPenalty=0)

v3file = open(outfile.replace('.fasta', '.v3prot'), 'w')
for i, (count, seq) in enumerate(intermed):
	"""
	if count == 1:
		break
	"""
	aquery, aref, ascore = pair_align(hyphy, g2p_ref, seq)
	
	if '?' in seq:
		v3file.write('>%s_variant_%d_count_%d_fpr_NA\n%s\n' % (prefix, i, count, aquery))
		continue
	
	g2p = calc_g2p(aquery, aref, None)
	
	fpr = convert2fpr(g2p)
	v3file.write('>%s_variant_%d_count_%d_fpr_%1.1f\n%s\n' % (prefix, i, count, fpr, aquery))

v3file.close()


