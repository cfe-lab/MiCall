"""
A library of functions for converting from a samtools SAM Format to a
FASTA formatted alignment.
"""

import re

cigar_re = re.compile('[0-9]+[MIDNSHPX=]')


def apply_cigar (cigar, seq, qual):
	"""
	Parse CIGAR info and apply it to the nucleotide sequence.
	"""
	newseq = ''
	newqual = ''
	tokens = cigar_re.findall(cigar)
	if len(tokens) == 0:
		return None, None, None
	
	shift = 0 # to account for removing soft clipped bases on left
	if tokens[0].endswith('S'):
		shift = int(tokens[0][:-1])
	
	left = 0
	for token in tokens:
		length = int(token[:-1])
		if token[-1] == 'M':
			# carry over matches
			newseq += seq[left:(left+length)]
			newqual += qual[left:(left+length)]
			left += length
			
		elif token[-1] == 'D':
			# deletion relative to reference, pad with gaps
			newseq += '-'*length
			newqual += 'A' # arbitrary high quality score
			
		elif token[-1] == 'I':
			# insertion relative to reference, excise
			left += length
			continue
			
		elif token[-1] == 'S':
			# ignore soft clipping
			left += length
			continue
			
		else:
			print 'how to handle CIGAR token %s?' % token
			sys.exit()
		
	return shift, newseq, newqual


def censor_bases (seq, qual, cutoff):
	"""
	For each base in a nucleotide sequence and quality string,
	replace a base with an ambiguous character 'N' if its associated
	quality score falls below a threshold value.
	"""
	newseq = ''
	for i, q in enumerate(qual):
		if ord(q)-33 >= cutoff:
			newseq += seq[i]
		else:
			newseq += 'N'
	return newseq


def merge_pairs (seq1, seq2):
	"""
	Merge two sequences that overlap over some portion (paired-end
	reads).  Using the positional information in the SAM file, we will
	know where the sequences lie relative to one another.  In the case 
	that the base in one read has no complement in the other read 
	(in partial overlap region), take that base at face value.
	"""
	mseq = ''
	for i, c2 in enumerate(seq2):
		if i < len(seq1):
			c1 = seq1[i]
			if c1 == c2:
				mseq += c1
			elif c1 in 'ACGT':
				if c2 in 'N-':
					mseq += c1
				else:
					mseq += 'N' # error
			elif c2 in 'ACGT':
				if c1 in 'N-':
					mseq += c2
				else:
					mseq += 'N'
			else:
				mseq += 'N'
		else:
			# past extent of seq1
			mseq += c2 
	return mseq


def sam2fasta (infile, cutoff):
	"""
	Parse SAM file contents and return FASTA string.  In the case of
	a matched set of paired-end reads, merge the reads together into
	a single sequence.
	"""
	fasta = ''
	lines = infile.readlines()
	
	# advance to first data line
	for start, line in enumerate(lines):
		if not line.startswith('@'):
			break
	
	# abort if file is empty
	if start == (len(lines)-1):
		return None
	
	i = start
	while i < len(lines):
		qname, flag, refname, pos, mapq, cigar, rnext, pnext, tlen, seq, qual = lines[i].strip('\n').split('\t')[:11]
		
		if refname == '*' or cigar == '*':
			# first read is unmapped
			i += 1
			continue
		
		pos1 = int(pos)
		shift, seq1, qual1 = apply_cigar(cigar, seq, qual)
		if not seq1:
			i += 1
			continue
		
		seq1 = '-'*pos1 + censor_bases(seq1, qual1, cutoff)
		
		if (i+1) == len(lines):
			# no more lines
			break
		
		# look ahead for matched pair
		qname2, flag, refname, pos, mapq, cigar, rnext, pnext, tlen, seq, qual = lines[i+1].strip('\n').split('\t')[:11]
		
		if qname2 == qname:
			# this is the second read
			
			if refname == '*' or cigar == '*': 
				# second read is unmapped
				fasta += '>%s\n%s\n' % (qname, seq1)
				i += 2
				continue
			
			pos2 = int(pos)
			shift, seq2, qual2 = apply_cigar(cigar, seq, qual)
			if not seq2:
				# failed to parse CIGAR string
				fasta += '>%s\n%s\n' % (qname, seq1)
				i += 2
				continue
			
			seq2 = '-'*pos2 + censor_bases(seq2, qual2, cutoff)
			
			mseq = merge_pairs(seq1, seq2)
			if mseq.count('N') / float(len(mseq)) < 0.5:
				# output only if sequence is good quality
				fasta += '>%s\n%s\n' % (qname, mseq)
			
			i += 2
			continue
		
		# ELSE no matched pair
		fasta += '>%s\n%s\n' % (qname, seq1)
		i += 1
	
	return fasta
