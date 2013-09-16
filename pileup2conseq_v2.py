"""
Generate a consensus sequence from samtools mpileup output.
Incorporate insertions where they occur reproducibly.

Apparently quality strings do not include scores for indels:

>>> qstr.strip('\n')[-20:]
'GFFDHHHFGHF5H1GH5GHC'
>>> ''.join(map(lambda x: chr(x+33), qlist[-20:]))
'GFFDHHHFGHF5H1GH5GHC'
"""

import sys
import re
from seqUtils import ambig_dict

qCutoff = int(sys.argv[2])
ambig_dict['ACGT'] = 'N'
misc = re.compile('\^.|\*+|\$+')
indels_re = re.compile('\+[0-9]+|-[0-9]+')	# Matches a '+' or '-' followed by 1+ numbers

# mpileup generates a massive file
# A pileup file stores each aligned read for each position with respect to reference
report_indels = True if '-indels' in sys.argv else False
infile = open(sys.argv[1], 'rU')
polyfile = open(sys.argv[1]+'.poly', 'w')
polyfile.write('Position,A,C,G,T,N\n')

# Get the prefix from the pileup file
prefix = sys.argv[1].split('/')[-1].split('_')[0]

if report_indels:
	indelfile = open(sys.argv[1]+'.indels', 'w')
	indelfile.write('Position,indel,count\n')

conseq = ''
freq_minor = [0.25, 0.1, 0.05, 0.01]
conseq_minor = ['', '', '', '']
to_skip = 0

# For each line in the pileup (For a given coordinate in the reference)
for ii, line in enumerate(infile):

	# Account for majority deletion in previous lines
	if to_skip > 0:
		to_skip -= 1
		continue

	# Extract out pileup features
	label, pos, en, depth, astr, qstr = line.strip('\n').split('\t')
	pos = int(pos)
	
	alist = []	# alist stores all bases at a given coordinate
	qlist = []
	i = 0		# Current index for astr
	j = 0		# Current indel for qstr

	# For each position in astr (The main feature list in the pileup)
	while i < len(astr):

		if astr[i] == '^':
			# '^' marks the start of a new read. Ex: "^7G" means a read starts
			# with the first base of 'G' with a quality character of '7'
			# ASCII code of the quality character minus 33 gives the Q-score
			q = ord(qstr[j])-33

			if q >= qCutoff:
				alist.append(astr[i+2])
			else:
				alist.append('N')
			qlist.append(q)

			# Traverse 3 characters in astr
			i += 3
			j += 1

		elif astr[i] in '*$':
			# '*' represents a deleted base
			# '$' indicates the end of a read
			i += 1

		else:
			# Look ahead for insertion/deletions relative to the reference in astr
			if i < len(astr)-1 and astr[i+1] in '+-':
				# returns match at start of string
				m = indels_re.match(astr[i+1:])
				
				# number of characters to look ahead
				indel_len = int(m.group().strip('+-'))
				left = i+1 + len(m.group())
				insertion = astr[left:(left+indel_len)]
				
				q = ord(qstr[j])-33
				base = astr[i].upper() if q >= qCutoff else 'N'
				
				token = base + m.group() + insertion
				alist.append(token)
				qlist.append(q)
				
				# update indices
				i += len(token)
				j += 1
				
			else:
				# no indel ahead
				q = ord(qstr[j])-33
				base = astr[i].upper() if q >= qCutoff else 'N'
				alist.append(base)
				qlist.append(q)
				j += 1
				i += 1
	
	# Is this dominated by an insertion or deletion?
	insertions = [x for x in alist if '+' in x]
	deletions = [x for x in alist if '-' in x]
	non_indel = sum([alist.count(nuc) for nuc in 'ACGT'])
	
	if len(insertions) > non_indel:
		intermed = [(insertions.count(token), token) for token in set(insertions)]
		intermed.sort(reverse=True)
		
		# add most frequent insertion to consensus
		count, token = intermed[0]
		m = indels_re.findall(token)[0] # \+[0-9]+
		
		conseq += token[0] + token[1+len(m):]
		continue
		
	if len(deletions) > non_indel:
		# skip this line and the next N lines as necessary
		intermed = [(deletions.count(token), token) for token in set(deletions)]
		intermed.sort(reverse=True)
		count, token = intermed[0]
		
		m = indels_re.findall(token)[0]
		to_skip = int(m.strip('-')) - 1 # omitting this line counts as one
		continue

	# For this coordinate (line in the pileup), alist now contains all characters that occured
	counts = [(nuc, alist.count(nuc)) for nuc in 'ACGTN']

	# Write results to the poly file
	polyfile.write("{},{},{},{},{},{}\n".format(
		ii,alist.count('A'), alist.count('C'), alist.count('G'), alist.count('T'),alist.count('N')))

	# Store in intermed so we can take the majority base
	intermed = [(v,k) for k, v in counts]
	intermed.sort(reverse=True)	
	conseq += intermed[0][1]

polyfile.close()

confile = open(sys.argv[1]+'.conseq', 'w')
confile.write(">{}\n{}\n".format(prefix, conseq)) 
confile.close()

infile.close()
