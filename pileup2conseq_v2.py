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

ambig_dict['ACGT'] = 'N'


misc = re.compile('\^.|\*+|\$+')
indels_re = re.compile('\+[0-9]+|-[0-9]+')

# mpileup generates a massive file
report_indels = True if '-indels' in sys.argv else False

infile = open(sys.argv[1], 'rU')
#polyfile = open(sys.argv[1]+'.poly', 'w')


prefix = sys.argv[1].split('/')[-1].split('_')[0]


if report_indels:
	indelfile = open(sys.argv[1]+'.indels', 'w')
	indelfile.write('Position,indel,count\n')


#polyfile.write('Position,count.A,count.C,count.G,count.T\n')

conseq = ''

freq_minor = [0.25, 0.1, 0.05, 0.01]
conseq_minor = ['', '', '', '']

#start = False
to_skip = 0
for ii, line in enumerate(infile):
	if to_skip > 0:
		# account for a majority deletion in previous lines
		to_skip -= 1
		continue
	
	label, pos, en, depth, astr, qstr = line.strip('\n').split('\t')
	pos = int(pos)
	
	# convert strings into lists
	alist = []
	qlist = []
	i = 0 # index astr
	j = 0 # indel qstr
	while i < len(astr):
		if astr[i] == '^':
			# '^' marks the start of a read, e.g. : "^7G"
			# the ASCII code of the next character minus 33 gives 
			# 	the READ mapping quality
			# the next character is the start base
			q = ord(qstr[j])-33
			if q >= 30:
				alist.append(astr[i+2])
			else:
				alist.append('N')
			qlist.append(q)
			i += 3
			j += 1
			
		elif astr[i] in '*$':
			# '*' represents a deleted base
			# '$' indicates the end of a read
			#qlist.append(None)
			i += 1
			
		else:
			# look ahead for indels at this position
			if i < len(astr)-1 and astr[i+1] in '+-':
				# returns match at start of string
				m = indels_re.match(astr[i+1:])
				
				# number of characters to look ahead
				indel_len = int(m.group().strip('+-'))
				left = i+1 + len(m.group())
				insertion = astr[left:(left+indel_len)]
				
				q = ord(qstr[j])-33
				base = astr[i].upper() if q >= 30 else 'N'
				
				token = base + m.group() + insertion
				alist.append(token)
				qlist.append(q)
				
				# update indices
				i += len(token)
				j += 1
				
			else:
				# no indel ahead
				q = ord(qstr[j])-33
				base = astr[i].upper() if q >= 30 else 'N'
				alist.append(base)
				qlist.append(q)
				j += 1
				i += 1
	
	
	# is this dominated by an indel?
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
	
	counts = [(nuc, alist.count(nuc)) for nuc in 'ACGT']
	
	intermed = [(v,k) for k, v in counts]
	intermed.sort(reverse=True)
	
	conseq += intermed[0][1]


confile = open(sys.argv[1]+'.conseq', 'w')
confile.write('>%s\n%s\n' % (prefix, conseq)) 
confile.close()



infile.close()





