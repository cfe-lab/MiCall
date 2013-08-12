import sys, HyPhy, re, math
import random

def convert_fasta (lines):	
	blocks = []
	sequence = ''
	for i in lines:
		if i[0] == '$': # skip h info
			continue
		elif i[0] == '>' or i[0] == '#':
			if len(sequence) > 0:
				blocks.append([h,sequence])
				sequence = ''	# reset containers
				h = i.strip('\n')[1:]
			else:
				h = i.strip('\n')[1:]
		else:
			sequence += i.strip('\n')
	try:
		blocks.append([h,sequence])	# handle last entry
	except:
		print lines
		raise
	return blocks


def fasta2phylip (fasta, handle):
	ntaxa = len(fasta)
	nsites = len(fasta[0][1])
	handle.write(str(ntaxa)+' '+str(nsites)+'\n')
	for row in fasta:
		# phylip format uses space delimiters
		header = regex.sub('',row[0]).replace(' ','_')
		handle.write(header+' '+row[1]+'\n')


def convert_phylip (lines):
	"""
	Convert line input from Phylip format file into
	Python list object.
	"""
	res = []
	try:
		ntaxa, nsites = lines[0].strip('\n').split()
	except:
		print lines[0]
		raise
		
	if len(lines) != int(ntaxa) + 1:
		raise AssertionError ('Number of taxa does not equal header')
	
	for line in lines[1:]:
		header, seq = line.strip('\n').split()
		res.append( [header, seq] )
	
	return res


def import_seqs (hyphy, path_to_in):
	# use HyPhy file handler to import sequences
	
	#dump = hyphy.ExecuteBF("DataSet ds = ReadDataFile("+os.getcwd()+'/'+path_to_in+");", False)
	dump = hyphy.ExecuteBF("DataSet ds = ReadDataFile("+path_to_in+");", False)
	#dump = hyphy.ExecuteBF("DataSet ds = ReadDataFile(\"/Users/art/wip/StanfordV3/_from_abi.fas\");", False)
	dump = hyphy.ExecuteBF("DataSetFilter dsf = CreateFilter(ds,1);", False)
	dump = hyphy.ExecuteBF("GetInformation(seqs, dsf)", False)
	dump = hyphy.ExecuteBF("GetString(hs, dsf, -1);", False)
	
	dump = hyphy.ExecuteBF('fprintf(stdout, seqs);', False)
	hyout = hyphy.GetStdout()
	seqs = hyout.sData.replace('\n','').replace('{','').replace('}','').replace('"','').replace('?','').split(',')
	
	dump = hyphy.ExecuteBF('fprintf(stdout, hs);', False)
	hyout = hyphy.GetStdout()
	#hs = hyout.sData.replace('\n','').replace('{','').replace('}','').replace('"','').split(',')
	try:
		hs = eval(hyout.sData.replace('{','[').replace('}',']'))
	except:
		print 'ERROR: File empty; did you specify an absolute path?'
		raise
	hs = hs[0]
	
	# re-merge hs and sequences as a dictionary
	sd = {}
	for i in range(len(hs)):
		if 'Standard' in hs[i]:	# drop the consensus sequence from Conan's pipeline
			continue
		try:
			# remove prefix and suffix gaps 
			# 	and in-frame deletions
			#sd.update({hs[i]:{'rawseq':seqs[i].strip('-').replace('---','')}})
			sd.update({hs[i]:{'rawseq':seqs[i]}})
		except:
			print "ERROR: Number of sequences does not match number of hs!"
			print len(hs)
			print hs
			print len(seqs)
			raise
	
	return sd

complement_dict = {'A':'T', 'C':'G', 'G':'C', 'T':'A', 
					'W':'S', 'R':'Y', 'K':'M', 'Y':'R', 'S':'W', 'M':'K',
					'B':'V', 'D':'H', 'H':'D', 'V':'B',
					'*':'*', 'N':'N', '-':'-'}

def reverse_and_complement(seq):
	rseq = seq[::-1]
	rcseq = ''
	for i in rseq:	# reverse order
		rcseq += complement_dict[i]
	return rcseq



codon_dict = {'TTT':'F', 'TTC':'F', 'TTA':'L', 'TTG':'L',
				'TCT':'S', 'TCC':'S', 'TCA':'S', 'TCG':'S',
				'TAT':'Y', 'TAC':'Y', 'TAA':'*', 'TAG':'*',
				'TGT':'C', 'TGC':'C', 'TGA':'*', 'TGG':'W',
				'CTT':'L', 'CTC':'L', 'CTA':'L', 'CTG':'L',
				'CCT':'P', 'CCC':'P', 'CCA':'P', 'CCG':'P',
				'CAT':'H', 'CAC':'H', 'CAA':'Q', 'CAG':'Q', 
				'CGT':'R', 'CGC':'R', 'CGA':'R', 'CGG':'R',
				'ATT':'I', 'ATC':'I', 'ATA':'I', 'ATG':'M',
				'ACT':'T', 'ACC':'T', 'ACA':'T', 'ACG':'T', 
				'AAT':'N', 'AAC':'N', 'AAA':'K', 'AAG':'K',
				'AGT':'S', 'AGC':'S', 'AGA':'R', 'AGG':'R',
				'GTT':'V', 'GTC':'V', 'GTA':'V', 'GTG':'V',
				'GCT':'A', 'GCC':'A', 'GCA':'A', 'GCG':'A',
				'GAT':'D', 'GAC':'D', 'GAA':'E', 'GAG':'E',
				'GGT':'G', 'GGC':'G', 'GGA':'G', 'GGG':'G',
				'---':'-', 'XXX':'?'}

mixture_regex = re.compile('[WRKYSMBDHVN-]')

mixture_dict = {'W':'AT', 'R':'AG', 'K':'GT', 'Y':'CT', 'S':'CG', 
				'M':'AC', 'V':'AGC', 'H':'ATC', 'D':'ATG', 
				'B':'TGC', 'N':'ATGC', '-':'ATGC'}

#mixture_dict_2 =  [ (set(v), k) for k, v in mixture_dict.iteritems() ]
ambig_dict = dict(("".join(sorted(v)), k) for k, v in mixture_dict.iteritems())


def translate_nuc (seq, offset, resolve=False):
	"""
	Translate nucleotide sequence into amino acid sequence.
		offset by X shifts sequence to the right by X bases
	Synonymous nucleotide mixtures are resolved to the corresponding residue.
	Nonsynonymous nucleotide mixtures are encoded with '?'
	"""
	
	seq = '-'*offset + seq
	
	aa_list = []
	aa_seq = ''	# use to align against reference, for resolving indels
	
	# loop over codon sites in nucleotide sequence
	for codon_site in xrange(0, len(seq), 3):
		codon = seq[codon_site:codon_site+3]
		
		if len(codon) < 3:
			break
		
		# note that we're willing to handle a single missing nucleotide as an ambiguity
		if codon.count('-') > 1 or '?' in codon:
			if codon == '---':	# don't bother to translate incomplete codons
				aa_seq += '-'
			else:
				aa_seq += '?'
			continue
		
		# look for nucleotide mixtures in codon, resolve to alternative codons if found
		num_mixtures = len(mixture_regex.findall(codon))
		
		if num_mixtures == 0:
			aa_seq += codon_dict[codon]
			
		elif num_mixtures == 1:
			resolved_AAs = []
			for pos in range(3):
				if codon[pos] in mixture_dict.keys():
					for r in mixture_dict[codon[pos]]:
						rcodon = codon[0:pos] + r + codon[(pos+1):]
						if codon_dict[rcodon] not in resolved_AAs:
							resolved_AAs.append(codon_dict[rcodon])
			if len(resolved_AAs) > 1:
				if resolve:
					# for purposes of aligning AA sequences
					# it is better to have one of the resolutions
					# than a completely ambiguous '?'
					aa_seq += resolved_AAs[0]
				else:
					aa_seq += '?'
			else:
				aa_seq += resolved_AAs[0]
				
		else:
			aa_seq += '?'
	
	return aa_seq


# =====================================

# when codon sequence contains non-synonymous mixtures (indicated
#	by '?' in translated sequence) then expand into all
#	possible residues in a list-object on which scores are applied
def expand (sd):
	for h in sd.iterkeys():
		aaseq = sd[h]['clipped_aa']
		nucseq = sd[h]['clipped_nuc']
		sd[h].update({'aa_list':expand_single (aaseq, nucseq)}) 
		
	return sd


def expand_single (aaseq, nucseq):
	aa_list = []
	for pos in range(len(aaseq)):
		if aaseq[pos] == '?':
			codon = nucseq[(3*pos):(3*(pos+1))]
			
			# leave in-frame codon gaps alone
			if codon == '---':
				aa_list.append('-')
				continue
			
			if len(codon) < 3:
				print 'WARNING: partial codon "'+codon+'" detected in sequence ' + nucseq + ' at codon ' + str(pos)
				print 'query = ' + aaseq
				print 'h = ' + h
				sys.exit()
			
			rcodons = [codon]
			while 1:
				ok_to_stop = True
				for rcodon in rcodons:
					for pos in range(3):
						if rcodon[pos] in mixture_dict.keys():
							rcodons.remove(rcodon)
							for r in mixture_dict[rcodon[pos]]:
								next_rcodon = rcodon[0:pos] + r + rcodon[(pos+1):]
								if next_rcodon not in rcodons:
									rcodons.append(next_rcodon)
							ok_to_stop = False
							break	# go to next item in list
				if ok_to_stop:
					break
			
			resolved_AAs = []
			for rcodon in rcodons:
				if codon_dict[rcodon] not in resolved_AAs:
					resolved_AAs.append(codon_dict[rcodon])
			
			if codon.count('-') > 0:
				aa_list.append('?')
			else:
				aa_list.append(resolved_AAs)
			"""
			if '-' in codon:
				if len(resolved_AAs) > 1:
					# this will have to be imputed
					# currently '?' contributes no score
					aa_list.append('?')
				else:
					aa_list.append(resolved_AAs[0])
			else:
				aa_list.append(resolved_AAs)
			"""
		else:
			aa_list.append(aaseq[pos])
			
	return aa_list



# =======================================================================
sg_regex = re.compile('[ACGT][N-][ACGT]')

def expand_clonal (sd):
	for h in sd.iterkeys():
		query_v3 = sd[h]['clipped_aa']
		codon_v3 = sd[h]['clipped_nuc']
		
		new_seq = ''
		for pos in range(len(query_v3)):
			if query_v3[pos] == '?':
				codon = codon_v3[(3*pos):(3*(pos+1))]
				
				# leave in-frame codon gaps alone
				if codon == '---':
					new_seq += '-'
					continue
				
				if len(codon) < 3:
					print 'WARNING: partial codon "'+codon+'" detected in sequence ' + codon_v3 + ' at codon ' + str(pos)
					print 'query = ' + query_v3
					print 'h = ' + h
					sys.exit()
				
				rcodons = [codon]
				while 1:
					ok_to_stop = True
					for rcodon in rcodons:
						for pos in range(3):
							if rcodon[pos] in mixture_dict.keys():
								rcodons.remove(rcodon)
								for r in mixture_dict[rcodon[pos]]:
									next_rcodon = rcodon[0:pos] + r + rcodon[(pos+1):]
									if next_rcodon not in rcodons:
										rcodons.append(next_rcodon)
								ok_to_stop = False
								break	# go to next item in list
					if ok_to_stop:
						break
				
				resolved_AAs = []
				for rcodon in rcodons:
					if codon_dict[rcodon] not in resolved_AAs:
						resolved_AAs.append(codon_dict[rcodon])
				
				if len(resolved_AAs) > 1:
					new_seq += '?'
				else:
					new_seq += resolved_AAs[0]
			else:
				new_seq += query_v3[pos]
		
		sd[h].update({'clipped_aa':new_seq})
	
	return sd

def patch_gaps (sd):
	# For clonal sequences (454), singleton 'N's or '-'s are common.
	# Rather than ignore these incomplete codons, it is better to resolve
	# them by one of the following procedures:
	#	- if the gap is at a synonymous site, then simply replace the 
	#		ambiguous character '?' with that residue
	#	- if the gap is at a nonsynonymous site, then resolve it into the
	#		majority nucleotide and the corresponding residue
	
	#  I'm not sure this is the best approach...
	
	# generate nucleotide frequency vector
	seqlen = len(sd.values()[0]['clipped_nuc'])
	nucfreqs = dict([(x,{'A':0, 'C':0, 'G':0, 'T':0}) for x in range(seqlen)])
	for h in sd.iterkeys():
		nucseq = sd[h]['clipped_nuc']
		for pos in range(seqlen):
			try:
				nucfreqs[pos][nucseq[pos]] += 1
			except:
				pass
	
	# generate majority consensus sequence
	major_seq = ''
	for pos in range(seqlen):
		major_seq += nucfreqs[pos].keys()[nucfreqs[pos].values().index(max(nucfreqs[pos].values()))]
	
	for h in sd.iterkeys():
		nucseq = sd[h]['clipped_nuc']
		
		if not sg_regex.findall(nucseq):
			continue
		
		nslist = [x for x in nucseq]
		
		# edit a.a. sequence based on resolution of broken codon
		bad_seq = False
		for sg in sg_regex.finditer(nucseq):
			nucpos = sg.start()+1
			try:
				nslist[nucpos] = major_seq[nucpos]
			except:
				# this is usually caused by a frameshift in the
				# original sequence that is not handled properly by
				# the align() function.
				bad_seq = True
				break
		
		if bad_seq:
			continue
		
		nucseq = ''.join(nslist)
		sd[h]['clipped_nuc'] = nucseq
		sd[h]['clipped_aa'] = translate_nuc(nucseq, 0)
	
	return sd


# =======================================================================
def aalist_to_str (aa_list):
	res = ''
	for pos in range(len(aa_list)):
		if type(aa_list[pos]) == list:
			res += '['
			for char in aa_list[pos]:
				res += char
			res += ']'
		else:
			res += aa_list[pos]
	return res


# =======================================================================

def gaps2ambig (nucseq):
	"""
	Convert all gap characters that are not a proper codon deletion
	into an ambiguous character.  Meant to operate on a nucleotide
	sequence that is in reading frame.
	"""
	aaseq = translate_nuc(nucseq, 0)

	

def consensus(column, alphabet='ACGT', resolve=False):
	"""
	Plurality consensus - nucleotide with highest frequency.
	In case of tie, report mixtures.
	"""
	freqs = {}
	for char in alphabet:
		freqs.update({char: 0})
	#freqs = {"A": 0, "T": 0, "C": 0, "G": 0, "-": 0}
	for char in column:
		if char in alphabet:
			freqs[char] += 1
		elif mixture_dict.has_key(char):
			# handled ambiguous nucleotides with equal weighting
			resolutions = mixture_dict[char]
			for char2 in resolutions:
				freqs[char2] += 1./len(resolutions)
		else:
			# unrecognized nucleotide character
			pass
			
	base = max(freqs, key=lambda n: freqs[n])
	max_count = freqs[base]
	possib = filter(lambda n: freqs[n] == max_count, freqs)
	if len(possib) == 1:
		return possib[0]
	elif "-" in possib:
		if resolve:
			possib.remove("-")
			if len(possib) == 0:
				return "-"
			elif len(possib) == 1:
				return possib[0]
			else:
				return ambig_dict["".join(sorted(possib))]
		else:
			# gap character overrides ties
			return "-"
	else:
		return ambig_dict["".join(sorted(possib))]


def majority_consensus (fasta, threshold = 0.5, alphabet='ACGT', ambig_char = 'N'):
	"""
	Return majority-rule consensus sequence.
	[threshold] = percentage of column that most common character must exceed
	[alphabet] = recognized character states
	"""
	
	"""
	res = ''
	if len(alphabet) == 0: alphabet = set(fasta[0][1])
	columns = transpose_fasta(fasta)
	for col in columns:
		cset = set(col)
		if len(cset) == 1:
			c = cset.pop()
			if c not in alphabet: res += ambig_char
			else: res += c
		else:
			counts = [(col.count(c), c) for c in cset if c in alphabet]
			if len(counts) == 0:
				res += ambig_char
				continue
			counts.sort(reverse=True) # descending order
			max_count, max_char = counts[0]
			if max_count / float(len(fasta)) > threshold: res += max_char
			else: res += ambig_char
	return res
	"""
	
	consen = []
	columns = transpose_fasta(fasta)
	seqs = [s for h, s in fasta]
	
	for column in columns:
		consen.append(consensus(column, alphabet=alphabet, resolve=False))
	
	newseq = "".join(consen)
	
	"""
	# Resolve missing data.
	# Proper indels start and end in-frame.
	indel_ptn = re.compile("(.{3})*?(?P<indel>(\?{3})+)")
	indels = []
	for match in indel_ptn.finditer(newseq):
		indels.extend(range(*match.span("indel")))
	
	for column in range(len(consen)):
		if consen[column] == "?" and column not in indels:
			consen[column] = consensus(column, resolve=True)
	
	return "".join(consen)
	"""
	return newseq
	

# =======================================================================
"""
transpose_fasta - return an array of alignment columns
"""
def transpose_fasta (fasta):
	# some checks to make sure the right kind of object is being sent
	if type(fasta) is not list:
		return None
	if type(fasta[0]) is not list or len(fasta[0]) != 2:
		return None
	
	n_columns = len(fasta[0][1])
	res = []
	for c in range(n_columns):
		res.append ( [ s[c] for h, s in fasta ] )
	
	return res

def untranspose_fasta(tfasta):
	nseq = len(tfasta[0])
	res = [ '' for s in range(nseq) ]
	for col in tfasta:
		for i in range(nseq):
			res[i] += col[i]
	return res



"""
entropy_from_fasta
	Calculate the mean entropy over columns of an alignment
	passed as a FASTA object (list of lists).
	Defaults to the nucleotide alphabet.
	If a vector of counts is passed, then entropy calculations will
	be weighted by the frequency of each sequence.  Otherwise each
	sequence will be counted as one instance.
	
	NOTE: Non-alphabet characters (e.g., mixtures, gaps) are being simply ignored!
	
	Test code:
	infile = open('/Users/apoon/wip/etoi/screened/ACT 60690_NEF_forward.pre2.screen1', 'rU')
	fasta = convert_fasta(infile.readlines())
	infile.close()
	counts = [int(h.split('_')[1]) for h, s in fasta]
	
"""
def entropy_from_fasta (fasta, alphabet = 'ACGT', counts = None):
	columns = transpose_fasta (fasta)
	ents = []
	for col in columns:
		ent = 0.
		
		# expand character count in vector if 'counts' argument is given
		if counts:
			new_col = []
			for i in range(len(col)):
				new_col.extend( [ col[i] for j in range(counts[i]) ] )
			col = new_col
		
		for char in alphabet:
			freq = float(col.count(char)) / len(col)
			if freq > 0:
				ent -= freq * math.log(freq, 2)
		
		ents.append(ent)
	
	mean_ent = sum(ents) / len(ents)
	return mean_ent



def bootstrap(fasta):
	"""
	Random sampling of columns with replacement from alignment.
	Returns a FASTA (list of lists)
	"""
	nsites = len(fasta[0][1])
	seqnames = [h for (h, s) in fasta]
	res = []
	tfasta = transpose_fasta(fasta)
	sample = []
	for j in range(nsites):
		sample.append(tfasta[random.randint(0, nsites-1)])
	
	seqs = untranspose_fasta(sample)
	for k, s in enumerate(seqs):
		res.append([seqnames[k], s])
	
	return res

