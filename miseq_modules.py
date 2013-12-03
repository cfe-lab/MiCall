def csf2counts (path,mode,mixture_cutoffs):
	"""
	Calculate nucleotide and amino acid counts from a FASTA or CSF file
	"""

	import csv,HyPhy,os,sys
	from hyphyAlign import change_settings, get_boundaries, pair_align
	from miseqUtils import ambig_dict, convert_csf, convert_fasta, mixture_dict, translate_nuc

	# FIXME: PARAMETERIZE THIS TO BE CONTROLLED FROM 1_MPI_WRAPPER
	amino_reference_sequence = "csf2counts_amino_sequences.csv"
	
	hyphy = HyPhy._THyPhy (os.getcwd(), 1)
	change_settings(hyphy)
	amino_alphabet = 'ACDEFGHIKLMNPQRSTVWY*'
	
	filename = path.split('/')[-1]
	sample, ref = filename.split('.')[:2]
	
	# Load HXB2 amino reference sequences (FIXME: FIX DOCUMENTATION TO INCLUDE HCV1A-H77)
	with open(amino_reference_sequence, "rb") as f:
		input_file = csv.reader(f)
		hxb2 = {}
		for row in input_file:
			region, amino = row
			hxb2[region] = amino
	
	# Make the output stem by removing the extension of the filename
	root = '/'.join(path.split('/')[:-1])
	if root == '':
		# in case script is executed on file in cwd
		root = '.'
	
	outpath = root + '/' + (filename.replace('.fasta', '') if mode == 'Amplicon' else filename.replace('.csf', ''))
	
	# output to files and compute consensus
	nucfile = open(outpath+'.nuc.csv', 'w')
	nucfile.write('query.nuc.pos,hxb2.nuc.pos,A,C,G,T\n')
	
	confile = open(outpath+'.conseq', 'w')
	indelfile = open(outpath+'.indels.csv', 'w')
	
	if not hxb2.has_key(ref):
		sys.exit()
	
	refseq = hxb2[ref]
	
	with open(path, 'rU') as infile:
		if mode == 'Nextera':
			fasta, lefts, rights = convert_csf(infile.readlines())
		elif mode == 'Amplicon':
			fasta = convert_fasta(infile.readlines())
		else:
			sys.exit()
	
	
	# Use the first read to determine reading frame
	max_score = 0
	best_frame = 0
	for frame in range(3):
		prefix = ('-'*lefts[fasta[0][0]] if mode == 'Nextera' else '')
		p = translate_nuc(prefix + fasta[0][1], frame)
		aquery, aref, ascore = pair_align(hyphy, refseq, p)
		if ascore > max_score:
			best_frame = frame # the reading frame of left = 0
			max_score = ascore
	
	# Iterate through reads and count WHAT?
	nucs = {}
	aminos = {}
	pcache = []	# Cache protein sequences
	
	# At this point, sequences are aligned against the sample-region
	# specific consensus. Thus, each read in the csf contains an offset
	# with respect to the sample-region specific consensus.
	
	# For each sequence in the fasta/csf
	for i, (h, s) in enumerate(fasta):
	
		# If this is Fasta, there is no offset
		left = lefts[h] if mode == 'Nextera' else 0
	
		# Headers contain read count in last entry
		count = 1 if mode == 'Nextera' else int(h.split('_')[-1])
	
		# Update nucleotide counts (with respect to self-coordinate?)
		for j, nuc in enumerate(s):
			pos = left + j
			if not nucs.has_key(pos):
				nucs.update({pos: {}})
			if not nucs[pos].has_key(nuc):
				nucs[pos].update({nuc: 0})
			nucs[pos][nuc] += count
		
		p = translate_nuc('-'*left + s, best_frame)
		pcache.append(p)
	
		# Update amino counts
		for pos, aa in enumerate(p):
			if aa == '-':
				continue
			if not aminos.has_key(pos):
				aminos.update({pos: {}})
			if not aminos[pos].has_key(aa):
				aminos[pos].update({aa: 0})
			aminos[pos][aa] += count
	
	# Generate AA plurality (max) consensus
	keys = aminos.keys()
	keys.sort()
	aa_max = ''
	
	for pos in keys:
		intermed = [(v, k) for k, v in aminos[pos].iteritems()]
		intermed.sort(reverse=True)
		aa_max += intermed[0][1]
	
	# Align consensus against HXB2
	aquery, aref, ascore = pair_align(hyphy, refseq, aa_max)
	
	# Ignore parts of query outside our reference
	left, right = get_boundaries(aref)
	qindex_to_hxb2 = {} # Maps query amino to HXB2 amino coordinates
	inserts = []		# Leep track of which aa positions are insertions
	qindex = 0
	rindex = 0
	for i in range(len(aref)):
		# ignore parts of query that do not overlap reference
		if i < left:
			qindex += 1
			continue
		if i >= right:
			break
		
		if aref[i] == '-':
			# insertion in query
			inserts.append(qindex)
			qindex += 1
		elif aquery[i] == '-':
			# deletion in query
			# do not increment qindex
			rindex += 1
			continue
		else:
			qindex_to_hxb2.update({qindex: rindex})
			qindex += 1
			rindex += 1
	
	# Reiterate through sequences to capture indels
	if len(inserts) > 0:
		indelfile.write('insertion,count\n')
		indel_counts = {}
		
		for p in pcache:
			ins_str = str(inserts[0])
			last_i = -1
			for i in inserts:
				if last_i > -1 and i - last_i > 1:
					# end of a contiguous indel
					ins_str += ',%d' % i
				try:
					ins_str += p[i]
				except IndexError:
					break
				last_i = i
			
			if not indel_counts.has_key(ins_str):
				indel_counts.update({ins_str: 0})
			
			indel_counts[ins_str] += 1
		
		for ins_str, count in indel_counts.iteritems():
			indelfile.write('%s,%d\n' % (ins_str, count))
	
	indelfile.close()
	
	
	# output nucleotide and amino acid counts in HXB2 coordinates
	# also output consensus sequences at varying thresholds
	
	keys = nucs.keys()	# nucs[self-pos][nuc] = count
	keys.sort()
	maxcon = ''
	conseqs = ['' for cut in mixture_cutoffs]
	
	codon_pos = 0
	
	# For each base coordinate in the query
	for pos in keys:
		try:
			aapos = pos/3				# Get the amino position
			codon_pos = pos % 3
			hxb2_pos = qindex_to_hxb2[aapos] + 1	# Get the HXB2 based on that amino position
		except KeyError:
			continue
	
		# hxb2_pos is the hxb2 amino position, pos is the hxb2 nucleotide position
		hxb2_nuc_pos = 3*hxb2_pos + codon_pos
		nucfile.write('%d,%d,%s\n' % (pos, hxb2_nuc_pos, ','.join(map(str, [nucs[pos].get(nuc, 0) for nuc in 'ACGT']))))
		
		# plurality consensus
		intermed = [(count, nuc) for nuc, count in nucs[pos].iteritems()]
		intermed.sort(reverse=True)
		maxcon += intermed[0][1]
		
		# consensuses with mixtures determined by frequency cutoff
		total_count = sum([count for count, nuc in intermed])
		
		for ci, cutoff in enumerate(mixture_cutoffs):
			mixture = []
			for count, nuc in intermed:
				if float(count) / total_count > cutoff:
					mixture.append(nuc)
			
			if 'N' in mixture:
				if len(mixture) > 1:
					mixture.remove('N')
				else:
					# completely ambiguous
					conseqs[ci] += 'N'
					continue
			
			if '-' in mixture:
				if len(mixture) > 1:
					mixture.remove('-')
				else:
					conseqs[ci] += '-'
					continue
			
			# encode nucleotide mixture
			if len(mixture) > 1:
				mixture.sort()
				conseqs[ci] += ambig_dict[''.join(mixture)]
			elif len(mixture) == 1:
				conseqs[ci] += mixture[0]
			else:
				# mixture of length zero, no bases exceed cutoff
				conseqs[ci] += 'N'
	nucfile.close()
	
	# output consensus sequences
	confile.write('>%s_MAX\n%s\n' % (sample, maxcon))
	for ci, cutoff in enumerate(mixture_cutoffs):
		confile.write('>%s_%1.3f\n%s\n' % (sample, cutoff, conseqs[ci]))
	confile.close()
	
	# Output amino acid counts
	with open(outpath+'.amino.csv', 'w') as aafile:
		aafile.write('query.aa.pos,hxb2.aa.pos,%s\n' % ','.join(list(amino_alphabet)))
		keys = aminos.keys()
		keys.sort()
		for aapos in keys:
			if aapos in inserts:
				continue
			try:
				hxb2_pos = qindex_to_hxb2[aapos] + 1
			except KeyError:
				continue
			aafile.write('%d,%d,%s\n' % (aapos,hxb2_pos, ','.join(map(str, [aminos[aapos].get(aa, 0) for aa in amino_alphabet]))))

def system_call(command):
	import logging, os
	logger = logging.getLogger()
	logger.debug(command)
	os.system(command)

def remap (R1_fastq, R2_fastq, samfile, ref, original_reference, conseq_qCutoff=30):
	""" 
	1) Generate sample-specific consensus from a samtools pileup.
	2) Remap everything to this consensus as a ref seq
	3) Returns paths to the SAM output and consensus sequence file
	"""
	import logging
	from miseq_modules import system_call

	logger = logging.getLogger()

	bamfile = samfile.replace('.sam', '.bam')
	confile = "{}.pileup.conseq".format(bamfile)

	# Overwrite previous iterations of the remapped sam (Only keep the original specific prelim sam)
	remapped_sam = samfile if samfile.endswith('.remap.sam') else samfile.replace('.sam', '.remap.sam')

	# If this is the first run, use the static reference
	if ref == original_reference:
		system_call('samtools view -bt {}.fasta.fai {} > {} 2>/dev/null'.format(ref, samfile, bamfile))
	else:
		# Make new samtools index file (Creates [ref].fai)
		system_call('samtools faidx {}'.format(ref))
		system_call('samtools view -bt {}.fai {} > {} 2>/dev/null'.format(ref, samfile, bamfile))

	# Sort the bam file by leftmost position on the reference assembly
	system_call('samtools sort {} {}.sorted'.format(bamfile, bamfile))

	# Make a pileup from the sorted bam
	system_call('samtools mpileup -A {}.sorted.bam > {}.pileup 2>/dev/null'.format(bamfile, bamfile))
	
	# Create new consensus sequence from the pileup
	system_call('python2.7 pileup2conseq_v2.py {}.pileup {}'.format(bamfile, conseq_qCutoff))
	
	# Convert consensus into bowtie2 reference files (Creates 6 files of *.bt2)
	system_call('bowtie2-build -f -q {} {}'.format(confile, confile))
	
	# Map original fastq reads to new reference
	cmd = 'bowtie2 --quiet -p 1 --local -x {} -1 {} -2 {} -S {} --no-unal --met-file {} --un {} --un-conc {}'.format(confile,
			R1_fastq, R2_fastq, remapped_sam, remapped_sam.replace('.sam', '.bt2_metrics'),
			remapped_sam.replace('.sam', '.bt2_unpaired_noalign.fastq'), 
			remapped_sam.replace('.sam', '.bt2_paired_noalign.fastq'))
	system_call(cmd)
	
	return remapped_sam, confile


def mapping(refpath, R1_fastq, conseq_qCutoff, mode, is_t_primer, REMAP_THRESHOLD, MAX_REMAPS):
	"""
	refpath		Absolute path to static reference sequences used during original mapping
	R1_fastq	Absolute path to R1 fastq file
	conseq_qCutoff	Q-cutoff for determining if a base contributes to the consensus in pileup2conseq
	is_t_primer	Used to detect contaminants from previous runs (FIXME: REMOVE?)
	REMAP_THRESHOLD	Efficiency of read mapping at which we stop bothering to perform additional remapping
	MAX_REMAPS	Number of extra attempts at remapping before giving up
	"""
	import logging, os, subprocess, sys
	from miseqUtils import samBitFlag
	from miseq_modules import system_call

	logger = logging.getLogger()
	original_reference = refpath		# Path to the original reference sequence

	# Store region codes of static reference fasta in refnames (ConB for HIV, H77 for HCV)
	with open(refpath+'.fasta', 'rU') as ref_fasta:
		refnames = []
		for line in ref_fasta:
			if line.startswith('>'): refnames.append(line.strip('>\n'))

	# Deduce R1/R2 file pairing
	root = '/'.join(R1_fastq.split('/')[:-1])
	filename = os.path.basename(R1_fastq)			# Filename of R1 fastq...
	file_prefix = filename.split('.')[0]			# Has a prefix containing...
	sample_name, sample_well = file_prefix.split('_')[:2]	# Sample name and well
	prefix = '_'.join([sample_name, sample_well])		# Change this prefix to be joined by _ (??)
	count_file = open(R1_fastq.replace('.fastq', '.counts'), 'w')
	R2_fastq = R1_fastq.replace('R1', 'R2')

	# Determine the number of reads in both (R1 + R2) fastq files, store in the .count file
	stdout, stderr = subprocess.Popen(['wc', '-l', R1_fastq], stdout = subprocess.PIPE, stderr = subprocess.PIPE).communicate()
	total_reads_R1 = int(stdout.split()[0])/4
	stdout, stderr = subprocess.Popen(['wc', '-l', R2_fastq], stdout = subprocess.PIPE, stderr = subprocess.PIPE).communicate()
	total_reads_R2 = int(stdout.split()[0])/4
	logger.debug("{} R1 and {} R2 reads in {} and {}".format(total_reads_R1, total_reads_R2, R1_fastq, R2_fastq))
	count_file.write('Raw FASTQ,R1,%d,R2,%d\n' % (total_reads_R1, total_reads_R2))

	# Initial consensus B mapping
	prelim_samfile = '{}/{}.prelim.sam'.format(root, prefix)
	system_call('bowtie2 --quiet -p 1 --local -x %s -1 %s -2 %s -S %s' % (refpath, R1_fastq, R2_fastq, prelim_samfile))

	# Define region-specific SAMs: refsams[refname] points to a nested dict which includes a file handle to each specific SAM
	refsams = {}
	for i, refname in enumerate(refnames):
		region_specific_sam = '%s/%s.%s.sam' % (root, prefix, refname)
		refsams.update({refname: {'sam_file_handle': open(region_specific_sam, 'w'),'count': [0,0]}})
	refsams.update({'*': {'sam_file_handle': open('%s/%s.unmapped.sam' % (root, prefix), 'w'),'count': [0,0]}})

	# Subdivide prelim SAMs into region-specific SAMs
	prelim_sam_infile = open(prelim_samfile, 'rU')
	line_counts = [0, 0]
	t_counts = [0, 0]
	contam_file = open(R1_fastq.replace('.fastq', '.Tcontaminants.fastq'), 'w')

	for line in prelim_sam_infile:

		# Copy the original SAM header into each region-specific SAM
		if line.startswith('@'):
			for refname in refnames: refsams[refname]['sam_file_handle'].write(line)
			continue

		# SAM documentation explains what these fields mean - http://samtools.sourceforge.net/SAMv1.pdf
		qname, flag, refname, pos, mapq, cigar, rnext, pnext, tlen, seq, qual = line.strip('\n').split('\t')[:11]

		# SAM bitwise flag variable specifies whether read is paired, successfully mapped, etc
		bitinfo = samBitFlag(flag)

		# Ignore unmapped reads
		if bitinfo['is_unmapped']:
			continue

		if mode == 'Amplicon':
			# filter T contaminants
			if is_t_primer:
				# Trim T(A) when known to be present due to primer
				if bitinfo['is_reverse'] and seq.endswith('A'):
					seq = seq[:-1]
				elif not bitinfo['is_reverse'] and seq.startswith('T'):
					seq = seq[1:]
				else:
					t_counts[bitinfo['is_reverse']] += 1 # missing A/T
					# Dump in fastq format to contaminant file
					contam_file.write('@%s\n%s\n+\n%s\n' % (qname, seq, qual))
					continue

			else:
				# look for contaminating T primed reads in non-T run
				if bitinfo['is_reverse'] and seq.endswith('A') and cigar.endswith('1S'):
					t_counts[bitinfo['is_reverse']] += 1 # unexpected A
					contam_file.write('@%s\n%s\n+\n%s\n' % (qname, seq, qual))
					continue
				elif not bitinfo['is_reverse'] and seq.startswith('T') and cigar.startswith('1S'):
					t_counts[bitinfo['is_reverse']] += 1 # unexpected T
					contam_file.write('@%s\n%s\n+\n%s\n' % (qname, seq, qual))
					continue
				else:
					pass
		
		# Write each line into respective region-specific SAMs
		try:
			refsams[refname]['sam_file_handle'].write(line)
			refsams[refname]['count'][bitinfo['is_reverse']] += 1
		except:
			print "{} appeared to be an incorrect refname".format(refname)

		# FIXME: Does this identify R1 reads from R2 in the prelim SAM...? 
		line_counts[bitinfo['is_reverse']] += 1

	prelim_sam_infile.close()
	contam_file.close()

	# Show the number of reads that made it to preliminary mapping in the count file
	count_file.write('Preliminary map,{},{}\n'.format(line_counts[0], line_counts[1]))
	if is_t_primer: error_type = "Missing T"
	else: error_type = "Unexpected T"
	count_file.write('{} (Detected at prelim map),R1,{},R2,{}\n'.format(error_type, t_counts[0], t_counts[1]))

	# Show the reads mapped at preliminary mapping with respect to region
	for refname in refnames:
		refsams[refname]['sam_file_handle'].close()
		count_file.write('prelim %s,R1,%d,R2,%d\n' % (refname, refsams[refname]['count'][0],refsams[refname]['count'][1]))

	# Track the number of reads that have been successfully remapped in total
	total_remap = 0

	# Remap the fastqs using sample/region specific conseqs
	for refname in refnames:

		# Ignore phiX, unmapped reads, and regions which had no mapping at the prelim mapping stage
		if sum(refsams[refname]['count']) == 0 or refname == 'phiX174' or refname == '*':
			continue

		# Run remap on the region-specific sam, and get the remapped sam and consensus pileup used to generate it
		samfile = refsams[refname]['sam_file_handle'].name

		logging.info("remap({},{},{},{},{},{})".format(R1_fastq, R2_fastq, samfile, refpath, original_reference, conseq_qCutoff))
		samfile, confile = remap(R1_fastq, R2_fastq, samfile, refpath, original_reference, conseq_qCutoff)
	
		# Track file paths
		refsams[refname].update({'samfile': samfile, 'confile': confile})
	
		# Track the number of mapped reads in each region
		stdout, stderr = subprocess.Popen(['wc', '-l', samfile], stdout = subprocess.PIPE, stderr = subprocess.PIPE).communicate()
		region_specific_count = (int(stdout.split()[0]) - 3) / 2. # First 3 lines contain comments
		total_remap += region_specific_count
		refsams[refname]['count'][0] = region_specific_count
		count_file.write('remap %s,%d\n' % (refname, int(region_specific_count)))

	# Continue to remap if we've failed to map enough total reads
	mapping_efficiency = total_remap / total_reads_R1
	if mapping_efficiency < REMAP_THRESHOLD:
		break_out = False
		logger.info("Poor mapping for {} (Mapping efficiency: {}) - trying additional remapping".format(sample_name, mapping_efficiency))

		# Repeat remapping up to the number of MAX_REMAPS permitted
		for iter in range(MAX_REMAPS):
			logger.debug("Additional remapping for {} (iteration #{})".format(sample_name, iter))

			total_remap = 0
			for refname in refnames:
				if refsams[refname]['count'][0] == 0 or refname == 'phiX174' or refname == '*':
					continue

				samfile = refsams[refname]['samfile']
				confile = refsams[refname]['confile']
				samfile, confile = remap(R1_fastq, R2_fastq, samfile, confile, original_reference, conseq_qCutoff)
				refsams[refname]['samfile'] = samfile
				refsams[refname]['confile'] = confile
			
				# Continue to determine the number of reads mapped in this region
				stdout, stderr = subprocess.Popen(['wc', '-l', samfile], stdout = subprocess.PIPE, stderr = subprocess.PIPE).communicate()
				region_specific_count = (int(stdout.split()[0]) - 3) / 2.
				
				if region_specific_count < refsams[refname]['count'][0]:
					logger.warn("Remapping for {} resulting in LESS reads - halting iterative remapping".format(sample_name))
					break_out = True
					break
			
				total_remap += region_specific_count
				refsams[refname]['count'][0] = region_specific_count

				# Write each remap iteration result to the file
				count_file.write('remap %d %s,%d\n' % (iter, refname, int(region_specific_count)))

			if break_out or total_remap / total_reads_R1 >= REMAP_THRESHOLD:
				break

	count_file.close()
	
def g2p_scoring(input_path, g2p_alignment_cutoff):
	"""
	Generate V3-specific nucleotide sequence from remapped env .fasta file,
	along with G2PFPR score in the header (Plus the count)
	"""

	import logging,os,sys
	from hyphyAlign import apply2nuc, change_settings, get_boundaries, HyPhy, pair_align, refSeqs
	from minG2P import conan_g2p
	from miseqUtils import convert_fasta, timestamp, translate_nuc

	logger = logging.getLogger()
	hyphy = HyPhy._THyPhy (os.getcwd(), 1)			# HyPhy is used for alignment
	change_settings(hyphy)					# Configure scoring matrix / gap penalties
	refseq = translate_nuc(refSeqs['V3, clinical'], 0)	# The V3 reference sequence is NON-STANDARD: talk to Guin

	if not input_path.endswith('.fasta'):
		logging.error("{} doesn't end in .fasta")
		sys.exit()
	
	with open(input_path, 'rU') as infile:
		try:
			fasta = convert_fasta(infile.readlines())
		except:
			logging.error('g2p_scoring(): {} not a valid fasta file'.format(input_path))
			sys.exit()
	
	# Determine offset from 1st sequence to correct frameshift induced by sample-specific remapping
	seq1 = fasta[0][1].strip("-")
	best_offset = 0
	best_score = -999
	for offset in range(3):
		aaEnvSeq = translate_nuc(seq1, offset)
		aquery, aref, ascore = pair_align(hyphy, refseq, aaEnvSeq)
		if ascore > best_score:
			best_offset = offset
			best_score = ascore

	# For each env sequence, extract the V3 nucleotide sequence
	badfile = open(input_path.replace('.fasta', '.badV3'), 'w')
	v3nucs = {}
	for header, seq in fasta:
		count = int(header.split('_')[-1])
		seq = seq.replace("-","")					# Strip dashes at flanking regions generated by alignment
		aaEnvSeq = translate_nuc(seq, best_offset)			# Translate env on correct ORF
		aquery, aref, ascore = pair_align(hyphy, refseq, aaEnvSeq)
		left, right = get_boundaries(aref)				# Get left/right boundaries of V3 protein
		v3prot = aquery[left:right]					# Extract V3 protein
		v3nuc = apply2nuc(seq[(3*left-best_offset):], v3prot,		# Use alignment to extract V3 nuc seq
				aref[left:right], keepIns=True, keepDel=False)
		
		# Drop V3 data that don't satisfy quality control
		if 'N' in v3nuc or not v3prot.startswith('C') or not v3prot.endswith('C') or '*' in v3prot or ascore < g2p_alignment_cutoff or len(v3prot) < 32 or len(v3prot) > 40:
			badfile.write('>%s_reason_%s\n%s\n' % (header,
				'|'.join(['stopcodon' if '*' in v3prot else '',					# V3 can't have internal stop codon
				'lowscore' if ascore < g2p_alignment_cutoff else '',				# The G2P alignment can't be poor
				'cystines' if not v3prot.startswith('C') or not v3prot.endswith('C') else '',	# V3 must start/end with C
				'ambig' if 'N' in v3nuc else '']),seq))						# There must be no unknown bases
		else:
			# Track the count of each v3 nucleotide sequence
			if v3nucs.has_key(v3nuc):
				v3nucs[v3nuc] += count
			else:
				v3nucs.update({v3nuc: count})
	badfile.close()
	
	# Calculate g2p scores for each v3 nuc sequence
	v3prots = {}
	for v3nuc, count in v3nucs.iteritems():
		g2p, fpr, aligned = conan_g2p(v3nuc)
	
		if g2p is None:
			continue
	
		# Track the count of each protein sequence
		if v3prots.has_key(aligned):
			v3prots[aligned]['count'] += count
		else:
			# Dict within dict - store count and fpr for each sequence
			v3prots.update({aligned: {'count': count, 'fpr': fpr}})
	
	# Collect identical V3 amino acid sequences and output
	# Extract the protein sequence and it's count
	# k: protein sequence, v: dict mapping to 'count' and 'fpr'
	intermed = [(v['count'], k) for k, v in v3prots.iteritems()]
	intermed.sort(reverse=True)
	
	# Write a file with the sample (prefix?), rank, count, fpr, and sequence
	v3prot_path = input_path.replace('.fasta', '.v3prot')
	with open(v3prot_path, 'w') as v3protfile:
		for i, (count, v3prot) in enumerate(intermed):
			fpr = v3prots[v3prot]['fpr']
			v3protfile.write('>%s_variant_%d_count_%d_fpr_%s\n%s\n' % (prefix, i, count, fpr, v3prot))

def sam2fasta_with_base_censoring(samfile, censoring_qCutoff, mapping_cutoff, mode, max_prop_N):
	"""
	From a SAM file, create a FASTA file for Amplicon runs, or CSF
	files for Nextera runs. Both output formats can be thought of as
	simplified representations of the original SAM.
	"""
	import logging, os, sys
	from miseqUtils import len_gap_prefix, sam2fasta

	logger = logging.getLogger()

	# Extract sample (prefix) and region
	filename = samfile.split('/')[-1]
	prefix, region = filename.split('.')[:2]

	# Convert SAM to fasta-structured variable
	with open(samfile, 'rU') as infile:
		logging.debug("sam2fasta({}, {}, {}, {})".format(infile, censoring_qCutoff, mapping_cutoff, max_prop_N))
		fasta = sam2fasta(infile, censoring_qCutoff, mapping_cutoff, max_prop_N)

	# Send warning to standard out if sam2fasta didn't return anything
	if fasta == None:
		logging.warn("{} likely empty or invalid - halting sam2fasta".format(samfile))
		sys.exit()

	# For Amplicon runs, generate a (compressed) fasta file
	if mode == 'Amplicon':

		# Store identical sequences as a single FASTA entry with count data in the header
		d = {}
		for h, s in fasta:
			if d.has_key(s):
				d[s] += 1
			else:
				d.update({s: 1})

		# Sort the fasta by read count and write the fasta to disk
		intermed = [(count, s) for s, count in d.iteritems()]
		intermed.sort(reverse=True)

		fasta_filename = '.'.join(map(str,[samfile.replace('.remap.sam', ''), censoring_qCutoff, 'fasta']))
		with open(fasta_filename, 'w') as outfile:
			for i, (count, seq) in enumerate(intermed):
				outfile.write('>%s_variant_%d_count_%d\n%s\n' % (prefix, i, count, seq))

	# For Nextera runs, write a csf file (Our proprietary format)
	elif mode == 'Nextera':

		# Sort csf by left-gap prefix: the offset of the read relative to the ref seq
		intermed = [(len_gap_prefix(s), h, s) for h, s in fasta]
		intermed.sort()
		csf_filename = '.'.join(map(str,[samfile.replace('.remap.sam', ''), censoring_qCutoff, 'csf']))
		with open(csf_filename, 'w') as outfile:
			for (gp, h, seq) in intermed:
				outfile.write('%s,%d,%s\n' % (h, gp, seq.strip('-')))

def slice_outputs(root, region_slices):
	"""
	Slice matrix output of previous step (*.nuc|amino.count.csv) into sub-regions.
	"""

	import os, sys, time
	from glob import glob
	
	# Coordinates are in nucleotide space: start/end are inclusive (Relative to HXB2 aligned sequences)
	# Ex: region_slices = [("PROTEASE", "HIV1B-pol", 1, 297), ("V3", "HIV1B-env", 887, 993)]
	
	# For each region slice rule
	for rule in region_slices:
		slice, region, start, end = rule
		files = glob(root + '/*.{}.*nuc.csv'.format(region))
		files += glob(root + '/*.{}.*amino.csv'.format(region))
	
		# Get all nuc/amino files containing the region to be sliced
		for path in files:
			fileName = os.path.basename(path)
			sample,old_region = fileName.split(".")[:2]
	
			with open(path, 'rU') as f:
				lines = f.readlines()
	
			newFileName = fileName.replace(region,slice)
			dirName = os.path.dirname(path)
			slice_filename = "{}/{}".format(dirName, newFileName)
			f = open(slice_filename, 'w')
			conseq_filename = "{}/{}".format(dirName, newFileName.replace(".csv",".conseq"))
			f_conseq = open(conseq_filename, 'w')
	
			conseq = ""
			is_empty = True
	
			# For each line in the frequency matrix file
			for i,line in enumerate(lines):
				line = line.rstrip("\n")
	
				# First, extract character dictionary from header
				if (i == 0):
					f.write("{}\n".format(line))
					dictionary = line.split(",")[2:]
					continue
	
				# For nuc.csv, hxb2_pos is in nucleotide space
				query_pos, hxb2_pos = map(int, line.split(",")[:2])
				if "nuc.csv" in path:
					if hxb2_pos < start or hxb2_pos > end+1: continue
					region_pos = hxb2_pos - start + 1
	
				# For amino.csv, hxb2_pos is in amino space
				elif "amino.csv" in path:
					if hxb2_pos < (start+2)/3 or hxb2_pos > end/3: continue
					region_pos = hxb2_pos - (start+2)/3 + 1
	
				# If we reached this point, the slice contains data
				is_empty = False
	
				# Generate consensus sequence
				counts = line.split(",")[2:]
				max_count = max(map(int,counts))
				max_char = filter(lambda x: int(x) == max_count, counts)
				index = counts.index(max_char[0])
				majority_char = dictionary[index]
				conseq += majority_char
				f.write("{},{},{}\n".format(query_pos,region_pos,",".join(counts)))
	
			# Close the new slice matrix, and write the consensus sequence of it
			f.close()
			f_conseq.write(">{}_{}\n{}".format(sample, slice, conseq))
			f_conseq.close()
	
			# If slice contains no data, delete it
			if is_empty:
				os.remove(slice_filename)
				os.remove(conseq_filename)
