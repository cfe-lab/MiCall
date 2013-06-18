"""
Automated processing of MiSeq data.  

Do a preliminary map of short reads in paired-end FASTQ files to a
large set of reference sequences.  Split SAM file output into
multiple SAM files by reference.  Remap all original FASTQ data
to new reference sequences generated as the consensus of each SAM
file.

Dependencies:
	pipeline2conseq_v2.py
	sam2fasta.py
"""

import sys
import os
from seqUtils import convert_fasta
from glob import glob
import subprocess
from sam2fasta import *

if len(sys.argv) != 3:
	print 'Usage: python pipeline3.py /path/to/*R1*.fastq 20'
	sys.exit()

root = sys.argv[1]
qCutoffForBaseCensor = sys.argv[2]

# Reference sequences for bowtie2-build
ref = '/Users/emartin/Desktop/Art_MiSeq_Processing_Scripts/cfe_reference_sequences/cfe'

# Get names of references
infile = open(ref+'.fasta', 'rU')
refnames = []
for line in infile:
	if line.startswith('>'):
		refnames.append(line.strip('>\n'))

infile.close()

def remap (f1, f2, samfile, ref):
	"""
	Generate a sample-specific consensus sequence from a samtools
	PILEUP file, and then remap everything to this consensus as a
	reference sequence.  Returns paths to the SAM output and
	consensus sequence file.
	"""
	bamfile = samfile.replace('.sam', '.bam')
	os.system('samtools view -bt %s.fasta.fai %s > %s 2>/dev/null' % (ref, samfile, bamfile))
	os.system('samtools sort %s %s.sorted' % (bamfile, bamfile))
	os.system('samtools mpileup -A %s.sorted.bam > %s.pileup 2>/dev/null' % (bamfile, bamfile))
	
	# generate consensus sequence from pileup
	os.system('python pileup2conseq_v2.py %s.pileup' % bamfile)
	confile = bamfile+'.pileup.conseq'
	os.system('bowtie2-build -q -f %s %s' % (confile, confile))
	
	# mapping against new reference
	samfile = samfile.replace('.sam', '.remap.sam')
	
	os.system('bowtie2 --quiet -p 6 --local -x %s -1 %s -2 %s -S %s\
				--no-unal --met-file %s --un %s --un-conc %s' % (confile, 
			f1, f2, samfile, samfile.replace('.sam', '.bt2_metrics'),
			samfile.replace('.sam', '.bt2_unpaired_noalign.fastq'), 
			samfile.replace('.sam', '.bt2_paired_noalign.fastq')))
	
	return samfile, confile




files = glob(root+'/*R1*.fastq')

for f in files:
	prefix = f.split('/')[-1].split('.')[0].split('_')[0]
	f1 = f
	f2 = f.replace('R1', 'R2')
	
	# initial mapping to construct reference
	print prefix, 'preliminary mapping'
	samfile = '%s/%s.prelim.sam' % (root, prefix)
	os.system('bowtie2 --quiet -p 6 --local -x %s -1 %s -2 %s -S %s' % (ref, 
			f1, 
			f2, 
			samfile))
	
	# prepare reference-specific SAM files
	refsams = {}
	for i, refname in enumerate(refnames):
		refsams.update({refname: {'handle': open('%s/%s.%s.sam' % (root, prefix, refname), 'w'),
									'count': 0}})
	
	refsams.update({'*': {'handle': open('%s/%s.unmapped.sam' % (root, prefix), 'w'),
							'count': 0}})
	
	# split SAM by mapping
	print 'splitting initial map by reference'
	infile = open(samfile, 'rU')
	
	for line in infile:
		if line.startswith('@'):
			# transfer header
			for refname in refnames:
				refsams[refname]['handle'].write(line)
			continue
		
		qname, flag, refname, pos, mapq, cigar, rnext, pnext, tlen, seq, qual = line.strip('\n').split('\t')[:11]
		if refname == '*':
			continue
		refsams[refname]['handle'].write(line)
		refsams[refname]['count'] += 1
	
	infile.close()
	for refname in refnames:
		refsams[refname]['handle'].close()
	
	
	# remap reads using sample- and target-specific consensus
	for refname in refnames:
		if refsams[refname]['count'] == 0 or refname == 'phiX174' or refname == '*':
			continue
		
		samfile = refsams[refname]['handle'].name
		
		if refsams[refname]['count'] > 10:
			print prefix, 'remapping on', refname
			samfile, confile = remap(f1, f2, samfile, ref)
		
		print prefix, 'generating FASTA for', refname
		infile = open(samfile, 'rU')
		fasta = sam2fasta(infile, cutoff=qCutoffForBaseCensor)
		infile.close()
		
		if fasta:
			outfile = open(samfile + '.fasta', 'w')
			outfile.write(fasta)
			outfile.close()
	

