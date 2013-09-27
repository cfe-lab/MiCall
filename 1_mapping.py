"""
Automated processing of MiSeq data.  

Do a preliminary map of short reads in paired-end FASTQ files to a
large set of reference sequences.  Split the SAM file output into
multiple SAM files by reference.  Remap all the original FASTQ data
to new reference sequences generated as the consensus of each SAM
file.

Dependencies:
	pileup2conseq_v2.py
	sam2fasta.py
"""

import os,sys
from seqUtils import convert_fasta
from glob import glob
import subprocess
from sam2fasta import *

def remap (f1, f2, samfile, ref, qCutoff=30):
	"""
	Generate a sample-specific consensus sequence from a samtools
	PILEUP file, and then remap everything to this consensus as a
	reference sequence.  Returns paths to the SAM output and
	consensus sequence file.
	"""
	bamfile = samfile.replace('.sam', '.bam')
	confile = bamfile+'.pileup.conseq'
	samfile = samfile.replace('.sam', '.remap.sam')

	os.system('samtools view -bt %s.fasta.fai %s > %s 2>/dev/null' % (ref, samfile, bamfile))
	os.system('samtools sort %s %s.sorted' % (bamfile, bamfile))
	os.system('samtools mpileup -A %s.sorted.bam > %s.pileup 2>/dev/null' % (bamfile, bamfile))
	os.system("python pileup2conseq_v2.py {}.pileup {}".format(bamfile, qCutoff))
	os.system('bowtie2-build -q -f %s %s' % (confile, confile))
	os.system('bowtie2 --quiet -p 1 --local -x %s -1 %s -2 %s -S %s --no-unal --met-file %s --un %s --un-conc %s' % (confile, 
			f1, f2, samfile, samfile.replace('.sam', '.bt2_metrics'),
			samfile.replace('.sam', '.bt2_unpaired_noalign.fastq'), 
			samfile.replace('.sam', '.bt2_paired_noalign.fastq')))
	return samfile, confile

refpath = sys.argv[1]
f = sys.argv[2]
qCutoff = int(sys.argv[3])

# Get names of references
infile = open(refpath+'.fasta', 'rU')
refnames = []
for line in infile:
	if line.startswith('>'):
		refnames.append(line.strip('>\n'))
infile.close()

# Path to R1 FASTQ
root = '/'.join(f.split('/')[:-1])
filename = f.split('/')[-1]
prefix = filename.split('.')[0].split('_')[0]

f1 = f
f2 = f.replace('R1', 'R2')

# Initial mapping to construct reference
samfile = '{}/{}.prelim.sam'.format(root, prefix)
os.system('bowtie2 --quiet -p 1 --local -x %s -1 %s -2 %s -S %s' % (refpath, f1, f2, samfile))

# Prepare reference-specific SAM files
refsams = {}
for i, refname in enumerate(refnames):
	refsams.update({refname: {'handle': open('%s/%s.%s.sam' % (root, prefix, refname), 'w'),'count': 0}})
refsams.update({'*': {'handle': open('%s/%s.unmapped.sam' % (root, prefix), 'w'),'count': 0}})

# Split SAM by mapping
infile = open(samfile, 'rU')
for line in infile:
	if line.startswith('@'):
		# transfer header
		for refname in refnames:
			refsams[refname]['handle'].write(line)
		continue
	
	qname, flag, refname, pos, mapq, cigar, rnext, pnext, tlen, seq, qual = line.strip('\n').split('\t')[:11]
	refsams[refname]['handle'].write(line)
	refsams[refname]['count'] += 1
infile.close()

for refname in refnames:
	refsams[refname]['handle'].close()

# Erase original prelim sam
os.remove(samfile)

# remap reads using sample- and target-specific consensus
for refname in refnames:
	if refsams[refname]['count'] == 0 or refname == 'phiX174' or refname == '*':
		continue
	samfile = refsams[refname]['handle'].name
	samfile, confile = remap(f1, f2, samfile, refpath, qCutoff)

# Just for now, keep all the files around!
sys.exit()

# Clean up intermediate bams and sams
bamfiles = glob(root+'/'+prefix+'.*bam*')
for bamfile in bamfiles:
	if bamfile.endswith('.conseq'): continue
	os.remove(bamfile)

# Clean up intermediate SAM files
for refname in refsams.iterkeys(): os.remove(refsams[refname]['handle'].name)
