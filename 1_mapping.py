"""
1) Map fastq reads to consensus B reference sequences
	--> Results stored in prelim.sam
2) Split the prelim.sam into region-specific SAMs
3) Make a sample/region specific consensus
4) Remap fastq reads against sample/region specific references
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
	remapped_sam = samfile.replace('.sam', '.remap.sam')

	os.system('samtools view -bt %s.fasta.fai %s > %s 2>/dev/null' % (ref, samfile, bamfile))
	os.system('samtools sort %s %s.sorted' % (bamfile, bamfile))
	os.system('samtools mpileup -A %s.sorted.bam > %s.pileup 2>/dev/null' % (bamfile, bamfile))
	os.system("python pileup2conseq_v2.py {}.pileup {}".format(bamfile, qCutoff))
	os.system('bowtie2-build -q -f %s %s' % (confile, confile))
	os.system('bowtie2 --quiet -p 1 --local -x %s -1 %s -2 %s -S %s --no-unal --met-file %s --un %s --un-conc %s' % (confile, 
			f1, f2, remapped_sam, remapped_sam.replace('.sam', '.bt2_metrics'),
			remapped_sam.replace('.sam', '.bt2_unpaired_noalign.fastq'), 
			remapped_sam.replace('.sam', '.bt2_paired_noalign.fastq')))
	return samfile, confile

refpath = sys.argv[1]
f = sys.argv[2]
qCutoff = 20

# Get region codes from the con B reference fasta headers
infile = open(refpath+'.fasta', 'rU')
refnames = []
for line in infile:
	if line.startswith('>'): refnames.append(line.strip('>\n'))
infile.close()

# Path to R1 FASTQ
root = '/'.join(f.split('/')[:-1])
filename = f.split('/')[-1]
prefix = filename.split('.')[0].split('_')[0]
f1 = f
f2 = f.replace('R1', 'R2')

# Initial consensus B mapping
samfile = '{}/{}.prelim.sam'.format(root, prefix)
command = 'bowtie2 --quiet -p 1 --local -x %s -1 %s -2 %s -S %s' % (refpath, f1, f2, samfile)
os.system(command)

# Create region-specific SAMs
refsams = {}
for i, refname in enumerate(refnames):
	refsams.update({refname: {'handle': open('%s/%s.%s.sam' % (root, prefix, refname), 'w'),'count': 0}})
refsams.update({'*': {'handle': open('%s/%s.unmapped.sam' % (root, prefix), 'w'),'count': 0}})

# Then, subdivide the prelim SAM into these region-specific SAMs
infile = open(samfile, 'rU')
for line in infile:
	# Copy the SAM header into each file
	if line.startswith('@'):
		for refname in refnames: refsams[refname]['handle'].write(line)
		continue
	qname, flag, refname, pos, mapq, cigar, rnext, pnext, tlen, seq, qual = line.strip('\n').split('\t')[:11]
	refsams[refname]['handle'].write(line)
	refsams[refname]['count'] += 1
infile.close()

for refname in refnames:
	refsams[refname]['handle'].close()

# FIXME: Q CUTOFFS NOT WORKING RIGHT NOW
# Remap fastqs using sample/region specific consensus reference
for refname in refnames:
	if refsams[refname]['count'] == 0 or refname == 'phiX174' or refname == '*':
		continue
	samfile = refsams[refname]['handle'].name
	samfile, confile = remap(f1, f2, samfile, refpath, 20)

# Clean up intermediate bams and sams
#bamfiles = glob(root+'/'+prefix+'.*bam*')
#for bamfile in bamfiles:
#	if bamfile.endswith('.conseq'): continue
#	os.remove(bamfile)

# Clean up intermediate SAM files
#for refname in refsams.iterkeys(): os.remove(refsams[refname]['handle'].name)
