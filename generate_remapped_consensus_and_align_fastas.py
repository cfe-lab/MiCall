"""
Generates final sam files.
1) From remap sam, generate remap consensus + bt files
2) From remap sam, generate persistent CSFs + temporary fastas
3) Generate final sam by aligning temporary fastas against remap consensus
4) Generate final poly file from final sam
5) Generate final consensus from the poly file

INPUT: Folder of interest containing remapped sam files
OUTPUT: Final sam files + final poly + final consensus
"""

import os
import sys
from glob import glob
from seqUtils import convert_csf

debug = 1
debug_seq = "35854A"

from mpi4py import MPI
my_rank = MPI.COMM_WORLD.Get_rank()
nprocs = MPI.COMM_WORLD.Get_size()

def consensus_from_remapped_sam(root,ref,samfile,qCutoff=30):
	"""Generate remap conseq + remap bt"""
	bamfile = samfile.replace('.sam', '.bam')

	command = 'samtools view -bt {}.fasta.fai {} > {} 2>/dev/null'.format(ref, samfile, bamfile)
	os.system(command)

	command = 'samtools sort {} {}.sorted'.format(bamfile, bamfile)
	os.system(command)

	command = 'samtools mpileup -A {}.sorted.bam > {}.pileup 2>/dev/null'.format(bamfile, bamfile)
	os.system(command)

	# Also generates a poly file corresponding to this remapped sam
	command = 'python pileup2conseq_v2.py {}.pileup {}'.format(bamfile, qCutoff)
	os.system(command)

	confile = bamfile+'.pileup.conseq'
	command = 'bowtie2-build -q -f {} {}'.format(confile, confile)
	os.system(command)

	print "Generated conseq for {}".format(samfile)


def poly_from_final_sam(root,ref,final_sam_file,qCutoff=30):
	bamfile = final_sam_file.replace('.sam', '.bam')

	command = 'samtools view -bt {}.fasta.fai {} > {} 2>/dev/null'.format(ref, final_sam_file, bamfile)
	os.system(command)

	command = 'samtools sort {} {}.sorted'.format(bamfile, bamfile)
	os.system(command)

	command = 'samtools mpileup -A {}.sorted.bam > {}.pileup 2>/dev/null'.format(bamfile, bamfile)
	os.system(command)

	# Generate a poly file from this final sam
	command = 'python pileup2conseq_v2.py {}.pileup {}'.format(bamfile, qCutoff)
	os.system(command)

	print "Generating poly from {}.pileup".format(bamfile)


refpath="/usr/local/share/miseq/refs/cfe"
root = sys.argv[1]

# For each remapped sam, generate an indexed consensus (bt file)
samfiles = glob(root + '/*.remap.sam')
if debug == 1:
	samfiles = glob('/data/miseq/130820_M01841_0018_000000000-A4V8D/{}-PR-RT.HIV1B-pol.remap.sam'.format(debug_seq))

for i in range(len(samfiles)):
	if i % nprocs != my_rank:
		continue
	samfile = samfiles[i]

	# FIXME: Currently refpath is wrong - it is using con B, but should be using sample specific .pileup.conseq from 1_mapping
	consensus_from_remapped_sam(root,refpath,samfile,20)
MPI.COMM_WORLD.Barrier()

# Convert each csf to fasta for alignment against remap consensus
fasta = ""
left_gap_position = {}
right_gap_position = {}

csf_files = glob(root + '/*.csf')

if debug == 1:
	csf_files = glob('/data/miseq/130820_M01841_0018_000000000-A4V8D/{}-PR-RT.HIV1B-pol.20.csf'.format(debug_seq))

# For each csf, generate fasta with left-padded sequences
for i in range(len(csf_files)):
	if i % nprocs != my_rank:
		continue
	f = csf_files[i]

	# convert_csf returns left-padded sequences
	fasta_filename = f.replace('.csf', '.csf.fasta')
	infile = open(f, 'rU')
	fasta, left_gap_position, right_gap_position = convert_csf(infile)
	infile.close()

	# Write this fasta to disk
	print "Writing fasta: {}".format(fasta_filename)
	outfile = open(fasta_filename, 'w')

	# Remove gaps so fasta is ready for alignment with remap conseq
	for j, (h, s) in enumerate(fasta):
		s = s.strip("-")
		outfile.write(">{}\n{}\n".format(h, s))
	outfile.close()

	# Determine filename of corresponding remap conseq + final.sam file
	prefix, gene = fasta_filename.split('.')[:2]
	conseq_filename = "{}.{}.remap.bam.pileup.conseq".format(prefix,gene)
	sam_filename = f.replace('.csf', '.csf.final.sam')

	# Generate final.sam: align fastas against remap conseq with bowtie2
	# -f 	Inputs are fasta	-p 1	Number of alignment threads
	# -x	bt2 index		-U	Fasta with unpaired reads
	# -S	samfile output		--un	Where stuff didn't map		--met-file	Metrics output
	command = 'bowtie2 -f --quiet -p 1 -x {} -U {} -S {} --no-unal --met-file {} --un {}'.format(
			conseq_filename,fasta_filename, sam_filename,
			f.replace('.csf', '.csf.bt2_metrics'),
			f.replace('.csf', '.csf.bt2_unpaired_noalign.fastq'))
	print "Aligning: {}".format(command)
	os.system(command)

	# Delete fasta intermediary
	#os.remove(fasta_filename)

MPI.COMM_WORLD.Barrier()

# For each prefix.regioncode.qScore.csf.final.sam, generate the poly file
final_sams = glob(root + '/*.final.sam')
if debug == 1:
	final_sams = glob(root + '/{}-PR-RT.HIV1B-pol.20.csf.final.sam'.format(debug_seq))

for i in range(len(final_sams)):
	if i % nprocs != my_rank:
		continue
	samfile = final_sams[i]

	poly_from_final_sam(root,refpath,samfile,20)

# Generate conseq files
poly_files = glob(root + '/*.remap.*.poly')

if debug == 1:
	poly_files = glob(root + '/{}-PR-RT.HIV1B-pol.20.csf.final.bam.pileup.poly'.format(debug_seq))

for j in range(len(poly_files)):
	if j % nprocs != my_rank:
		continue
	f = poly_files[j]

	# Open the poly file, generate the conseq
	infile = open(f, 'rU')

	consensus_sequence = ""

	for i, line in enumerate(infile):

		# Ignore the header
		if i == 0:
			continue

		line = line.strip("\n")
		freqs = {}
		coord, freqs['A'], freqs['C'], freqs['G'], freqs['T'], freqs['N'] = [int(n) for n in line.split(',')]
		base = max(freqs, key=lambda n: freqs[n])
		max_count = freqs[base]
		# Return all bases (elements of freqs) such that the freq(b) = max_count
		possib = filter(lambda n: freqs[n] == max_count, freqs)
		consensus_sequence += base

	print "{}".format(consensus_sequence)
	infile.close()

MPI.COMM_WORLD.Barrier()
