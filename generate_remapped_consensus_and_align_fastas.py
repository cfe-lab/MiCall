"""
Generate final sams, poly, and consensus sequence files.

1A) From remap sam, generate remap consensus + bt files
1B) From remap sam, generate persistent CSFs + temporary fastas
1C) Generate final sam by aligning temp fastas against remap conseq
2) Generate final poly file from final sam
3) Generate final conseq from the poly file

INPUT: Folder containing remapped sams
OUTPUT: Final sam/poly/conseq
"""
import os
import sys
from glob import glob
from seqUtils import convert_csf, mixture_dict, ambig_dict
from mpi4py import MPI
my_rank = MPI.COMM_WORLD.Get_rank()
nprocs = MPI.COMM_WORLD.Get_size()

# Takes sam as input, generates pileup, calls pileup2conseq to generate consensus
# Then generates an indexed consensus using bowtie2-build
def consensus_from_remapped_sam(root,ref,samfile,qCutoff=30):
	bamfile = samfile.replace('.sam', '.bam')
	confile = bamfile+'.pileup.conseq'
	os.system('samtools view -bt {}.fasta.fai {} > {} 2>/dev/null'.format(ref, samfile, bamfile))
	os.system('samtools sort {} {}.sorted'.format(bamfile, bamfile))
	os.system('samtools mpileup -A {}.sorted.bam > {}.pileup 2>/dev/null'.format(bamfile, bamfile))
	os.system('python pileup2conseq_v2.py {}.pileup {}'.format(bamfile, qCutoff))	# From pileup, generate poly + conseq
	os.system('bowtie2-build -q -f {} {}'.format(confile, confile))
	print "Generated conseq for {}".format(samfile)

debug = 0
#debug_seq = "35759A"

refpath="/usr/local/share/miseq/refs/cfe"
root = sys.argv[1]

# Step 1A: From remapped sams, generate remap conseqs
remapped_sams = glob(root + '/*.remap.sam')
if debug == 1:
	remapped_sams = glob('/data/miseq/0_testing/{}-PR-RT.HIV1B-pol.remap.sam'.format(debug_seq))

for i in range(len(remapped_sams)):
	if i % nprocs != my_rank:
		continue
	remapped_sam = remapped_sams[i]

	# From remapped sams, generate remap conseq
	remap_conseq = remapped_sam.replace('.sam', '.bam.pileup.conseq')
	consensus_from_remapped_sam(root,remap_conseq,remapped_sam,20)
MPI.COMM_WORLD.Barrier()

# Step 1B: From csfs, create fastas to be aligned against remap conseq
csf_files = glob(root + '/*.csf')
if debug == 1:
	csf_files = glob('/data/miseq/0_testing/{}-PR-RT.HIV1B-pol.20.csf'.format(debug_seq))

for i in range(len(csf_files)):
	if i % nprocs != my_rank:
		continue
	f = csf_files[i]

	# convert_csf returns left-padded sequences
	fasta_filename = f.replace('.csf', '.csf.fasta')
	infile = open(f, 'rU')
	fasta, left_gap_position, right_gap_position = convert_csf(infile)
	infile.close()

	# Remove gaps so fasta is ready for alignment with remap conseq
	outfile = open(fasta_filename, 'w')
	for j, (h, s) in enumerate(fasta):
		outfile.write(">{}\n{}\n".format(h, s.strip("-")))
	outfile.close()

	prefix, gene = fasta_filename.split('.')[:2]
	conseq_filename = "{}.{}.remap.bam.pileup.conseq".format(prefix,gene)
	sam_filename = f.replace('.csf', '.csf.final.sam')

	# Step 1C: Generate final sam by aligning fastas against remap conseq

	# BOWTIE2 PARAMETERS
	# -f 	Inputs are fasta	-p 1	Number of alignment threads
	# -x	bt2 index		-U	Fasta with unpaired reads
	# -S	samfile output		--un	Where stuff didn't map		--met-file	Metrics output
	command = 'bowtie2 -f --quiet -p 1 -x {} -U {} -S {} --no-unal --met-file {} --un {}'.format(
			conseq_filename,fasta_filename, sam_filename,
			f.replace('.csf', '.csf.bt2_metrics'),
			f.replace('.csf', '.csf.bt2_unpaired_noalign.fastq'))
	print "Generating final sam: {}".format(command)
	os.system(command)

MPI.COMM_WORLD.Barrier()

def poly_from_final_sam(root,ref,final_sam_file,qCutoff=30):
	bamfile = final_sam_file.replace('.sam', '.bam')
	os.system('samtools view -bt {}.fasta.fai {} > {} 2>/dev/null'.format(ref, final_sam_file, bamfile))
	os.system('samtools sort {} {}.sorted'.format(bamfile, bamfile))
	os.system('samtools mpileup -A {}.sorted.bam > {}.pileup 2>/dev/null'.format(bamfile, bamfile))
	os.system('python pileup2conseq_v2.py {}.pileup {}'.format(bamfile, qCutoff))   # Generates poly

# Step 2: For each final sam, generate poly needed for making final conseq
final_sams = glob(root + '/*.final.sam')
if debug == 1:
	final_sams = glob('/data/miseq/0_testing/{}-PR-RT.HIV1B-pol.20.csf.final.sam'.format(debug_seq))

for i in range(len(final_sams)):
	if i % nprocs != my_rank:
		continue
	final_sam = final_sams[i]
	poly_from_final_sam(root,refpath,final_sam,20)

# Step 3: For each poly, generate final conseq
poly_files = glob(root + '/*.remap.*.poly')

if debug == 1:
	poly_files = glob('/data/miseq/0_testing/{}-PR-RT.HIV1B-pol.20.csf.final.bam.pileup.poly'.format(debug_seq))

for j in range(len(poly_files)):
	if j % nprocs != my_rank:
		continue
	poly_file = poly_files[j]
	infile = open(poly_file, 'rU')

	seq = ""
	for i, line in enumerate(infile):

		if i == 0:	# Ignore poly header
			continue

		line = line.strip("\n")
		freqs = {}
		coord, freqs['A'], freqs['C'], freqs['G'], freqs['T'], freqs['N'] = [int(n) for n in line.split(',')]
		base = max(freqs, key=lambda n: freqs[n])
		max_count = freqs[base]

		# possib stores all bases having the highest count
		possib = filter(lambda n: freqs[n] == max_count, freqs)

		# Convert base to correct mixture character using ambig_dict
		if len(possib) == 1:
			base = possib[0]
		else:
			base = ambig_dict["".join(sorted(possib))]
		seq += base

	infile.close()

	final_conseq = poly_files[j] + ".final_conseq"
	f = open(final_conseq, 'w')
	prefix, gene, q = final_conseq.split('.')[:3]
	fasta_header = "{}.{}.{}".format(prefix,gene,q)
	f.write(">{}\n{}".format(fasta_header,seq))
	f.close()
	print "Writing final conseq: {}".format(final_conseq)

MPI.COMM_WORLD.Barrier()

# Aggregate conseq files into _final_conseqs.fasta
if my_rank == 0:
	print "Aggregating *.final_conseq into _final_conseqs.fasta"
	f = open(root + "/_final_conseqs.fasta", 'w')
	final_conseqs = glob(root + '/*.final_conseq')
	for j in range(len(final_conseqs)):
		conseq_contents = open(final_conseqs[j]).read()
		f.write(conseq_contents)
	f.close()
