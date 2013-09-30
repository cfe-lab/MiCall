"""
Distribute pipeline processes across the cluster.
"""

import os,sys
from glob import glob
from generate_hxb2_poly_files import generate_nuc_counts, generate_amino_counts
from mpi4py import MPI
from sam2fasta import apply_cigar
from seqUtils import ambig_dict, convert_csf, convert_fasta, mixture_dict
from datetime import datetime

my_rank = MPI.COMM_WORLD.Get_rank()
nprocs = MPI.COMM_WORLD.Get_size()

refpath = '/usr/local/share/miseq/refs/cfe'     # Consensus B refs
root = sys.argv[1]                              # Path to folder containing fastqs
mode = sys.argv[2]                              # Nextera or Amplicon - from SampleSheet.csv
qCutoff = int(sys.argv[3])                      # INCORRECT - FIXME

HXB2_mapping_cutoff = 10

amino_alphabet = 'ACDEFGHIKLMNPQRSTVWY*-'
conB_refpath="/usr/local/share/miseq/refs/cfe"
hxb2refpath="/usr/local/share/miseq/refs/"      # hxb2_pol_ref.fasta

# Set debug to 1 and modify debug parameters to test development on this pipeline
debug = 1
debug_root = "/data/miseq/0_testing"
debug_seq = "35759A"

def timestamp(msg): print '[{}] (rank={}/{}) {}'.format(datetime.now().strftime('%Y-%m-%d %H:%M:%S'),my_rank,nprocs,msg)

# Map and remap each fastq to consensus B
files = glob(root + '/*R1*.fastq')
for i in range(len(files)):
	if i % nprocs != my_rank: continue
	filename = files[i].split('/')[-1]
	timestamp('1_mapping {}'.format(filename))
	os.system("python 1_mapping.py {} {} {}".format(refpath, files[i], qCutoff))
MPI.COMM_WORLD.Barrier()
if my_rank == 0: timestamp('Barrier #1 passed\n')

sys.exit()

# For each remapped SAM, generate CSFs (done)
files = glob(root + '/*.remap.sam')
for i in range(len(files)):
	if i % nprocs != my_rank: continue
	filename = files[i].split('/')[-1]

	# For each q-cutoff, generate fasta for amplicon, OR a csf file for Nextera
	for qcut in [0, 10, 15, 20]:
		timestamp('2_sam2fasta {} {} {}'.format(filename, qcut, mode))
		os.system('python 2_sam2fasta_with_base_censoring.py {} {} {}'.format(files[i], qcut, mode))
MPI.COMM_WORLD.Barrier()
if my_rank == 0: timestamp('Barrier #2 passed\n')

# Fpr amplicon sequencing runs, compute g2p scores for env-mapped FASTAs
if mode == 'Amplicon':
	files = glob(root + '/*env*.fasta')
	for i, file in enumerate(files):
		if i % nprocs != my_rank: continue
		filename = files[i].split('/')[-1]
		timestamp('3_g2p_scoring {}'.format(filename))
		os.system('python 3_g2p_scoring.py {}'.format(file))
MPI.COMM_WORLD.Barrier()
if my_rank == 0: timestamp('Barrier #3 passed\n')



# Take remapped sam, generate pileup, call pileup2conseq to generate remap conseq, convert to an indexed .bt2
def consensus_from_remapped_sam(root,ref,samfile,qCutoff=30):
	bamfile = samfile.replace('.sam', '.bam')
	confile = bamfile+'.pileup.conseq'
	os.system('samtools view -bt {}.fasta.fai {} > {} 2>/dev/null'.format(ref, samfile, bamfile))
	os.system('samtools sort {} {}.sorted'.format(bamfile, bamfile))
	os.system('samtools mpileup -A {}.sorted.bam > {}.pileup 2>/dev/null'.format(bamfile, bamfile))
	os.system('python pileup2conseq_v2.py {}.pileup {}'.format(bamfile, qCutoff))	# From pileup, generate poly + conseq
	os.system('bowtie2-build -q -f {} {}'.format(confile, confile))
	print "\n\nMaking remap conseq (indexed .bt2) for {}".format(samfile)

# Create remap conseq from remap sam
remapped_sams = glob(root + '/*.remap.sam')
if debug == 1: remapped_sams = glob('{}/{}-PR-RT.HIV1B-pol.remap.sam'.format(debug_root,debug_seq))

for i in range(len(remapped_sams)):
	if i % nprocs != my_rank: continue
	remapped_sam = remapped_sams[i]
	remap_conseq = remapped_sam.replace('.sam', '.bam.pileup.conseq')
	consensus_from_remapped_sam(root,remap_conseq,remapped_sam,20)
MPI.COMM_WORLD.Barrier()
if my_rank == 0: timestamp('Barrier #4 passed\n')



# From remap csf, create remap fastas to be aligned against remap conseq
csf_files = glob(root + '/*.csf')
if debug == 1: csf_files = glob('{}/{}-PR-RT.HIV1B-pol.20.csf'.format(debug_root,debug_seq))

for i in range(len(csf_files)):
	if i % nprocs != my_rank: continue
	f = csf_files[i]
	fasta_filename = f.replace('.csf', '.csf.fasta')
	prefix, gene = fasta_filename.split('.')[:2]
	conseq_filename = "{}.{}.remap.bam.pileup.conseq".format(prefix,gene)
	sam_filename = f.replace('.csf', '.final.sam')

	# Write temp remap fasta (Remove gaps so remap fasta can be aligned with remap conseq)
	infile = open(f, 'rU')
	fasta, left_gaps, right_gaps = convert_csf(infile)	# Returns left-padded sequences
	infile.close()
	print "Making (temp) remap fasta from {}".format(f)
	outfile = open(fasta_filename, 'w')
	for j, (h, s) in enumerate(fasta): outfile.write(">{}\n{}\n".format(h, s.strip("-")))
	outfile.close()

	# Generate final sam by aligning fastas against remap conseq, then delete the temp fasta
	# -f 	Inputs are fasta	-p 1	Number of alignment threads
	# -x	bt2 index		-U	Fasta with unpaired reads
	# -S	samfile output		--un	Where stuff didn't map		--met-file	Metrics output
	command = 'bowtie2 -f --quiet -p 1 -x {} -U {} -S {} --no-unal --met-file {} --un {}'.format(
			conseq_filename,fasta_filename, sam_filename,
			f.replace('.csf', '.csf.bt2_metrics'),
			f.replace('.csf', '.csf.bt2_unpaired_noalign.fastq'))
	print "\nMaking final sam: {}".format(command)
	os.system(command)
	os.remove(fasta_filename)
MPI.COMM_WORLD.Barrier()
if my_rank == 0: timestamp('Barrier #5 passed\n')



# From final sam, extract final CSFs
final_sams = glob(root + '/*.final.sam')
if debug == 1: final_sams = glob('{}/{}-PR-RT.HIV1B-pol.20.final.sam'.format(debug_root,debug_seq))

for i in range(len(final_sams)):
	if i % nprocs != my_rank: continue
	filename  = final_sams[i].split('/')[-1]
	sample, region, qCut = filename.split('.')[:3]
	command = 'python 2_sam2fasta_with_base_censoring.py {} {} {}'.format(final_sams[i], qCut, "Nextera")
	print "Making final csf: {}".format(command)
	os.system(command)
MPI.COMM_WORLD.Barrier()
if my_rank == 0: timestamp('Barrier #6 passed\n')

# HXB2 align final csfs
final_csfs = glob(root + '/*.final.sam.20.csf')
if debug == 1: final_sams = glob('{}/{}-PR-RT.HIV1B-pol.20.final.sam.20.csf'.format(debug_root,debug_seq))

for i in range(len(final_csfs)):
	if i % nprocs != my_rank: continue
	f = final_csfs[i]
	sample, region, qCut = f.split('.')[:3]
	fasta_filename = f.replace('.csf', '.fasta')
	sam_filename = f.replace('.csf', '.HXB2.sam')
	hxb2_reference = "{}{}_hxb2.fasta".format(hxb2refpath,region)   # /usr/local/share/miseq/refs/HIV1B-pol_hxb2.fasta

	# From final csf, generate temp final fasta with gaps removed (IE, suitable for HXB2 alignment)
	print "Making (temp) final fasta from {}".format(f)
	infile = open(f, 'rU')
	fasta, left_gaps, right_gaps = convert_csf(infile)
	infile.close()
	outfile = open(fasta_filename, 'w')
	for j, (h, s) in enumerate(fasta): outfile.write(">{}\n{}\n".format(h, s.strip("-")))
	outfile.close()

	# Generate HXB2 sam by aligning final fasta against HXB2, delete the temp final fasta
	command = 'bowtie2 -f --quiet -p 1 -x {} -U {} -S {} --no-unal --met-file {} --un {}'.format(
		hxb2_reference,fasta_filename, sam_filename,
		f.replace('.csf', '.csf.bt2_metrics'),
		f.replace('.csf', '.csf.bt2_unpaired_noalign.fastq'))
        print "\nMaking HXB2 sam: {}".format(command)
	os.system(command)
	os.remove(fasta_filename)
MPI.COMM_WORLD.Barrier()
if my_rank == 0: timestamp('Barrier #7 passed\n')


# Extract HXB2 aligned sequences from the HXB2.sam and generate amino + nuc poly files
HXB2_sams = glob(root + '/*.HXB2.sam')
if debug == 1: HXB2_sams = glob('{}/{}-PR-RT.HIV1B-pol.20.final.sam.20.HXB2.sam'.format(debug_root,debug_seq))
for i in range(len(HXB2_sams)):
	if i % nprocs != my_rank: continue
	generate_nuc_counts(HXB2_sams[i],HXB2_mapping_cutoff)
	generate_amino_counts(HXB2_sams[i],HXB2_mapping_cutoff)

MPI.COMM_WORLD.Barrier()
if my_rank == 0: timestamp('Barrier #8 passed\n')


# Step 4B: Generate an amino summary file
if my_rank == 0:
	amino_polys = glob(root + '*.HXB2.amino_poly')
	amino_summary = root + "/HXB2.amino_poly.summary"
	print "Dumping amino matrix data into {}".format(amino_summary)
	os.system("echo \"Sample,AA.pos,A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,*,-\" > {}".format(amino_summary))
	for i in range(len(amino_polys)): os.system("tail -n +2 {} >> {}".format(amino_polys[i], amino_summary))
MPI.COMM_WORLD.Barrier()
if my_rank == 0: timestamp('Barrier #9 passed\n')

MPI.Finalize()
