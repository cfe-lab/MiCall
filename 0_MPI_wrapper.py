"""
Distribute pipeline processes across the cluster.
"""

import os,sys
from glob import glob
from generate_hxb2_poly_files import generate_amino_counts, generate_nuc_counts
from mpi4py import MPI
from sam2fasta import apply_cigar
from seqUtils import ambig_dict, convert_csf, convert_fasta, mixture_dict, timestamp

my_rank = MPI.COMM_WORLD.Get_rank()
nprocs = MPI.COMM_WORLD.Get_size()

## settings
refpath = '/usr/local/share/miseq/refs/cfe'	# Consensus B refs
hxb2refpath="/usr/local/share/miseq/refs/" 	# location of HXB2 reference files
HXB2_mapping_cutoff = 10			# Read mapping quality cutoff
amino_conseq_cutoff = 10			# Minimum occurences of char before it counts for conseq generation

## arguments
root = sys.argv[1]				# Path to folder containing fastqs
mode = sys.argv[2]				# Nextera or Amplicon - from SampleSheet.csv
qCutoff = int(sys.argv[3])			# INCORRECT - FIXME




def consensus_from_remapped_sam(root, ref, samfile, qCutoff=30):
	"""
	Take remapped sam, generate pileup, call pileup2conseq to generate 
	remap conseq, convert to an indexed .bt2 so that it can be used as 
	a reference for bowtie2.
	"""
	bamfile = samfile.replace('.sam', '.bam')
	confile = bamfile+'.pileup.conseq'
	
	os.system('samtools view -bt {}.fasta.fai {} > {} 2>/dev/null'.format(ref, samfile, bamfile))
	os.system('samtools sort {} {}.sorted'.format(bamfile, bamfile))
	os.system('samtools mpileup -A {}.sorted.bam > {}.pileup 2>/dev/null'.format(bamfile, bamfile))
	
	# From pileup, generate poly + conseq
	os.system('python pileup2conseq_v2.py {}.pileup {}'.format(bamfile, qCutoff))
	os.system('bowtie2-build -q -f {} {}'.format(confile, confile))


#######################################

# Map and remap each fastq to consensus B
files = glob(root + '/*R1*.fastq')
for i, file in enumerate(files):
	if i % nprocs != my_rank:
		continue
	filename = file.split('/')[-1]
	timestamp('1_mapping {} {} {}'.format(refpath, file, qCutoff), my_rank, nprocs)
	os.system("python 1_mapping.py {} {} {}".format(refpath, file, qCutoff))

MPI.COMM_WORLD.Barrier()
if my_rank == 0: timestamp('Barrier #1 (Prelim + remapping)\n', my_rank, nprocs)



# For each remapped SAM, generate CSFs (done)
files = glob(root + '/*.remap.sam')
for i in range(len(files)):
	if i % nprocs != my_rank:
		continue
	filename = files[i].split('/')[-1]
	for qcut in [0, 10, 15, 20]:
		timestamp('2_sam2fasta {} {} {}'.format(filename, qcut, mode), my_rank, nprocs)
		os.system('python 2_sam2fasta_with_base_censoring.py {} {} {}'.format(files[i], qcut, mode))

MPI.COMM_WORLD.Barrier()
if my_rank == 0: timestamp('Barrier #2 (Remap CSF from remap SAM)\n', my_rank, nprocs)



# For amplicon sequencing runs, compute g2p scores for env-mapped FASTAs
if mode == 'Amplicon':
	files = glob(root + '/*env*.fasta')
	for i, file in enumerate(files):
		if i % nprocs != my_rank:
			continue
		filename = files[i].split('/')[-1]
		timestamp('3_g2p_scoring {}'.format(filename), my_rank, nprocs)
		os.system('python 3_g2p_scoring.py {}'.format(file))

MPI.COMM_WORLD.Barrier()

if mode == 'Amplicon' and my_rank == 0: timestamp('Barrier #3 (G2P - amplicon only)\n', my_rank, nprocs)
elif mode != 'Amplicon' and my_rank == 0: timestamp('Barrier #3 (Not an amplicon run - skipped G2P)\n', my_rank, nprocs)



# Create remap conseq from remap sam
remapped_sams = glob(root + '/*.remap.sam')

for i, remapped_sam in enumerate(remapped_sams):
	if i % nprocs != my_rank:
		continue
		
	remap_conseq = remapped_sam.replace('.sam', '.bam.pileup.conseq')
	timestamp('consensus_from_remapped_sam({},{},{},{})'.format(root, remap_conseq, 
		remapped_sam, qCutoff), my_rank, nprocs)
	consensus_from_remapped_sam(root, remap_conseq, remapped_sam, qCutoff)

MPI.COMM_WORLD.Barrier()
if my_rank == 0: timestamp('Barrier #4 (Remap conseq)\n', my_rank, nprocs)



# From remap csf, create remap fastas to be aligned against remap conseq
if mode == 'Nextera':
	csf_files = glob(root + '/*.csf')

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
		outfile = open(fasta_filename, 'w')
		for j, (h, s) in enumerate(fasta): outfile.write(">{}\n{}\n".format(h, s.strip("-")))
		outfile.close()

		# Generate final sam by aligning fastas against remap conseq, then delete the temp fasta
		# -f 	Inputs are fasta	-p 1	Number of alignment threads
		# -x	bt2 index		-U	Fasta with unpaired reads
		# -S	samfile output		--un	Where stuff didn't map		--met-file	Metrics output
		command = 'bowtie2 -f --local --quiet -p 1 -x {} -U {} -S {} --no-unal --met-file {} --un {}'.format(
				conseq_filename,fasta_filename, sam_filename,
				f.replace('.csf', '.csf.bt2_metrics'),
				f.replace('.csf', '.csf.bt2_unpaired_noalign.fastq'))
		
		timestamp("Making final sam: {}".format(command), my_rank, nprocs)
		os.system(command)
		os.remove(fasta_filename)
	
	MPI.COMM_WORLD.Barrier()
	if my_rank == 0: timestamp('Barrier #5 (Final SAM)\n', my_rank, nprocs)



# From final sam, extract final CSFs
final_sams = glob(root + '/*.final.sam')

for i, final_sam in enumerate(final_sams):
	if i % nprocs != my_rank:
		continue
	filename = final_sam.split('/')[-1]
	sample, region, qCut = filename.split('.')[:3]
	
	command = 'python 2_sam2fasta_with_base_censoring.py {} {} {}'.format(final_sams[i], qCut, "Nextera")
	timestamp("Making final csf: {}".format(command), my_rank, nprocs)
	os.system(command)

MPI.COMM_WORLD.Barrier()
if my_rank == 0: timestamp('Barrier #6 (CSF from final SAM)\n', my_rank, nprocs)



# HXB2 align final csfs
final_csfs = glob(root + '/*.final.sam.%d.csf' % qCutoff)

for i in range(len(final_csfs)):
	if i % nprocs != my_rank:
		continue
	
	f = final_csfs[i]
	sample, region, qCut = f.split('.')[:3]
	
	fasta_filename = f.replace('.csf', '.fasta')
	sam_filename = f.replace('.csf', '.HXB2.sam')
	hxb2_reference = "{}{}_hxb2.fasta".format(hxb2refpath,region)   # /usr/local/share/miseq/refs/HIV1B-pol_hxb2.fasta

	# From final csf, generate temp final fasta with gaps removed (Make suitable for HXB2 alignment)
	infile = open(f, 'rU')
	fasta, left_gaps, right_gaps = convert_csf(infile)
	infile.close()
	
	# write out an unaligned FASTA file for re-mapping with bowtie2
	outfile = open(fasta_filename, 'w')
	for h, s in fasta:
		outfile.write(">{}\n{}\n".format(h, s.strip("-")))
	outfile.close()

	# Generate HXB2 sam by aligning final fasta against HXB2
	command = 'bowtie2 -f --local -N 1 -D 20 -L 15 --quiet -p 1 -x {} -U {} -S {} --no-unal --met-file {} --un {}'.format(
		hxb2_reference,fasta_filename, sam_filename,
		f.replace('.csf', '.csf.bt2_metrics'),
		f.replace('.csf', '.csf.bt2_unpaired_noalign.fastq'))
		
	timestamp("Making HXB2 sam: {}".format(command), my_rank, nprocs)
	os.system(command)
	
	# delete the temp final fasta
	os.remove(fasta_filename)

MPI.COMM_WORLD.Barrier()
if my_rank == 0: timestamp('Barrier #7 (HXB2 SAM)\n', my_rank, nprocs)



# Extract HXB2 aligned sequences from the HXB2.sam, generate amino + nuc poly files
HXB2_sams = glob(root + '/*.HXB2.sam')
for i in range(len(HXB2_sams)):
	if i % nprocs != my_rank:
		continue
	timestamp("Making nuc_poly + amino_poly from {}".format(HXB2_sams[i]), my_rank, nprocs)
	generate_nuc_counts(HXB2_sams[i],HXB2_mapping_cutoff)
	generate_amino_counts(HXB2_sams[i],HXB2_mapping_cutoff)

MPI.COMM_WORLD.Barrier()
if my_rank == 0: timestamp('Barrier #8 (Amino/nuc poly)\n', my_rank, nprocs)


# Concatenate nuc + amino poly files into aggregate summary files
if my_rank == 0:
	nuc_polys = glob(root + '/*.HXB2.nuc_poly')
	nuc_summary = root + "/HXB2.nuc_poly.summary"
	os.system("echo \"Sample,region,nuc.pos,G,A,T,C,N,-\" > {}".format(nuc_summary))
	for i in range(len(nuc_polys)): os.system("tail -n +2 {} >> {}".format(nuc_polys[i], nuc_summary))
	amino_polys = glob(root + '/*.HXB2.amino_poly')
	amino_summary = root + "/HXB2.amino_poly.summary"
	os.system("echo \"Sample,region,AA.pos,A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,*,-\" > {}".format(amino_summary))
	for i in range(len(amino_polys)): os.system("tail -n +2 {} >> {}".format(amino_polys[i], amino_summary))
MPI.COMM_WORLD.Barrier()
if my_rank == 0: timestamp('Barrier #9 (Summary files)\n', my_rank, nprocs)


# Make HXB2.amino_poly.summary_conseq from from HXB2.amino_poly.summary
if my_rank == 0:
	timestamp("Making HXB2.amino_poly.summary")
	os.system("python amino_conseq_from_summary.py {} {}".format(root + "HXB2.amino_poly.summary", amino_conseq_cutoff))


# Clean up intermediary files
if my_rank == 0:
	timestamp("Deleting intermediary files")
	files_to_delete = []
	files_to_delete += glob(root+'/.*bam')
	for file in files_to_delete:
		os.remove(file)
