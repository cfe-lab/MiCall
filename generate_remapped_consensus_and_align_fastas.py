"""
INPUT: Folder path containing remapped sams
OUTPUT: *.final.sam, *.HXB2.sam, *.amino_poly, amino_poly.summary
"""
import os,sys
from glob import glob
from seqUtils import ambig_dict, convert_csf, mixture_dict, translate_nuc
from sam2fasta import apply_cigar
from mpi4py import MPI

my_rank = MPI.COMM_WORLD.Get_rank()
nprocs = MPI.COMM_WORLD.Get_size()

amino_alphabet = 'ACDEFGHIKLMNPQRSTVWY*-'
conB_refpath="/usr/local/share/miseq/refs/cfe"
hxb2refpath="/usr/local/share/miseq/refs/"      # hxb2_pol_ref.fasta
root = sys.argv[1]

# Set debug to 1 and modify debug parameters to test development on this pipeline
debug = 0
debug_root = "/data/miseq/0_testing"
debug_seq = "35759A"

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

# Step 1A: From remapped sams, generate remap conseqs
remapped_sams = glob(root + '/*.remap.sam')
if debug == 1: remapped_sams = glob('{}/{}-PR-RT.HIV1B-pol.remap.sam'.format(debug_root,debug_seq))

for i in range(len(remapped_sams)):
	if i % nprocs != my_rank: continue
	remapped_sam = remapped_sams[i]
	remap_conseq = remapped_sam.replace('.sam', '.bam.pileup.conseq')
	consensus_from_remapped_sam(root,remap_conseq,remapped_sam,20)

MPI.COMM_WORLD.Barrier()

# Step 1B: From remap csfs, create remap fastas to be aligned against remap conseq
csf_files = glob(root + '/*.csf')
if debug == 1: csf_files = glob('{}/{}-PR-RT.HIV1B-pol.20.csf'.format(debug_root,debug_seq))

for i in range(len(csf_files)):
	if i % nprocs != my_rank: continue
	f = csf_files[i]

	# convert_csf returns left-padded sequences
	fasta_filename = f.replace('.csf', '.csf.fasta')
	infile = open(f, 'rU')
	fasta, left_gap_position, right_gap_position = convert_csf(infile)
	infile.close()

	# Write remap temp fasta: remove gaps so fasta is ready for alignment with remap conseq
	print "Making (temp) remap fasta from {}".format(f)
	outfile = open(fasta_filename, 'w')
	for j, (h, s) in enumerate(fasta): outfile.write(">{}\n{}\n".format(h, s.strip("-")))
	outfile.close()

	prefix, gene = fasta_filename.split('.')[:2]
	conseq_filename = "{}.{}.remap.bam.pileup.conseq".format(prefix,gene)
	sam_filename = f.replace('.csf', '.final.sam')

	# Step 1C: Generate final sam by aligning fastas against remap conseq, then delete the temp fasta
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


# Step 2: From final sam, extract final CSFs
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

# Step 3: HXB2 alignment of final csf
final_csfs = glob(root + '/*.final.sam.20.csf')
if debug == 1: final_sams = glob('{}/{}-PR-RT.HIV1B-pol.20.final.sam.20.csf'.format(debug_root,debug_seq))

for i in range(len(final_csfs)):
	if i % nprocs != my_rank: continue
	f = final_csfs[i]

	sample, region, qCut = f.split('.')[:3]
	fasta_filename = f.replace('.csf', '.fasta')
	sam_filename = f.replace('.csf', '.HXB2.sam')
	hxb2_reference = "{}{}_hxb2.fasta".format(hxb2refpath,region)   # /usr/local/share/miseq/refs/HIV1B-pol_hxb2.fasta

	# Step 3A: From final csf, generate temp final fasta for HXB2 alignment
	print "Making (temp) final fasta from {}".format(f)
	infile = open(f, 'rU')
	fasta, left_gap_position, right_gap_position = convert_csf(infile)
	infile.close()
	outfile = open(fasta_filename, 'w')
	for j, (h, s) in enumerate(fasta): outfile.write(">{}\n{}\n".format(h, s.strip("-")))
	outfile.close()

	# Step 3B: Generate HXB2 sam by aligning the final fasta against HXB2, then delete the temp fasta
	command = 'bowtie2 -f --quiet -p 1 -x {} -U {} -S {} --no-unal --met-file {} --un {}'.format(
		hxb2_reference,fasta_filename, sam_filename,
		f.replace('.csf', '.csf.bt2_metrics'),
		f.replace('.csf', '.csf.bt2_unpaired_noalign.fastq'))
        print "\nMaking HXB2 sam: {}".format(command)
	os.system(command)
	os.remove(fasta_filename)
MPI.COMM_WORLD.Barrier()

# Step 4A: Extract HXB2 aligned sequences from the HXB2.sam and generate poly counts
HXB2_sams = glob(root + '/*.HXB2.sam')
if debug == 1: HXB2_sams = glob('{}/{}-PR-RT.HIV1B-pol.20.final.sam.20.HXB2.sam'.format(debug_root,debug_seq))

for i in range(len(HXB2_sams)):
	if i % nprocs != my_rank: continue
	f = HXB2_sams[i]
	prefix, region, qCut = f.split('.')[:3]
	sample = prefix.split("/")[-1]
	infile = open(f, 'rU')
	lines = infile.readlines()

	# Determine where the SAM header ends
	for start, line in enumerate(lines):
		if not line.startswith('@'): break

	# For each read (Each line in the SAM), generate the amino sequence and tally amino counts
	aminos = {}
	j = start
	while j < len(lines):
		qname, flag, refname, pos, mapq, cigar, rnext, pnext, tlen, seq, qual = lines[j].strip('\n').split('\t')[:11]
		shift, seq_aligned, qual_aligned = apply_cigar(cigar, seq, qual)

		# Check the bitwise flag and filter out unmapped sequences
		if int(flag) & 4 == 4: continue

		# Determine ORF offset, translate on the correct ORF, add a left offset to the amino seq
		ORF = (int(pos)-1) % 3
		p = translate_nuc(seq_aligned, ORF)
		amino_offset = (int(pos)-1)/3
		p = "{}{}".format("?"*amino_offset,p)

		# Iterate through the amino seq and tally the (coordinate,aminos) counts
		for coord, aa in enumerate(p):
			if aa in '?': continue
			if not aminos.has_key(coord): aminos.update({coord:{}})
			if not aminos[coord].has_key(aa): aminos[coord].update({aa: 0})
			aminos[coord][aa] += 1

		# Move to next SAM line
		j += 1

	# For each sample, write amino frequency data to a .amino_poly file
	print "Making amino matrix: {}".format(f.replace('.HXB2.sam','.HXB2.amino_poly'))
	amino_poly = open(f.replace('.HXB2.sam','.HXB2.amino_poly'), 'w')
	keys = aminos.keys()
	keys.sort()
	amino_poly.write('Sample,region,AA.pos,' + ','.join(amino_alphabet) + '\n')
	for k in keys:
		amino_poly.write("{},{},{}".format(sample,region,str(k+1)))
		for aa in amino_alphabet: amino_poly.write(",{}".format(aminos[k].get(aa,0)))
		amino_poly.write('\n')
	amino_poly.close()

	# Close the SAM
	infile.close()
MPI.COMM_WORLD.Barrier()

# Step 4B: Generate an amino summary file
if my_rank == 0:
	amino_polys = glob(root + '*.HXB2.amino_poly')
	amino_summary = root + "HXB2.amino_poly.summary"
	print "Dumping amino matrix data into {}".format(amino_summary)
	for i in range(len(amino_polys)): os.system("tail -n +2 {} >> {}".format(amino_polys[i], amino_summary))
MPI.COMM_WORLD.Barrier()

MPI.Finalize()
