"""
Distribute pipeline processes across the cluster.
"""

import os,sys
from glob import glob
from generate_hxb2_poly_files import generate_amino_counts, generate_nuc_counts
from mpi4py import MPI
from proportion_x4 import prop_x4
from sam2fasta import apply_cigar
from miseqUtils import ambig_dict, convert_csf, convert_fasta, mixture_dict, timestamp, sampleSheetParser

if sys.version_info[:2] != (2, 7):
	print "MPI wrapper requires python 2.7"
	sys.exit()

my_rank = MPI.COMM_WORLD.Get_rank()
nprocs = MPI.COMM_WORLD.Get_size()


## settings
refpath = '/usr/local/share/miseq/refs/cfe'	# Consensus B refs
hxb2refpath="/usr/local/share/miseq/refs/" 	# location of HXB2 reference files
HXB2_mapping_cutoff = 10			# Read mapping quality cutoff
amino_conseq_cutoff = 10			# Minimum occurences of char before it counts for conseq generation
fpr_cutoffs = [3.5]				# Cutoffs for assuming X4-ness
mincounts = [3]					# Minimum occurences of char before it counts in amplicon/v3
g2p_alignment_cutoff = 50			# Minimum score of alignment during g2p scoring


## arguments
root = sys.argv[1]				# Path to folder containing fastqs
qCutoff = int(sys.argv[2])			# INCORRECT - FIXME


# parse sample sheet
ssfile = open(root+'/SampleSheet.csv', 'rU')
run_info = sampleSheetParser(ssfile)
ssfile.close()

mode = run_info['Description'] # Nextera or Amplicon


# prepare logs
log_file = open(root + "/pipeline_output.MPI_rank_{}.log".format(my_rank), "w")
sys.stderr = log_file



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


# exclude files that were generated to record contaminants
files = [f for f in files if not f.endswith('.Tcontaminants.fastq')]


for i, file in enumerate(files):
	if i % nprocs != my_rank:
		continue
	filename = file.split('/')[-1]
	
	sample_name = filename.split('_')[0]
	assert run_info['Data'].has_key(sample_name), 'ERROR: sample name %s not in SampleSheet.csv' % sample_name
	is_t_primer = run_info['Data'][sample_name]['is_T_primer']
	
	timestamp('1_mapping {} {} {}'.format(refpath, file, qCutoff), my_rank, nprocs)
	os.system("python 1_mapping.py {} {} {} {}".format(refpath, file, qCutoff, int(is_t_primer)))

MPI.COMM_WORLD.Barrier()
if my_rank == 0: timestamp('Barrier #1 (Prelim clade B + sample specific remapping)\n', my_rank, nprocs)
MPI.COMM_WORLD.Barrier()


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
MPI.COMM_WORLD.Barrier()


# For amplicon sequencing runs, compute g2p scores for env-mapped FASTAs
if mode == 'Amplicon':
	files = glob(root + '/*.HIV1B-env.*.fasta')
	for i, file in enumerate(files):
		if i % nprocs != my_rank:
			continue
		command = 'python 3_g2p_scoring.py {} {}'.format(file, g2p_alignment_cutoff)
		timestamp(command, my_rank, nprocs)
		os.system(command)

MPI.COMM_WORLD.Barrier()




# If this is amplicon, make final proportion X4 summary, clean up files, and exit
if mode == 'Amplicon':
	if my_rank == 0:
		timestamp('Barrier #3 G2P calculations (Amplicon)\n', my_rank, nprocs)

		# Read each v3prot file and generate summary of v3 data in v3prot.summary
		summary_file = open(root + '/v3prot.summary', 'w')
		summary_file.write("sample,qCutoff,fpr_cutoff,mincount,proportion_x4\n")
		v3_files = glob(root + '/*.v3prot')
		for i, file in enumerate(v3_files):
			prefix, gene, qCutoff = file.split('.')[:3]
			for fpr_cutoff in fpr_cutoffs:
				for mincount in mincounts:
					try:
						total_x4_count, total_count = prop_x4(file, fpr_cutoff, mincount)
						proportion_x4 = (float(total_x4_count) / float(total_count)) * 100
						sample = prefix.split('/')[-1]
						summary_file.write("{},{},{},{},{}\n".format(sample, qCutoff, fpr_cutoff, mincount, proportion_x4))
					except:
						continue
		summary_file.close()
		timestamp('Generated v3prot.summary file - cleaning up intermediary files...\n', my_rank, nprocs)



# generate counts and consensus from FASTA files
files = glob(root + '/*.fasta' if mode == 'Amplicon' else '/*.csf')
for i, file in enumerate(files):
	if i % nprocs != my_rank:
		continue
	command = 'python 4_csf2counts.py %s %s' % (file, mode)
	timestamp(command, my_rank, nprocs)
	os.system(command)


# Clean up intermediary files
if my_rank == 0:
	timestamp("Deleting intermediary files")
	files_to_delete = []
	files_to_delete += glob(root+'/.*bam')
	files_to_delete += glob(root + '/*.bt2')
	files_to_delete += glob(root + '/*.bt2_metrics')
	#files_to_delete += glob(root + '/*.fastq') # why delete original file?
	files_to_delete += glob(root + '/*.pileup')
	files_to_delete += glob(root + '/*.poly')
	for file in files_to_delete: os.remove(file)

log_file.close()
MPI.Finalize()
