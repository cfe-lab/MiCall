"""
Distribute pipeline processes across the cluster.
"""

import os,subprocess,sys
from glob import glob
from mpi4py import MPI
from miseqUtils import ambig_dict, convert_csf, convert_fasta, mixture_dict, prop_x4, sampleSheetParser, timestamp

if sys.version_info[:2] != (2, 7):
	print "MPI wrapper requires python 2.7"
	sys.exit()

my_rank = MPI.COMM_WORLD.Get_rank()
nprocs = MPI.COMM_WORLD.Get_size()

## Reference sequences
conbrefpath ='/usr/local/share/miseq/refs/cfe'
hxb2refpath='/usr/local/share/miseq/refs/'

## General parameters
read_mapping_cutoff = 10	# Read mapping cutoff
consensus_q_cutoff = 20		# Min Q for base to contribute to conseq generation (1_mapping/pileup2conseq)
REMAP_THRESHOLD = 0.95		# Fraction of fastq reads successfully mapped to be considered acceptable
MAX_REMAPS = 3			# Number of remapping attempts if below REMAP_THRESHOLD
sam2fasta_q_cutoffs = [0,10,15]	# q cutoff for base censoring (sam2fasta)
max_prop_N = 0.5		# Max proportion of censored bases allowed in a sequence (sam2fasta)

## V3 specific parameters
g2p_alignment_cutoff = 50	# Minimum alignment score during g2p scoring
g2p_fpr_cutoffs = [3.5]		# FPR cutoff for G2P X4
mincounts = [0,50,100,1000]	# Min read counts before counting in amplicon/v3 pipeline

## Arguments
root = sys.argv[1]		# Path to root of folder containing fastqs


def consensus_from_remapped_sam(root, ref, samfile, qCutoff=15):
	"""
	Take remapped sam, generate pileup, call pileup2conseq to generate 
	remap conseq, and convert to indexed .bt2 so that it can be used as 
	a reference for bowtie2.
	"""
	bamfile = samfile.replace('.sam', '.bam')
	confile = bamfile+'.pileup.conseq'
	os.system('samtools view -bt {}.fasta.fai {} > {} 2>/dev/null'.format(ref, samfile, bamfile))
	os.system('samtools sort {} {}.sorted'.format(bamfile, bamfile))
	os.system('samtools mpileup -A {}.sorted.bam > {}.pileup 2>/dev/null'.format(bamfile, bamfile))
	os.system('python pileup2conseq_v2.py {}.pileup {}'.format(bamfile, qCutoff))
	os.system('bowtie2-build -q -f {} {}'.format(confile, confile))

# Parse sample sheet
with open(root+'/SampleSheet.csv', 'rU') as ssfile:
	run_info = sampleSheetParser(ssfile)
	mode = run_info['Description']

# Print pipeline parameters to the pipeline log
if my_rank == 0:
	timestamp("PIPELINE PARAMETERS:")
	timestamp("Minimum mapping score while parsing SAM files: {}".format(read_mapping_cutoff))
	timestamp("Min Q for base to contribute to conseq generation: {}".format(consensus_q_cutoff))
	timestamp("Fraction of fastq reads successfully mapped to be considered acceptable: {}".format(REMAP_THRESHOLD))
	timestamp("Max number of remapping attempts: {}".format(MAX_REMAPS))
	timestamp("Base censoring cutoff (sam2fasta): {}".format(sam2fasta_q_cutoffs))
	timestamp("Max proportion censored bases allowed (sam2fasta): {}".format(max_prop_N))

	if mode == "Amplicon":
		timestamp("Mininum G2P alignment score: {}".format(g2p_alignment_cutoff))
		timestamp("G2P FPR cutoff for X4/R5 tropism: {}".format(g2p_fpr_cutoffs))
		timestamp("Minimum read count to contribute to proportion X4: {}".format(mincounts))

# Map and remap raw fastqs to con B (exclude files generated to record contaminants)
fastq_files = glob(root + '/*R1*.fastq')
fastq_files = [f for f in fastq_files if not f.endswith('.Tcontaminants.fastq')]

for i, fastq in enumerate(fastq_files):
	if i % nprocs != my_rank: continue
	fastq_filename = fastq.split('/')[-1]
	sample_name = fastq_filename.split('_')[0]

	if not run_info['Data'].has_key(sample_name):
		sys.stderr.write('ERROR: sample name {} (derived from fastq filename) not in SampleSheet.csv'.format(sample_name))
		sys.exit()

	is_t_primer = run_info['Data'][sample_name]['is_T_primer']

	command = "python 1_mapping.py {} {} {} {} {} {} {}".format(conbrefpath, fastq, consensus_q_cutoff, mode,
				int(is_t_primer), REMAP_THRESHOLD, MAX_REMAPS)
	timestamp(command)
	p = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	stdout, stderr = p.communicate()
	if stdout != "": sys.stdout.write("{}\n".format(stdout))
	if stderr != "": sys.stderr.write("{}\n".format(stderr))

timestamp("Arrived at barrier #1")
MPI.COMM_WORLD.Barrier()
if my_rank == 0: timestamp('All processes reached barrier #1\n')
MPI.COMM_WORLD.Barrier()

# For each remapped SAM, generate CSFs (done)
files = glob(root + '/*.remap.sam')
for i in range(len(files)):
	if i % nprocs != my_rank: continue
	filename = files[i].split('/')[-1]

	for qcut in sam2fasta_q_cutoffs:
		command = 'python 2_sam2fasta_with_base_censoring.py {} {} {} {} {}'.format(files[i], qcut, read_mapping_cutoff, mode, max_prop_N)
		timestamp(command)
		p = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		stdout, stderr = p.communicate()
		if stdout != "": sys.stdout.write("{}\n".format(stdout))
		if stderr != "": sys.stderr.write("{}\n".format(stderr))

timestamp("Arrived at barrier #2")
MPI.COMM_WORLD.Barrier()
if my_rank == 0: timestamp('All processes reached barrier #2\n')
MPI.COMM_WORLD.Barrier()

# For amplicon sequencing runs, compute g2p scores for env-mapped FASTAs
if mode == 'Amplicon':
	files = glob(root + '/*.HIV1B-env.*.fasta')
	for i, file in enumerate(files):
		if i % nprocs != my_rank: continue
		command = 'python 3_g2p_scoring.py {} {}'.format(file, g2p_alignment_cutoff)
		timestamp(command)
		p = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		stdout, stderr = p.communicate()
		if stdout != "": sys.stdout.write("{}\n".format(stdout))
		if stderr != "": sys.stderr.write("{}\n".format(stderr))

timestamp("Arrived at barrier #3")
MPI.COMM_WORLD.Barrier()



# If this is amplicon, make final proportion X4 summary, clean up files, and exit
if mode == 'Amplicon':
	if my_rank == 0:
		timestamp('Barrier #3 G2P calculations (Amplicon)\n')

		# Read each v3prot file and generate summary of v3 data in v3prot.summary
		summary_file = open(root + '/v3prot.summary', 'w')
		summary_file.write("sample,q_cutoff,fpr_cutoff,mincount, x4, reads, proportion_x4\n")
		v3_files = glob(root + '/*.v3prot')
		for i, file in enumerate(v3_files):
			prefix, gene, sam2fasta_q_cutoff = file.split('.')[:3]
			for fpr_cutoff in g2p_fpr_cutoffs:
				for mincount in mincounts:
					try:
						total_x4_count, total_count = prop_x4(file, fpr_cutoff, mincount)
						proportion_x4 = (float(total_x4_count) / float(total_count)) * 100
						sample = prefix.split('/')[-1]
						summary_file.write("{},{},{},{},{},{},{}\n".format(sample, sam2fasta_q_cutoff,
								fpr_cutoff, mincount, total_x4_count, total_count, proportion_x4))
					except:
						continue
		summary_file.close()
		timestamp('Generated v3prot.summary file - cleaning up intermediary files...\n')
MPI.COMM_WORLD.Barrier()

# Generate counts + consensus from FASTA/CSFs
files = glob(root + ('/*.fasta' if mode == 'Amplicon' else '/*.csf'))
for i, file in enumerate(files):
	if i % nprocs != my_rank:
		continue
	command = 'python 4_csf2counts.py %s %s' % (file, mode)
	timestamp(command)
	p = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	stdout, stderr = p.communicate()
	if stdout != "": sys.stdout.write("{}\n".format(stdout))
	if stderr != "": sys.stderr.write("{}\n".format(stderr))
MPI.COMM_WORLD.Barrier()

# Slice outputs + delete intermediary files
if my_rank == 0:
	command = 'python 5_slice_outputs.py {}'.format(root)
	timestamp(command)
	p = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	stdout, stderr = p.communicate()
	if stdout != "": sys.stdout.write("{}\n".format(stdout))
	if stderr != "": sys.stderr.write("{}\n".format(stderr))

	timestamp("Deleting intermediary files")
	files_to_delete = []
	files_to_delete += glob(root+'/.*bam')
	files_to_delete += glob(root + '/*.bt2')
	files_to_delete += glob(root + '/*.bt2_metrics')
	files_to_delete += glob(root + '/*.pileup')
	files_to_delete += glob(root + '/*.poly')
	for file in files_to_delete: os.remove(file)

MPI.COMM_WORLD.Barrier()
MPI.Finalize()
