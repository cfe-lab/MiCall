"""
1_MPI_Wrapper essentially performs the full MiSeq pipeline - both the general Amplicon and Nextera,
along with the V3-specific HIV co-receptor tropism prediction.

First, the sample sheet is parsed to determine if this is an Amplicon or Nextera run. Then, fastq
reads are mapped against HIV/HCV/HLA genes to determine what regions these reads come from.

Frequency by coordinate data are then generated in the context of a standard coordinate system: for
example, HXB2 in HIV and H77 for HCV.

HIV1B-env sequences are also subject to a viral tropism scoring using the G2P algorithm.
"""

import logging, miseq_logging, os, subprocess, sys
from glob import glob
from miseqUtils import ambig_dict, convert_csf, convert_fasta, mixture_dict, prop_x4, sampleSheetParser
from mpi4py import MPI
from miseq_modules import csf2counts, g2p_scoring, mapping, remap, sam2fasta_with_base_censoring, slice_outputs

my_rank = MPI.COMM_WORLD.Get_rank()
nprocs = MPI.COMM_WORLD.Get_size()

## Reference sequences
mapping_ref_path = '/usr/local/share/miseq/refs/cfe'

## General parameters
read_mapping_cutoff = 10		# Minimum read mapping quality
consensus_q_cutoff = 20			# Min Q for base to contribute to conseq generation (1_mapping/pileup2conseq)
REMAP_THRESHOLD = 0.95			# Fraction of fastq reads successfully mapped to be considered acceptable
MAX_REMAPS = 3				# Number of remapping attempts if below REMAP_THRESHOLD
sam2fasta_q_cutoffs = [0,10,15]		# q cutoff for base censoring (sam2fasta)
max_prop_N = 0.5			# Max proportion of censored bases allowed in a sequence (sam2fasta)
consensus_mixture_cutoffs = [0.01,0.02,0.05,0.1,0.2,0.25]	# Proportion thresholds for mixtures in the consensus

## V3 specific parameters
g2p_alignment_cutoff = 50		# Minimum alignment score during g2p scoring
g2p_fpr_cutoffs = [3.0,3.5,4.0,5.0]	# FPR cutoff for G2P X4
mincounts = [0,50,100,1000]		# Min read counts before counting in amplicon/v3 pipeline

# Define slices with (Label, region to slice, start, end)
# Coordinates are in nucleotide space: start/end are inclusive (Relative to HXB2 aligned sequences)
region_slices = [("PROTEASE", "HIV1B-pol", 1, 297), ("PRRT", "HIV1B-pol", 1, 1617), ("INTEGRASE", "HIV1B-pol", 1978, 2844),
		 ("P17", "HIV1B-gag", 1, 396), ("P24", "HIV1B-gag", 397, 1089), ("P2P7P1P6","HIV1B-gag", 1090, 1502),
		 ("GP120","HIV1B-env", 1, 1533), ("GP41", "HIV1B-env", 1534, 2570), ("V3", "HIV1B-env", 887, 993)]

## Arguments
root = sys.argv[1]			# Path on local cluster containing fastqs/SampleSheet.csv

# Each MPI process will have a separate log file which contains highly detailed (DEBUG-level) information
log_file = "{}/rank_{}.log".format(root, my_rank)
logger = miseq_logging.init_logging(log_file, file_log_level=logging.DEBUG, console_log_level=logging.INFO)

# Parse sample sheet
with open(root+'/SampleSheet.csv', 'rU') as ssfile:
	logger.debug("sampleSheetParser({})".format(ssfile))
	run_info = sampleSheetParser(ssfile)
	mode = run_info['Description']

# Print pipeline parameters to the log
if my_rank == 0:
	logger.info("===== Pipeline parameters =====")
	logger.info("Minimum read mapping score while parsing SAM files: {}".format(read_mapping_cutoff))
	logger.info("Min Q for base to contribute to conseq generation: {}".format(consensus_q_cutoff))
	logger.info("Fraction of fastq reads successfully mapped to be considered acceptable: {}".format(REMAP_THRESHOLD))
	logger.info("Max number of remapping attempts: {}".format(MAX_REMAPS))
	logger.info("Base censoring cutoff (sam2fasta): {}".format(sam2fasta_q_cutoffs))
	logger.info("Max proportion censored bases allowed (sam2fasta): {}".format(max_prop_N))
	logger.info("Base proportion needed to count towards a consensus sequence mixture: {}".format(consensus_mixture_cutoffs))
	if mode == "Amplicon":
		logger.info("Mininum G2P alignment score: {}".format(g2p_alignment_cutoff))
		logger.info("G2P FPR cutoff for X4/R5 tropism: {}".format(g2p_fpr_cutoffs))
		logger.info("Minimum read count to contribute to proportion X4: {}".format(mincounts))
	logger.info("Region slices: {}".format(region_slices))
	logger.info("===== End of pipeline parameters =====")

logger.debug("Arrived at barrier #1")
MPI.COMM_WORLD.Barrier()

# Map and remap raw fastqs to con B (exclude files generated to record contaminants)
fastq_files = glob(root + '/*R1*.fastq')
fastq_files = [f for f in fastq_files if not f.endswith('.Tcontaminants.fastq')]

# Align reads in fastq files to generate sam files
for i, fastq in enumerate(fastq_files):
	if i % nprocs != my_rank: continue
	fastq_filename = fastq.split('/')[-1]
	sample_name = fastq_filename.split('_')[0]

	if not run_info['Data'].has_key(sample_name):
		logger.error('ERROR: sample name {} (derived from fastq filename) not in SampleSheet.csv'.format(sample_name))
		sys.exit()

	is_t_primer = run_info['Data'][sample_name]['is_T_primer']

	# Map fastq reads to a static reference, generate sample specific consensus, remap reads to specific consensus
	logging.info("mapping({}, {}, {}, {}, {}, {}, {})".format(mapping_ref_path, fastq, consensus_q_cutoff, mode, is_t_primer, REMAP_THRESHOLD, MAX_REMAPS))
	mapping(mapping_ref_path, fastq, consensus_q_cutoff, mode, is_t_primer, REMAP_THRESHOLD, MAX_REMAPS)

logger.debug("Arrived at barrier #2")
MPI.COMM_WORLD.Barrier()

# For each remapped SAM
files = glob(root + '/*.remap.sam')
for i in range(len(files)):
	if i % nprocs != my_rank: continue
	filename = files[i].split('/')[-1]

	# Generate a fasta or csf file under different q cutoff base censoring rules (Store the q cutoff used in the filename)
	for qcut in sam2fasta_q_cutoffs:
		logger.info("sam2fasta_with_base_censoring({}, {}, {}, {}, {})".format(files[i], qcut, read_mapping_cutoff, mode, max_prop_N))
		sam2fasta_with_base_censoring(files[i], qcut, read_mapping_cutoff, mode, max_prop_N)

logger.debug("Arrived at barrier #3")
MPI.COMM_WORLD.Barrier()

# For amplicon runs, compute V3 g2p scores from HIV1B-env fasta files
if mode == 'Amplicon':
	files = glob(root + '/*.HIV1B-env.*.fasta')
	for i, file in enumerate(files):
		if i % nprocs != my_rank: continue

		# For each HIV1B-env fasta, generate a v3prot file
		logger.info("g2p_scoring({},{})".format(file, g2p_alignment_cutoff))
		g2p_scoring(file, g2p_alignment_cutoff)

logging.debug("Arrived at barrier #4")
MPI.COMM_WORLD.Barrier()

# For amplicon runs, determine the proportion X4 of each sample
if mode == 'Amplicon' and my_rank == 0:

	# For each v3prot file (Denoting the base censoring q cutoff used in the filename)
	with open(root + '/v3prot.summary', 'w') as summary_file:
		summary_file.write("sample,q_cutoff,fpr_cutoff,mincount, x4, reads, proportion_x4\n")
		v3_files = glob(root + '/*.v3prot')
		for i, file in enumerate(v3_files):
			prefix, gene, sam2fasta_q_cutoff = file.split('.')[:3]

			# Determine the proportion x4 under different FPR cutoffs and min counts
			for fpr_cutoff in g2p_fpr_cutoffs:
				for mincount in mincounts:
					try:
						total_x4_count, total_count = prop_x4(file, fpr_cutoff, mincount)
						logging.debug("{}, {} = prop_x4({},{},{})".format(total_x4_count, total_count, file, fpr_cutoff, mincount))
						proportion_x4 = (float(total_x4_count) / float(total_count)) * 100
						sample = prefix.split('/')[-1]
						summary_file.write("{},{},{},{},{},{},{}\n".format(sample, sam2fasta_q_cutoff,
								fpr_cutoff, mincount, total_x4_count, total_count, proportion_x4))
					except:
						continue

	logging.info('Generated v3prot.summary file')
logging.debug("Arrived at barrier #5")
MPI.COMM_WORLD.Barrier()

# From each fasta/csf file
files = glob(root + ('/*.fasta' if mode == 'Amplicon' else '/*.csf'))

# Determine the nucleotide/amino counts, along with the consensus, in HXB2/H77 space
for i, file in enumerate(files):
	if i % nprocs != my_rank: continue
	logging.info("csf2counts({},{})".format(file,mode))
	csf2counts(file,mode,consensus_mixture_cutoffs)
logging.debug("Arrived at barrier #6")
MPI.COMM_WORLD.Barrier()

# Represent count output with respect to different region codes (Ex: V3 in env)
if my_rank == 0:
	logging.info("slice_outputs({})".format(root))
	slice_outputs(root, region_slices)

# Delete files on the local cluster that don't need to be kept
if my_rank == 0:
	logging.info("Deleting intermediary files")
	files_to_delete = []
	files_to_delete += glob(root+'/.*bam')
	files_to_delete += glob(root + '/*.bt2')
	files_to_delete += glob(root + '/*.bt2_metrics')
	files_to_delete += glob(root + '/*.pileup')
	files_to_delete += glob(root + '/*.poly')
	for file in files_to_delete:
		logging.debug("os.remove({})".format(file))
		os.remove(file)

logging.debug("Arrived at barrier #7")
MPI.COMM_WORLD.Barrier()
MPI.Finalize()
