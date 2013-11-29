import logging, miseq_logging, os, subprocess, sys
from glob import glob
from miseqUtils import ambig_dict, convert_csf, convert_fasta, mixture_dict, prop_x4, sampleSheetParser, timestamp
from mpi4py import MPI
from miseq_modules import remap, mapping, sam2fasta_with_base_censoring, g2p_scoring, csf2counts, slice_outputs

my_rank = MPI.COMM_WORLD.Get_rank()
nprocs = MPI.COMM_WORLD.Get_size()

## Reference sequences
conbrefpath ='/usr/local/share/miseq/refs/cfe'
hxb2refpath='/usr/local/share/miseq/refs/'

## General parameters
read_mapping_cutoff = 10	# Minimum read mapping quality
consensus_q_cutoff = 20		# Min Q for base to contribute to conseq generation (1_mapping/pileup2conseq)
REMAP_THRESHOLD = 0.95		# Fraction of fastq reads successfully mapped to be considered acceptable
MAX_REMAPS = 3			# Number of remapping attempts if below REMAP_THRESHOLD
sam2fasta_q_cutoffs = [0,10,15]	# q cutoff for base censoring (sam2fasta)
max_prop_N = 0.5		# Max proportion of censored bases allowed in a sequence (sam2fasta)
consensus_mixture_cutoffs = [0.01, 0.02, 0.05, 0.1, 0.2, 0.25]	# Determines thresholds for mixture characters in the consensus

## V3 specific parameters
g2p_alignment_cutoff = 50	# Minimum alignment score during g2p scoring
g2p_fpr_cutoffs = [3.5]		# FPR cutoff for G2P X4
mincounts = [0,50,100,1000]	# Min read counts before counting in amplicon/v3 pipeline

## Arguments
root = sys.argv[1]		# Path on local cluster containing fastqs/SampleSheet.csv

# Setup MPI process-specific logging
log_file = "{}/rank_{}.log".format(root, my_rank)
logger = miseq_logging.init_logging(log_file, file_log_level=logging.DEBUG, console_log_level=logging.INFO)

# Parse sample sheet
with open(root+'/SampleSheet.csv', 'rU') as ssfile:
	logger.debug("sampleSheetParser({})".format(ssfile))
	run_info = sampleSheetParser(ssfile)
	mode = run_info['Description']

# Print pipeline parameters to the pipeline log
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
	logger.info("===== End of pipeline parameters =====")

logger.debug("Arrived at barrier #1")
MPI.COMM_WORLD.Barrier()

# Map and remap raw fastqs to con B (exclude files generated to record contaminants)
fastq_files = glob(root + '/*R1*.fastq')
fastq_files = [f for f in fastq_files if not f.endswith('.Tcontaminants.fastq')]

for i, fastq in enumerate(fastq_files):
	if i % nprocs != my_rank: continue
	fastq_filename = fastq.split('/')[-1]
	sample_name = fastq_filename.split('_')[0]

	if not run_info['Data'].has_key(sample_name):
		logger.error('ERROR: sample name {} (derived from fastq filename) not in SampleSheet.csv'.format(sample_name))
		sys.exit()

	is_t_primer = run_info['Data'][sample_name]['is_T_primer']
	logging.info("mapping({}, {}, {}, {}, {}, {}, {})".format(conbrefpath, fastq, consensus_q_cutoff, mode, int(is_t_primer), REMAP_THRESHOLD, MAX_REMAPS))
	mapping(conbrefpath, fastq, consensus_q_cutoff, mode, int(is_t_primer), REMAP_THRESHOLD, MAX_REMAPS)

logger.debug("Arrived at barrier #2")
MPI.COMM_WORLD.Barrier()
if my_rank == 0: logger.debug('All processes reached barrier #2')
MPI.COMM_WORLD.Barrier()

# For each remapped SAM, generate CSFs (done)
files = glob(root + '/*.remap.sam')
for i in range(len(files)):
	if i % nprocs != my_rank: continue
	filename = files[i].split('/')[-1]

	for qcut in sam2fasta_q_cutoffs:
		logger.info("sam2fasta_with_base_censoring({}, {}, {}, {}, {})".format(files[i], qcut, read_mapping_cutoff, mode, max_prop_N))
		sam2fasta_with_base_censoring(files[i], qcut, read_mapping_cutoff, mode, max_prop_N)

logger.debug("Arrived at barrier #3")
MPI.COMM_WORLD.Barrier()
if my_rank == 0: logger.debug('All processes reached barrier #3')
MPI.COMM_WORLD.Barrier()

# For amplicon sequencing runs, compute g2p scores for env-mapped FASTAs
if mode == 'Amplicon':
	files = glob(root + '/*.HIV1B-env.*.fasta')
	for i, file in enumerate(files):
		if i % nprocs != my_rank: continue
		logger.info("g2p_scoring({},{})".format(file, g2p_alignment_cutoff))
		g2p_scoring(file, g2p_alignment_cutoff)
logging.debug("Arrived at barrier #4")
MPI.COMM_WORLD.Barrier()

# If this is amplicon, make final proportion X4 summary, clean up files, and exit
if mode == 'Amplicon' and my_rank == 0:

	# Read each v3prot file and generate summary of v3 data in v3prot.summary
	with open(root + '/v3prot.summary', 'w') as summary_file:
		summary_file.write("sample,q_cutoff,fpr_cutoff,mincount, x4, reads, proportion_x4\n")
		v3_files = glob(root + '/*.v3prot')
		for i, file in enumerate(v3_files):
			prefix, gene, sam2fasta_q_cutoff = file.split('.')[:3]
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

	logging.info('Generated v3prot.summary file - cleaning up intermediary files...')
logging.debug("Arrived at barrier #5")
MPI.COMM_WORLD.Barrier()

# Generate counts + consensus from FASTA/CSFs
files = glob(root + ('/*.fasta' if mode == 'Amplicon' else '/*.csf'))
for i, file in enumerate(files):
	if i % nprocs != my_rank:
		continue
	logging.info("csf2counts({},{})".format(file,mode))
	csf2counts(file,mode,consensus_mixture_cutoffs)
logging.debug("Arrived at barrier #6")
MPI.COMM_WORLD.Barrier()

# Slice outputs + delete intermediary files
if my_rank == 0:
	logging.info("slice_outputs({})".format(root))
	slice_outputs(root)
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
