import logging, miseq_logging, miseq_modules, miseqUtils, os, subprocess, sys, time
sys.path.append('/usr/local/share/fifo_scheduler')
from fifo_scheduler import Factory
from glob import glob

# Input parameters
root = "/data/miseq/130711_M01841_0010_000000000-A3TCY"		# root = sys.argv[1]
mode = "Amplicon"						# mode = sys.argv[2]
is_t_primer = False						# is_t_primer = sys.argv[3]

# Logging parameters
log_file = "{}/pipeline_output.log".format(root)
logger = miseq_logging.init_logging(log_file, file_log_level=logging.DEBUG, console_log_level=logging.INFO)

# Mapping parameters
mapping_ref_path = "/usr/local/share/miseq/refs/cfe"
bowtie_threads = 8                      # Bowtie performance roughly scales with number of threads
min_mapping_efficiency = 0.95		# Fraction of fastq reads mapped needed
max_remaps = 3				# Number of remapping attempts if mapping efficiency unsatisfied
consensus_q_cutoff = 20                 # Min Q for base to contribute to conseq (pileup2conseq)
mapping_factory_resources = [("bpsh -1", 1), ("bpsh 0", 1), ("bpsh 1", 2), ("bpsh 2", 2)]

# sam2csf parameters
sam2csf_q_cutoffs = [0,10,15]		# Q-cutoff for base censoring
max_prop_N = 0.5			# Drop reads with more censored bases than this proportion
read_mapping_cutoff = 10		# Minimum bowtie read mapping quality

# g2p parameters (Amplicon only)
g2p_alignment_cutoff = 50		# Minimum alignment score during g2p scoring
g2p_fpr_cutoffs = [3.0,3.5,4.0,5.0]	# FPR cutoff to determine R5/X4 tropism
v3_mincounts = [0,50,100,1000]		# Min number of reads to contribute to %X4 calculation

# csf2counts parameters
consensus_mixture_cutoffs = "0.01,0.02,0.05,0.1,0.2,0.25"
final_alignment_ref_path = "/usr/local/share/miseq/development/miseqpipeline/csf2counts_amino_sequences.csv"

# File extensions to delete at the end of the run
file_extensions_to_delete = ['bam', 'bt2', 'bt2_metrics', 'fastq', 'pileup', 'pileup.conseq', 'poly']

def factory_barrier(my_factory):
	"""Wait for factory to complete all queued jobs"""
	while True:
		processes_spawned = my_factory.supervise()
		if processes_spawned:
			for p, command in processes_spawned:
				logger.info("{}".format(command))
		if my_factory.completely_idle(): break
		time.sleep(1)
	return

# Mapping factory suitable for multi-thread jobs (4 processes * 8 threads / job = 32 cores allocated)
mapping_factory = Factory(mapping_factory_resources)

### Begin Mapping
fastq_files = glob(root + '/*R1*.fastq')
fastq_files = [f for f in fastq_files if not f.endswith('.Tcontaminants.fastq')]
for fastq in fastq_files:
	fastq_filename = os.path.basename(fastq)
	sample_name = fastq_filename.split('_')[0]
	command = "python STEP_1_MAPPING.py {} {} {} {} {} {} {} {}".format(mapping_ref_path,
			fastq, consensus_q_cutoff, mode, is_t_primer, min_mapping_efficiency, max_remaps, bowtie_threads)
	queue_request = mapping_factory.queue_work(command, log_file, log_file)
	if queue_request:
		p, command = queue_request
		logger.info("{}".format(command))
factory_barrier(mapping_factory)

# Make factory more suitable for single thread jobs (64 cores allocated)
single_thread_resources = [("bpsh -1", 16), ("bpsh 0", 16), ("bpsh 1", 16), ("bpsh 2", 16)] 
single_thread_factory = Factory(single_thread_resources)

### Begin sam2csf
for file in glob(root + '/*.remap.sam'):
	filename = file.split('/')[-1]

	# Generate csf with different q cutoff censoring rules
	for qcut in sam2csf_q_cutoffs:

		# CSFs are sorted by read prevalence for Amplicon and left-offset for Nextera
		command = "python STEP_2_SAM2CSF.py {} {} {} {} {}".format(file, qcut, read_mapping_cutoff, mode, max_prop_N)
		queue_request = single_thread_factory.queue_work(command, log_file, log_file)
		if queue_request:
			p, command = queue_request
			logger.info("{}".format(command))
factory_barrier(single_thread_factory)

### Begin g2p (For Amplicon)
if mode == 'Amplicon':

	# Compute g2p V3 tropism scores from HIV1B-env csf files and store in v3prot files
	for env_csf_file in glob(root + '/*.HIV1B-env.*.csf'):
		command = "python STEP_3_G2P.py {} {}".format(env_csf_file, g2p_alignment_cutoff)
		queue_request = single_thread_factory.queue_work(command, log_file, log_file)
		if queue_request:
			p, command = queue_request
			logger.info("Starting: {}".format(command))
	factory_barrier(single_thread_factory)

	# Summarize the v3prot files in v3_tropism_summary.txt
	with open("{}/v3_tropism_summary.txt".format(root), 'wb') as summary_file:
		summary_file.write("sample,q_cutoff,fpr_cutoff,min_count,total_x4,total_reads,proportion_x4\n")
		for file in glob(root + '/*.v3prot'):
			prefix, gene, sam2csf_q_cutoff = file.split('.')[:3]

			# Determine proportion x4 under different FPR cutoffs / min counts
			for fpr_cutoff in g2p_fpr_cutoffs:
				for mincount in v3_mincounts:
					try:
						sample = prefix.split('/')[-1]
						proportion_x4, total_x4_count, total_count = miseqUtils.prop_x4(file, fpr_cutoff, mincount)
						summary_file.write("{},{},{},{},{},{},{:.3}\n".format(sample, sam2csf_q_cutoff,
								fpr_cutoff, mincount, total_x4_count, total_count, proportion_x4))
					except:
						continue

### Begin csf2counts
for csf_file in glob(root + '/*.csf'):

	# Determine nucleotide/amino counts, along with the consensus, in HXB2/H77 space
	command = "python STEP_4_CSF2COUNTS.py {} {} {} {}".format(csf_file,mode,consensus_mixture_cutoffs,final_alignment_ref_path)
	queue_request = single_thread_factory.queue_work(command, log_file, log_file)
	if queue_request:
		p, command = queue_request
		logger.info("{}".format(command))
factory_barrier(single_thread_factory)

# Remove empty files generated by csf2counts
files = glob(root + '/*.csv')
files += glob (root + '/*.conseq')
for file in files:
	stdout, stderr = subprocess.Popen(['wc', '-l', file], stdout = subprocess.PIPE, stderr = subprocess.PIPE).communicate()
	num_lines = int(stdout.split()[0])
	if num_lines == 0:
		logging.debug("Removing empty file {}".format(file))
		os.remove(file)

# Represent count output with respect to different region codes (Ex: V3 in env)
#logging.info("slice_outputs({})".format(root))
#miseq_modules.slice_outputs(root, region_slices)

# Delete files on the local cluster that don't need to be kept
logging.info("Deleting intermediary files")
for extension in file_extensions_to_delete:
	for file in glob("{}/*.{}".format(root, extension)):
		logging.debug("os.remove({})".format(file))
		os.remove(file)
