"""
To make pipeline portable, allow user to specify local paths and thread counts.
"""

### Core settings ###
pipeline_version = '7.0'        # Change for each release
are_temp_folders_deleted = True # Should FIFO worker clean up working folders?

## Mapping parameters
bowtie_version = '2.2.1'        # version of bowtie2, used for version control
bowtie_path = 'bowtie2-align-s'         # path to executable, so you can install more than one version
bowtie_build_path = 'bowtie2-build-s'
bowtie_threads = 4              # Bowtie performance roughly scales with number of threads
samtools_version = None #'1.1'
samtools_path = None #'samtools-1.1'
min_mapping_efficiency = 0.95   # Fraction of fastq reads mapped needed
max_remaps = 3                  # Number of remapping attempts if mapping efficiency unsatisfied
consensus_q_cutoff = 20         # Min Q for base to contribute to conseq (pileup2conseq)

## sam2aln parameters
sam2aln_q_cutoffs = [15]  # Q-cutoff for base censoring
max_prop_N = 0.5                 # Drop reads with more censored bases than this proportion
read_mapping_cutoff = 10         # Minimum bowtie read mapping quality

## aln2counts parameters
amino_alphabet = 'ACDEFGHIKLMNPQRSTVWY*'
conseq_mixture_cutoffs = [0.01, 0.02, 0.05, 0.1, 0.2, 0.25]

# Read and reference gap open/extension penalties.
read_gap_open_prelim = 10
read_gap_extend_prelim = 3
ref_gap_open_prelim = 10
ref_gap_extend_prelim = 3
read_gap_open_remap = read_gap_open_prelim
read_gap_extend_remap = read_gap_extend_prelim
ref_gap_open_remap = ref_gap_open_prelim
ref_gap_extend_remap = ref_gap_extend_prelim
gap_open_coord = 40
gap_extend_coord = 10


### Monitor settings ###
instrument_number = 'M01841'  # for Illumina MiSeq, second item in run folder name
production = False  # set this to True to push results to NAS

# Local path for writing data
home = '/data/miseq/'
nruns_to_store = 20  # protect X most recent runs from cleaning up intermediate files

# Scheduling processes: these should be a multiple of the total number of slots
# in your hostfile.
mapping_processes = 36
counting_processes = 144

rawdata_mount = '/media/RAW_DATA/'  # NAS
delay = 3600  # Delay (seconds) for polling NAS for unprocessed runs
NEEDS_PROCESSING = 'needsprocessing'  # File flags
ERROR_PROCESSING = 'errorprocessing'
DONE_PROCESSING = 'doneprocessing'
QC_UPLOADED = 'qc_uploaded'

# Connection to QAI RESTful API for uploading results
qai_user = "FILLINUSERNAME"
qai_password = "****"
qai_path = "http://192.168.X.Y:port"

# Connection to QAI RESTful API for dumping project configuration (read-only)
qai_project_user = "FILLINUSERNAME"
qai_project_password = "****"
qai_project_path = "http://192.168.X.Y:port"
