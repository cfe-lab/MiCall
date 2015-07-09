"""
To make pipeline portable, allow user to specify local paths and thread counts.
"""
import os

### Core settings ###
pipeline_version = '6.8'        # Change for each release
base_path = os.path.dirname(os.path.realpath(__file__))
projects_json = os.path.join(base_path, 'projects.json')
are_temp_folders_deleted = True # Should FIFO worker clean up working folders?

## Mapping parameters
bowtie_version = '2.2.1'        # version of bowtie2, used for version control
bowtie_path = 'bowtie2'         # path to executable, so you can install more than one version
bowtie_build_path = 'bowtie2-build'
bowtie_threads = 4              # Bowtie performance roughly scales with number of threads
samtools_version = '0.1.18'
samtools_path = 'samtools'
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


### Monitor settings ###
instrument_number = 'M01841'  # for Illumina MiSeq, second item in run folder name
production = False  # set this to True to push results to NAS

# Local path for writing data
home = '/data/miseq/'
nruns_to_store = 20  # protect X most recent runs from cleaning up intermediate files

# Scheduling processes: these should be a multiple of the total number of slots
# in your hostfile.
mapping_processes = 38
counting_processes = 152

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
