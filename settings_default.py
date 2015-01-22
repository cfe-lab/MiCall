"""
To make pipeline portable, allow user to specify local paths and thread counts.
"""

pipeline_version = '6.6'

production = False  # set this to True to push results to NAS
filter_cross_contaminants = False

## Modify this path to your pipeline directory
base_path = '/usr/local/share/miseq/'
base_path += 'production/' if production else 'development/'

# Local path for writing data
home = '/data/miseq/'


are_temp_folders_deleted = True # Should FIFO worker clean up working folders?

# Scheduling processes: these should be a multiple of the total number of slots
# in your hostfile.
mapping_processes = 38
counting_processes = 152


## MISEQ_MONITOR settings

rawdata_mount = '/media/RAW_DATA/'  # NAS

delay = 3600  # Delay (seconds) for polling NAS for unprocessed runs

NEEDS_PROCESSING = 'needsprocessing'  # File flags
ERROR_PROCESSING = 'errorprocessing'
DONE_PROCESSING = 'doneprocessing'


## Mapping parameters

projects_json = base_path + 'projects.json'

bowtie_threads = 4              # Bowtie performance roughly scales with number of threads
min_mapping_efficiency = 0.95   # Fraction of fastq reads mapped needed
max_remaps = 3                  # Number of remapping attempts if mapping efficiency unsatisfied
consensus_q_cutoff = 20         # Min Q for base to contribute to conseq (pileup2conseq)


## sam2aln parameters
sam2aln_q_cutoffs = [15]  # Q-cutoff for base censoring
max_prop_N = 0.5                 # Drop reads with more censored bases than this proportion
read_mapping_cutoff = 10         # Minimum bowtie read mapping quality
insert_qcutoff = 20              # minimum Q score for an insertion polymorphism to pass

## g2p parameters (Amplicon only)
alignment_lib = 'alignment.so'
g2p_alignment_cutoff = 50               # Minimum alignment score during g2p scoring
g2p_fpr_cutoffs = [3.0, 3.5, 4.0, 5.0]  # FPR cutoff to determine R5/X4 tropism
v3_mincounts = [0, 50, 100, 1000]       # Min number of reads to contribute to %X4 calculation

## aln2counts parameters
conseq_mixture_cutoffs = [0.01, 0.02, 0.05, 0.1, 0.2, 0.25]

# Connection to QAI RESTful API for uploading results
qai_user = "FILLINUSERNAME"
qai_password = "****"
qai_path = "http://192.168.X.Y:port"

# Connection to QAI RESTful API for dumping project configuration (read-only)
qai_project_user = "FILLINUSERNAME"
qai_project_password = "****"
qai_project_path = "http://192.168.X.Y:port"
