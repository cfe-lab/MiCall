import miseq_modules,sys

mapping_ref_path = sys.argv[1]
fastq = sys.argv[2]
consensus_q_cutoff = int(sys.argv[3])
mode = sys.argv[4]
is_t_primer = sys.argv[5]
REMAP_THRESHOLD = float(sys.argv[6])
MAX_REMAPS = int(sys.argv[7])
num_threads = int(sys.argv[8])

miseq_modules.mapping(mapping_ref_path, fastq, consensus_q_cutoff, mode, is_t_primer, REMAP_THRESHOLD, MAX_REMAPS, num_threads)
