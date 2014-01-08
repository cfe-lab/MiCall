import logging, miseq_logging, miseq_modules, os, sys

# python STEP_4_CSF2COUNTS.py /data/miseq/A5F9R/blah.0.csf Nextera 0.2 /usr/local/share/miseq/refs/csf2counts_amino_refseqs.csv

logger = miseq_logging.init_logging_console_only(logging.DEBUG)

file = sys.argv[1]
mode = sys.argv[2]
consensus_mixture_cutoffs = map(float,(sys.argv[3].split(",")))
final_alignment_ref_path = sys.argv[4]

logger.info("pid {}: miseq_modules.csf2counts({},{},{},{})".format(os.getpid(), file, mode, consensus_mixture_cutoffs, final_alignment_ref_path))
miseq_modules.csf2counts(file,mode,consensus_mixture_cutoffs,final_alignment_ref_path)
logging.shutdown()
