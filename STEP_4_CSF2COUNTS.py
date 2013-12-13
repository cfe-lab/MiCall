import logging, miseq_logging, miseq_modules, os, sys

file = sys.argv[1]
mode = sys.argv[2]
consensus_mixture_cutoffs = map(float,(sys.argv[3].split(",")))
final_alignment_ref_path = sys.argv[4]

logger = miseq_logging.init_logging_console_only(logging.INFO)
logger.info("pid {}: miseq_modules.csf2counts({},{},{},{})".format(os.getpid(), file, mode, consensus_mixture_cutoffs, final_alignment_ref_path))
miseq_modules.csf2counts(file,mode,consensus_mixture_cutoffs,final_alignment_ref_path)
logging.shutdown()
