import logging, miseq_logging, miseq_modules, os, sys

csf_file = sys.argv[1]
g2p_alignment_cutoff = float(sys.argv[2])

logger = miseq_logging.init_logging_console_only(logging.INFO)
logger.info("pid {}: miseq_modules.g2p_scoring({}, {})".format(os.getpid(), csf_file, g2p_alignment_cutoff))
miseq_modules.g2p_scoring(csf_file, g2p_alignment_cutoff)
logging.shutdown()
