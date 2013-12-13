import logging, miseq_logging, miseq_modules, os, sys

file = sys.argv[1]
qcut = int(sys.argv[2])
read_mapping_cutoff = int(sys.argv[3])
mode = sys.argv[4]
max_prop_N = float(sys.argv[5])

logger = miseq_logging.init_logging_console_only(logging.INFO)
logger.info("pid {}: miseq_modules.sam2csf_with_base_censoring({}, {})".format(os.getpid(),file, qcut))
miseq_modules.sam2csf_with_base_censoring(file, qcut, read_mapping_cutoff, mode, max_prop_N)
logging.shutdown()
