import logging, miseq_logging, miseq_modules, os, sys

logger = miseq_logging.init_logging_console_only(logging.DEBUG)

file = sys.argv[1]
mode = sys.argv[2]
ref_path = sys.argv[3]

logger.info("pid {}: miseq_modules.csf2nuc({},{},{})".format(os.getpid(), file, mode, ref_path))
miseq_modules.csf2nuc(file, mode, ref_path)
logging.shutdown()
