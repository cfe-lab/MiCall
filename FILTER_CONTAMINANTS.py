import logging
import miseq_logging
import os
import miseq_modules
import sys


logger = miseq_logging.init_logging_console_only(logging.DEBUG)

# arguments
root = sys.argv[1] # folder with data
qcut = int(sys.argv[2]) # base quality score cutoff
num_threads = int(sys.argv[3]) # for bowtie2

logger.info("pid {}: miseq_modules.filter_cross_contaminants({},{},{})".format(os.getpid(), root, qcut, num_threads))
miseq_modules.filter_cross_contaminants(root=root, qcutoff=qcut, num_threads=num_threads)
logging.shutdown()
