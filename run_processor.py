#!/usr/bin/env python

import argparse
from glob import glob
import logging
import os
import subprocess
import sys
import traceback

import miseq_logging
import settings
import shutil

def parseOptions():
    parser = argparse.ArgumentParser(
        description='Process all the samples in a single run folder.')
    
    parser.add_argument('run_folder',
                        help='Path to sample fastq files and SampleSheet.csv')
    parser.add_argument('mode',
                        help='Amplicon or Nextera, default from sample sheet',
                        nargs=argparse.OPTIONAL)
    
    return parser.parse_args()



def check_mpi_version(prefix):
    try:
        return subprocess.check_output(prefix + 'mpirun -V',
                                       shell=True,
                                       stderr=subprocess.STDOUT)
    except:
        etype, value, _tb = sys.exc_info()
        return traceback.format_exception_only(etype, value)

def main():
    args = parseOptions()
    log_file = "{}/run.log".format(args.run_folder)
    logger = miseq_logging.init_logging(log_file,
                                        file_log_level=logging.DEBUG,
                                        console_log_level=logging.INFO)
    
    logger.info('Start processing run %s', args.run_folder)
    logger.info('Removing old working files')
    old_files = glob(args.run_folder+'/*')
    for f in old_files:
        if not f.endswith('.fastq') and not f.endswith('SampleSheet.csv'):
            if os.path.isdir(f):
                shutil.rmtree(f)
            else:
                os.remove(f)
        
    # Check for Open MPI
    prefix = ''
    expected_version = 'Open MPI'
    version = check_mpi_version(prefix)
    if not expected_version in version:
        prefix = 'module load openmpi/gnu && '
        version = check_mpi_version(prefix)
        if not expected_version in version:
            sys.exit("Couldn't find Open MPI:\n{}".format(version))
    
    mapping_args = ['mpirun', 
                    '-np', 
                    str(settings.mapping_processes), 
                    '--hostfile', 
                    'hostfile', 
                    './sample_pipeline.py',
                    args.run_folder]
    if args.mode is not None:
        mapping_args.append(args.mode)
    mapping_args.append('--phase')
    mapping_args.append('mapping')
    mapping_command = prefix + ' '.join(mapping_args)
    
    counting_args = mapping_args[:]
    counting_args[2] = str(settings.counting_processes)
    counting_args[-1] = 'counting'
    counting_command = prefix + ' '.join(counting_args)
    
    summarizing_args = mapping_args[:]
    summarizing_args[2] = '1' # Only one process does the summary
    summarizing_args[-1] = 'summarizing'
    summarizing_command = prefix + ' '.join(summarizing_args)
    
    subprocess.check_call(mapping_command, shell=True)

    subprocess.check_call(counting_command, shell=True)

    subprocess.check_call(summarizing_command, shell=True)
    
    logger.info('Finished processing run %s', args.run_folder)

if __name__ == '__main__':
    main()
