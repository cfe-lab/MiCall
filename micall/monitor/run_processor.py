#!/usr/bin/env python

import argparse
from glob import glob
import logging
import os
import shutil
import subprocess
import sys
import traceback

from micall.core import miseq_logging
from micall import settings

def parseOptions():
    parser = argparse.ArgumentParser(
        description='Process all the samples in a single run folder.')
    
    parser.add_argument('run_folder',
                        help='Path to sample fastq files and SampleSheet.csv')
    parser.add_argument('mode',
                        help='Amplicon or Nextera, default from sample sheet',
                        nargs=argparse.OPTIONAL)
    parser.add_argument('--clean',
                        help='Remove intermediate files after run is complete.',
                        action='store_true')
    
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
    if args.clean:
        logger.info('Clean mode ON')
    else:
        logger.info('Clean mode OFF')

    try:
        logger.info('Removing old working files')
        excluded_files = ('.fastq',
                          'SampleSheet.csv',
                          '.launch',
                          'MISEQ_MONITOR_OUTPUT.log',
                          'run.log',
                          'quality.csv')
        old_files = glob(args.run_folder+'/*')
        for f in old_files:
            is_excluded = False
            for ending in excluded_files:
                if f.endswith(ending):
                    is_excluded = not (f.endswith('unmapped1.fastq') or
                                       f.endswith('unmapped2.fastq') or
                                       f.endswith('censored1.fastq') or
                                       f.endswith('censored2.fastq'))
                    break
            if not is_excluded:
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
        monitor_path = os.path.abspath(os.path.dirname(__file__))
        
        base_args =    ['mpirun', 
                        '-np', 
                        '1', 
                        '--hostfile', 
                        os.path.join(monitor_path, 'hostfile'), 
                        os.path.join(monitor_path, 'sample_pipeline.py'),
                        args.run_folder]

        if args.mode is not None:
            base_args.append(args.mode)
        base_args.append('--phase')
        
        filter_args = base_args[:]
        filter_args.append('filter')
        filter_command = prefix + ' '.join(filter_args)
        
        mapping_args = base_args[:]
        mapping_args[2] = str(settings.mapping_processes)
        mapping_args.append('mapping')
        mapping_command = prefix + ' '.join(mapping_args)
        
        counting_args = base_args[:]
        counting_args[2] = str(settings.counting_processes)
        counting_args.append('counting')
        counting_command = prefix + ' '.join(counting_args)
        
        summarizing_args = base_args[:]
        summarizing_args.append('summarizing')
        summarizing_command = prefix + ' '.join(summarizing_args)
        
        subprocess.check_call(filter_command, shell=True)
        
        subprocess.check_call(mapping_command, shell=True)
    
        subprocess.check_call(counting_command, shell=True)
    
        subprocess.check_call(summarizing_command, shell=True)

        if args.clean:
            # remove intermediate files
            logger.info('Removing large working files, clean mode')
            files_to_remove = glob(args.run_folder+'/*.prelim.csv')
            files_to_remove += glob(args.run_folder+'/*.remap.csv')
            files_to_remove += glob(args.run_folder+'/*.aligned.csv')
            files_to_remove += glob(args.run_folder+'/*.censored?.fastq')
            for f in files_to_remove:
                os.remove(f)

        logger.info('Finished processing run %s', args.run_folder)
    except:
        logger.error('Failed to process run %s', args.run_folder, exc_info=True)
        exit(1)

if __name__ == '__main__':
    main()
