"""
MISEQ_MONITOR.py
1) For runs flagged 'needsprocessing' that have not yet been processed, copy fastqs to local disk
2) Process the run by calling MISEQ_PIPELINE.py
3) Upload results to the network drive
"""

from glob import glob
import logging
import os
import subprocess
import sys
import tarfile
import time

import miseq_logging
from sample_sheet_parser import sample_sheet_parser
from settings import delay, ERROR_PROCESSING, home, NEEDS_PROCESSING,\
    pipeline_version, production, rawdata_mount, base_path
import update_oracle


if sys.version_info[:2] != (2, 7):
    raise Exception("Python 2.7 not detected")

def init_logging(log_file):
    try:
        logger = miseq_logging.init_logging(log_file, file_log_level=logging.DEBUG, console_log_level=logging.INFO)
    except Exception as e:
        raise Exception("Couldn't setup logging (init_logging() threw exception '{}') - HALTING NOW!".format(str(e)))
    return logger

def mark_run_as_disabled(root, message):
    """ Mark a run that failed, so it won't be processed again. """
    logger.error(message + " - skipping run " + root)
    if production:
        with open(root + ERROR_PROCESSING, 'w') as f:
            f.write(message)

def is_marked_as_disabled(run):
    return os.path.exists(run.replace(NEEDS_PROCESSING, ERROR_PROCESSING))

def execute_command(command):
    logger.info(" ".join(command))
    stdout, stderr = subprocess.Popen(command, stdout = subprocess.PIPE, stderr = subprocess.PIPE).communicate()
    if stderr != "": logging.warn(stderr)
    return stdout

def post_files(files, destination):
    for f in files:
        execute_command(['rsync', '-a', f, '{}/{}'.format(destination, os.path.basename(f))])

processed_runs = set()
logger = None

# Process runs flagged for processing not already processed by this version of the pipeline
while True:
    if logger is None:
        logger = init_logging(home + '/MISEQ_MONITOR_OUTPUT.log')
    # flag indicates that Illumina MiseqReporter has completed pre-processing, files available on NAS
    runs = glob(rawdata_mount + 'MiSeq/runs/*/{}'.format(NEEDS_PROCESSING))
    #runs = glob(rawdata_mount + 'MiSeq/runs/131119_M01841_0041_000000000-A5EPY/{}'.format(NEEDS_PROCESSING))

    runs_needing_processing = []
    for run in runs:
        result_path = '{}/version_{}'.format(run.replace(NEEDS_PROCESSING, 'Results'), pipeline_version)

        if is_marked_as_disabled(run):
            continue

        # if version-matched Results folder already exists, then do not re-process
        if production:
            if os.path.exists(result_path):
                continue
        else:
            if run in processed_runs:
                continue

        runs_needing_processing.append(run)

    if not runs_needing_processing:
        logger.info('No runs need processing')
        if not production:
            break
        time.sleep(delay)
        continue

    # Process most recently generated run and work backwards
    runs_needing_processing.sort(reverse=True)
    curr_run = runs_needing_processing[0]
    root = curr_run.replace(NEEDS_PROCESSING, '')
    run_name = root.split('/')[-2]

    # Make folder on the cluster for intermediary files
    if not os.path.exists(home+run_name):
        os.mkdir(home+run_name)

    # Record standard input / output of monitor
    log_file = home + run_name + '/MISEQ_MONITOR_OUTPUT.log'
    logger = init_logging(log_file)
    logger.info('===== Processing {} with pipeline version {} ====='.format(root, pipeline_version))

    # transfer SampleSheet.csv
    remote_file = curr_run.replace(NEEDS_PROCESSING, 'SampleSheet.csv')
    local_file = home + run_name + '/SampleSheet.csv'
    execute_command(['rsync', '-a', remote_file, local_file])

    # Delete raw tiffs from this run
    # logger.info("Deleting tiffs from {}Images/Focus/L001/*".format(root))
    #for tiff in glob("{}Images/Focus/L001/*/*.tif".format(root)):
    #    log.info("os.remove({})".format(tiff))
    #    os.remove(tiff)

    try:
        with open(local_file, 'rU') as sample_sheet:
            # parse run information from SampleSheet
            run_info = sample_sheet_parser(sample_sheet)
            mode = run_info['Description']
    except Exception as e:
        mark_run_as_disabled(
            root,
            "Parsing sample sheet failed: '{}'".format(e))
        continue

    if mode not in ['Nextera', 'Amplicon']:
        mark_run_as_disabled(root, "{} not a valid mode".format(mode))
        continue

    # Copy fastq.gz files to the cluster and unzip them
    gz_files = glob(root + 'Data/Intensities/BaseCalls/*R?_001.fastq.gz')
    if not gz_files:
        mark_run_as_disabled(root, "No data files found")
        continue
    
    for gz_file in gz_files:
        filename = os.path.basename(gz_file)

        # Report number of reads failing to demultiplex to the log
        if filename.startswith('Undetermined'):

            # Second file has the same number of lines - don't bother
            if filename.endswith('_L001_R2_001.fastq.gz'):
                continue

            # Do word count directly on stream redirected from gunzip
            p1 = subprocess.Popen(['gunzip', '-c', gz_file], stdout=subprocess.PIPE)
            p2 = subprocess.Popen(['wc', '-l'], stdin=p1.stdout, stdout=subprocess.PIPE)
            output = p2.communicate()[0]
            failed_demultiplexing = output.strip(' \n')
            logger.info("{} reads failed to demultiplex in {} (removing file)".format(failed_demultiplexing, filename))
            continue


        local_file = home + run_name + '/' + filename

        if os.path.exists(local_file.replace('.gz', '')):
            # If a local copy of the unzipped fastq exists, skip to next file
            continue

        # otherwise, transfer fastq.gz file and unzip
        execute_command(['rsync', '-a', gz_file, local_file])
        execute_command(['gunzip', '-f', local_file])


    # Standard out/error concatenates to the log
    logger.info("Launching pipeline for %s%s", home, run_name)
    try:
        subprocess.check_call([os.path.join(base_path, 'run_processor.py'),
                               home+run_name])
        logger.info("===== {} successfully completed! =====".format(run_name))
    except Exception as e:
        mark_run_as_disabled(root, "MISEQ_PIPELINE.py failed: '{}'".format(e))

    if not production:
        processed_runs.add(curr_run)
    else:
        try:
            # Determine output paths
            result_path = curr_run.replace(NEEDS_PROCESSING, 'Results')
            result_path_final = '{}/version_{}'.format(result_path, pipeline_version)
    
            # Post files to appropriate sub-folders
            logger.info("Posting results to {}".format(result_path_final))
            file_sets = [# (pattern, destination)
                         ('results/*', None),
                         ('*.log', 'logs'),
                         ('*.unmapped?.fastq', 'unmapped')]
            
            for pattern, destination in file_sets:
                destination_path = (
                    os.path.join(result_path_final, destination)
                    if destination
                    else result_path_final)
                full_pattern = os.path.join(home, run_name, pattern)
                if not os.path.isdir(destination_path):
                    os.mkdir(destination_path)
                post_files(glob(full_pattern), destination_path)
                    
            for tar_path in glob(home + run_name + '/*.coverage_maps.tar'):
                with tarfile.open(tar_path) as tar:
                    tar.extractall(result_path_final)
            
            update_oracle.process_folder(result_path_final, logger)
            
            # Close the log and copy it to rawdata
            logger.info("===== %s file transfer completed =====", run_name)
        except:
            logger.error('Failed to post pipeline results for %r.', 
                         run_name, 
                         exc_info=True)
            mark_run_as_disabled(root, "Failed to post pipeline results")

    logging.shutdown()
    logger = None
