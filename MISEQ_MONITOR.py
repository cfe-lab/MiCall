"""
MISEQ_MONITOR.py
1) For runs flagged 'needsprocessing' that have not yet been processed, copy fastqs to local disk
2) Process the run by calling MISEQ_PIPELINE.py
3) Upload results to the network drive
"""

import csv
from glob import glob
import logging
import os
import shutil
import subprocess
import sys
import tarfile
import time
from xml.etree import ElementTree

import miseq_logging
import qai_helper
from sample_sheet_parser import sample_sheet_parser
from settings import delay, DONE_PROCESSING, ERROR_PROCESSING, home,\
    NEEDS_PROCESSING, pipeline_version, production, rawdata_mount, base_path,\
    qai_path, qai_user, qai_password, instrument_number, nruns_to_store
import update_qai
import itertools
import operator


if sys.version_info[:2] != (2, 7):
    raise Exception("Python 2.7 not detected")

def init_logging(log_file):
    try:
        logger = miseq_logging.init_logging(log_file, file_log_level=logging.DEBUG, console_log_level=logging.INFO)
    except Exception as e:
        raise Exception("Couldn't setup logging (init_logging() threw exception '{}') - HALTING NOW!".format(str(e)))
    return logger

def mark_run_as_disabled(root, message, exc_info=None):
    """ Mark a run that failed, so it won't be processed again.
    
    @param root: path to the run folder that had an error
    @param message: a description of the error
    @param exc_info: details about the error's exception in the standard tuple,
        True to look up the current exception, or None if there is no exception
        to report
    """
    failure_message = message + " - skipping run " + root
    logger.error(failure_message, exc_info=exc_info)
    if production:
        with open(root + ERROR_PROCESSING, 'w') as f:
            f.write(message)
    
    return failure_message

def is_marked_as_disabled(run):
    return os.path.exists(run.replace(NEEDS_PROCESSING, ERROR_PROCESSING))

def mark_run_as_done(results_folder):
    """ Mark a run that has completed its processing.
    
    @param results_folder: path to the results folder
    """
    if production:
        with open(os.path.join(results_folder, DONE_PROCESSING), 'w'):
            pass # Leave the file empty

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
failure_message = None

# Process runs flagged for processing not already processed by this version of the pipeline


def download_quality(run_info_path, destination, read_lengths):
    """ Download quality control data for the run.
    
    @return the QC run id as a string
    """
    runInfoTree = ElementTree.parse(run_info_path)
    runInfoRoot = runInfoTree.getroot()
    qcRunId = runInfoRoot[0].attrib['Id']
    direction_params = [(read_lengths[0], 1), (read_lengths[1], -1)]
    with qai_helper.Session() as session:
        session.login(qai_path, qai_user, qai_password)
        metrics = session.get_json('/miseqqc_errormetrics?runid=' + qcRunId)
    with open(destination, 'w') as f:
        writer = csv.DictWriter(f, ['tile', 'cycle', 'errorrate'])
        writer.writeheader()
        for tile, tile_metrics in itertools.groupby(metrics,
                                                    operator.itemgetter('tile')):
            expected_cycle = 0
            metric = next(tile_metrics)
            for cycles, sign in direction_params:
                report_cycle = sign
                for _ in range(1, cycles+1):
                    expected_cycle += 1
                    cycle = int(metric['cycle'])
                    if cycle == expected_cycle:
                        errorrate = metric['errorrate']
                        try:
                            metric = next(tile_metrics)
                        except StopIteration:
                            metric = dict(cycle=-1)
                    else:
                        errorrate = ''
                    writer.writerow(dict(tile=tile,
                                         cycle=report_cycle,
                                         errorrate=errorrate))
                    report_cycle += sign
                expected_cycle += 16
    return qcRunId

while True:
    if failure_message is not None:
        # We're still logging to a run folder that failed.
        logging.shutdown()
        logger = None
        
    if logger is None:
        logger = init_logging(home + '/MISEQ_MONITOR_OUTPUT.log')
    
    if failure_message is not None:
        # Log it here again, now logging to home folder's log file.
        logger.error(failure_message)
        failure_message = None
        
    # flag indicates that Illumina MiseqReporter has completed pre-processing, files available on NAS
    runs = glob(rawdata_mount + 'MiSeq/runs/*/{}'.format(NEEDS_PROCESSING))
    #runs = glob(rawdata_mount + 'MiSeq/runs/131119_M01841_0041_000000000-A5EPY/{}'.format(NEEDS_PROCESSING))

    runs_needing_processing = []
    for run in runs:
        root = os.path.dirname(run)
        result_path = os.path.join(root,
                                   'Results',
                                   'version_{}'.format(pipeline_version))
        done_path = os.path.join(result_path, DONE_PROCESSING)

        if is_marked_as_disabled(run):
            continue

        # if doneprocessing file already exists, then do not re-process
        if production:
            if os.path.exists(done_path):
                continue
            if os.path.exists(result_path):
                # Not done, but results folder exists. Assume it's incomplete,
                # delete it, and rerun.
                shutil.rmtree(result_path)
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

    all_runs = map(lambda x: x.split('/')[-1], glob('%s/*/*%s*' % (home, instrument_number)))
    all_runs.sort()
    runs_to_keep = all_runs[-nruns_to_store:]  # X most recent runs
    do_cleanup = (run_name in runs_to_keep)


    # Record standard input / output of monitor
    logger.info('Starting run %s', root)
    logging.shutdown()
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
            read_lengths = run_info['Reads']
    except Exception as e:
        failure_message = mark_run_as_disabled(
            root,
            "Parsing sample sheet failed",
            exc_info=True)
        continue

    if mode not in ['Nextera', 'Amplicon']:
        failure_message = mark_run_as_disabled(
            root,
            "{} not a valid mode".format(mode))
        continue

    # Copy fastq.gz files to the cluster and unzip them
    gz_files = glob(root + 'Data/Intensities/BaseCalls/*R?_001.fastq.gz')
    if not gz_files:
        failure_message = mark_run_as_disabled(root, "No data files found")
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

    try:
        download_quality(os.path.join(root, 'RunInfo.xml'),
                         os.path.join(home, run_name, 'quality.csv'),
                         read_lengths)
    except StandardError as e:
        failure_message = mark_run_as_disabled(root,
                                               "Quality could not be downloaded.",
                                               exc_info=True)
        continue

    # Standard out/error concatenates to the log
    logger.info("Launching pipeline for %s%s", home, run_name)
    try:
        subprocess.check_call([os.path.join(base_path, 'run_processor.py'),
                               home+run_name,
                               '-clean' if do_cleanup else ''])
        logger.info("===== {} successfully processed! =====".format(run_name))
    except Exception as e:
        failure_message = mark_run_as_disabled(
            root,
            "MISEQ_PIPELINE.py failed: '{}'".format(e))
        continue

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
                         ('bad_cycles.csv', None),
                         ('*.log', 'logs'),
                         ('*.unmapped?.fastq', 'unmapped')]
            
            if not os.path.isdir(result_path):
                os.mkdir(result_path)
            for pattern, destination in file_sets:
                destination_path = (
                    os.path.join(result_path_final, destination)
                    if destination
                    else result_path_final)
                full_pattern = os.path.join(home, run_name, pattern)
                if not os.path.isdir(destination_path):
                    os.mkdir(destination_path)
                post_files(glob(full_pattern), destination_path)
                    
            untar_path = os.path.join(result_path_final, 'untar')
            coverage_source_path = os.path.join(untar_path, 'coverage_maps')
            coverage_dest_path = os.path.join(result_path_final, 'coverage_maps')
            os.mkdir(untar_path)
            os.mkdir(coverage_dest_path)
            for tar_path in glob(home + run_name + '/*.coverage_maps.tar'):
                basename = os.path.basename(tar_path)
                map_name = os.path.splitext(basename)[0]
                sample_name = os.path.splitext(map_name)[0]
                with tarfile.open(tar_path) as tar:
                    tar.extractall(untar_path)
                for image_filename in os.listdir(coverage_source_path):
                    source = os.path.join(coverage_source_path, image_filename)
                    destination = os.path.join(coverage_dest_path,
                                               sample_name+'.'+image_filename)
                    os.rename(source, destination)
            os.rmdir(coverage_source_path)
            os.rmdir(untar_path)
            
            update_qai.process_folder(result_path_final, logger)
            
            mark_run_as_done(result_path_final)
            
            logger.info("===== %s file transfer completed =====", run_name)
        except:
            failure_message = mark_run_as_disabled(
                root,
                "Failed to post pipeline results",
                exc_info=True)

    logging.shutdown()
    logger = None
