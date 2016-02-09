#! /usr/bin/env python

"""
MISEQ_MONITOR.py
1) For runs flagged 'needsprocessing' that have not yet been processed, copy fastqs to local disk
2) Process the run by calling MISEQ_MONITOR.py
3) Upload results to the network drive
"""

import csv
from glob import glob
import itertools
import logging
import operator
import os
import sched
import shutil
import subprocess
import sys
import time
from xml.etree import ElementTree

from micall.core import miseq_logging
from micall.monitor import qai_helper
from micall.utils.sample_sheet_parser import sample_sheet_parser
import micall.settings as settings
from micall.monitor import update_qai
from micall.monitor.kive_download import download_results, kive_login
import hashlib


if sys.version_info[:2] != (2, 7):
    raise Exception("Python 2.7 not detected")

MAX_RUN_NAME_LENGTH = 60
logger = None


def init_logging(log_file):
    try:
        logger = miseq_logging.init_logging(log_file,
                                            file_log_level=logging.INFO,
                                            console_log_level=logging.INFO)
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
    if settings.production:
        with open(root + settings.ERROR_PROCESSING, 'w') as f:
            f.write(message)
    else:
        # in development mode - exit the monitor if a run fails
        sys.exit()
    return failure_message


def is_marked_as_disabled(run):
    return os.path.exists(run.replace(settings.NEEDS_PROCESSING, settings.ERROR_PROCESSING))


def is_quality_control_uploaded(run):
    return os.path.exists(run.replace(settings.NEEDS_PROCESSING, settings.QC_UPLOADED))


def mark_run_as_done(results_folder):
    """ Mark a run that has completed its processing.

    @param results_folder: path to the results folder
    """
    if settings.production:
        with open(os.path.join(results_folder, settings.DONE_PROCESSING), 'w'):
            pass  # Leave the file empty


def execute_command(command):
    logger.info(" ".join(command))
    stdout, stderr = subprocess.Popen(command,
                                      stdout=subprocess.PIPE,
                                      stderr=subprocess.PIPE).communicate()
    if stderr != "":
        logging.warn(stderr)
    return stdout


def post_files(files, destination):
    for f in files:
        execute_command(['rsync', '-a', f, '{}/{}'.format(destination, os.path.basename(f))])


def download_quality(run_info_path, destination, read_lengths, index_lengths):
    """ Download quality control data for the run.

    @return the QC run id as a string
    """
    runInfoTree = ElementTree.parse(run_info_path)
    runInfoRoot = runInfoTree.getroot()
    qcRunId = runInfoRoot[0].attrib['Id']
    direction_params = [(read_lengths[0], 1), (read_lengths[1], -1)]
    with qai_helper.Session() as session:
        session.login(settings.qai_path, settings.qai_user, settings.qai_password)
        metrics = session.get_json('/miseqqc_errormetrics?runid=' + qcRunId)
        if not metrics:
            raise StandardError(
                'No quality control metrics found for run ' + qcRunId)

    with open(destination, 'w') as f:
        writer = csv.DictWriter(f,
                                ['tile', 'cycle', 'errorrate'],
                                lineterminator=os.linesep)
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
                expected_cycle += sum(index_lengths)
    return qcRunId


#####################################################################################

logger = init_logging(settings.home + '/MISEQ_MONITOR_OUTPUT.log')

kive = kive_login(settings.kive_server_url,
                  settings.kive_user,
                  settings.kive_password)

# retrieve Pipeline object based on version
pipeline = kive.get_pipeline(settings.pipeline_version_kive_id)

# retrieve quality.csv compound data type
quality_cdt = kive.get_cdt(settings.quality_cdt_kive_id)


def check_kive(scheduler, active_runs):
    """
    :param scheduler: instance of Python sched
    :param active_runs: a List of Kive pipeline run objects
    :return:
    """
    sample_name, run = active_runs[0]
    if run.is_complete():
        logger.info('completed sample %s.', sample_name)
        active_runs.pop(0)
    if active_runs:
        # reschedule a check on runs at subsequent time
        scheduler.enter(5, 1, check_kive, (scheduler, active_runs))


def find_runs(processed_runs):
    # flag indicates that Illumina MiseqReporter has completed pre-processing, files available on NAS
    runs = glob(settings.rawdata_mount +
                'MiSeq/runs/*/{}'.format(settings.NEEDS_PROCESSING))
    runs_needing_processing = []
    for run in runs:
        root = os.path.dirname(run)
        result_path = os.path.join(root,
                                   'Results',
                                   'version_{}'.format(settings.pipeline_version))
        done_path = os.path.join(result_path, settings.DONE_PROCESSING)
        if is_marked_as_disabled(run):
            continue
        if not is_quality_control_uploaded(run):
            continue
        # if doneprocessing file already exists, then do not re-process
        if settings.production:
            if os.path.exists(done_path):
                continue
            if os.path.exists(result_path):
                # Not done, but results folder exists. Assume it's incomplete,
                # delete it, and rerun.
                shutil.rmtree(result_path)
            # running in development mode - do all runs even if already processed
            # note that results will not be uploaded
        elif run in processed_runs:
            continue
        runs_needing_processing.append(run)

    return runs_needing_processing


def format_run_name(root):
    return os.path.basename(os.path.normpath(root))


def trim_run_name(run_name):
    # run name is YYMMDD_MACHINEID_SEQ_RANDOMNUMBER.
    # Date and machine id should be unique together.
    return '_'.join(run_name.split('_')[:2])


def upload_dataset(filename, description):
    dataset_name = os.path.basename(filename)
    CHUNK_SIZE = 4096
    hash = hashlib.md5()
    with open(filename, 'rb') as f:
        for chunk in iter(lambda: f.read(CHUNK_SIZE), b""):
            hash.update(chunk)
        checksum = hash.hexdigest()
        datasets = kive.find_datasets(dataset_name=dataset_name, md5=checksum)
        needed_groups = set(settings.kive_groups_allowed)
        for dataset in datasets:
            missing_groups = needed_groups - set(dataset.groups_allowed)
            if not missing_groups:
                return datasets[0]

        f.seek(0)
        dataset = kive.add_dataset(name=dataset_name,
                                   description=description,
                                   handle=f,
                                   groups=settings.kive_groups_allowed)
    return dataset


def upload_data(root, run_folder):
    """ Upload FASTQ files and quality to Kive.

    :param root: the path to the root folder of the run that holds all the
        FASTQ files
    :param run_folder: the local folder that will hold working files.
    :returns: a list of sample names and FASTQ datasets, in tuples, a
        failure message string, and a dataset of quality data
        [(sample_name, read1_dataset, read2_dataset)], failure_message, quality
    """
    fastqs = []
    failure_message = None
    quality_input = None
    run_name = format_run_name(root)
    trimmed_run_name = trim_run_name(run_name)

    # transfer SampleSheet.csv
    remote_file = os.path.join(root, 'SampleSheet.csv')
    local_file = run_folder + '/SampleSheet.csv'
    execute_command(['rsync', '-a', remote_file, local_file])
    try:
        with open(local_file, 'rU') as sample_sheet:
            # parse run information from SampleSheet
            run_info = sample_sheet_parser(sample_sheet)
            read_lengths = run_info['Reads']
            entry = run_info['Data'].values()[0]
            indexes = [entry.get('index1', ''), entry.get('index2', '')]
            while 'X' in indexes:
                indexes.remove('X')
            index_lengths = map(len, indexes)
    except Exception:
        failure_message = mark_run_as_disabled(root,
                                               "Parsing sample sheet failed",
                                               exc_info=True)
        return fastqs, failure_message, quality_input

    gz_files = glob(root + 'Data/Intensities/BaseCalls/*R?_001.fastq.gz')
    if not gz_files:
        failure_message = mark_run_as_disabled(root, "No data files found")
        return fastqs, failure_message, quality_input

    # use Kive-API to transfer fastq.gz files
    R1_files = sorted(filter(lambda x: '_R1_' in os.path.basename(x), gz_files))
    for R1_file in R1_files:
        R2_file = R1_file.replace('_R1_', '_R2_')
        if R2_file not in gz_files:
            logger.warn("Detected an unpaired R1 FASTQ file: {}".format(R1_file))
            continue

        filename = os.path.basename(R1_file)
        sample, snum = filename.split('_')[:2]
        # Report number of reads failing to demultiplex to the log
        if filename.startswith('Undetermined'):
            # Second file has the same number of lines - don't bother
            if filename.endswith('_L001_R2_001.fastq.gz'):
                continue

            # Do word count directly on stream redirected from gunzip
            p1 = subprocess.Popen(['gunzip', '-c', R1_file], stdout=subprocess.PIPE)
            p2 = subprocess.Popen(['wc', '-l'], stdin=p1.stdout, stdout=subprocess.PIPE)
            output = p2.communicate()[0]
            failed_demultiplexing = output.strip(' \n')
            logger.info("%s reads failed to demultiplex in %s (removing file)", failed_demultiplexing, filename)
            continue

        logger.info(filename)
        description = 'FASTQ for sample %s (%s) from MiSeq run %s' % (
            sample,
            snum,
            trimmed_run_name)
        R1_obj = upload_dataset(R1_file, 'R1 ' + description)
        R2_obj = upload_dataset(R2_file, 'R2 ' + description)
        fastqs.append(((sample + '_' + snum), R1_obj, R2_obj))

    # generate quality.csv
    try:
        quality_csv = os.path.join(settings.home,
                                   run_name,
                                   '{}_quality.csv'.format(trimmed_run_name))
        download_quality(run_info_path=os.path.join(root, 'RunInfo.xml'),
                         destination=quality_csv,
                         read_lengths=read_lengths,
                         index_lengths=index_lengths)
    except StandardError:
        failure_message = mark_run_as_disabled(root,
                                               "Quality could not be downloaded.",
                                               exc_info=True)
        return fastqs, failure_message, quality_input

    # transfer quality.csv with Kive API
    quality_input = upload_dataset(
        quality_csv,
        'phiX174 quality per tile and cycle for run ' + run_name)
    return fastqs, failure_message, quality_input


def launch_runs(fastqs, quality_input, run_name):
    """ Launch runs on Kive, then wait for them to finish.

    :param fastqs: a list of sample names and FASTQ datasets, in tuples
        [(sample_name, read1_dataset, read2_dataset)]
    :param quality_input: a dataset with quality data
    :param run_name: the name of the run for labeling the Kive runs
    :return: a list of sample names and RunStatus objects, in tuples
        [(sample_name, run_status)]
    """
    trimmed_run_name = trim_run_name(run_name)
    # Standard out/error concatenates to the log
    logger.info("Launching pipeline for %s%s", settings.home, run_name)
    kive_runs = []
    # push all samples into the queue
    for sample_name, fastq1, fastq2 in fastqs:
        name = '{} - {} ({})'.format(pipeline.family,
                                     sample_name,
                                     trimmed_run_name)
        name = name[:MAX_RUN_NAME_LENGTH]
        # note order of inputs is critical
        status = kive.run_pipeline(pipeline=pipeline,
                                   inputs=[quality_input, fastq1, fastq2],
                                   name=name,
                                   groups=settings.kive_groups_allowed)
        kive_runs.append((sample_name, status))

    # initialize progress monitoring
    sc = sched.scheduler(time.time, time.sleep)
    active_runs = kive_runs[:]
    sc.enter(5, 1, check_kive, (sc, active_runs))
    sc.run()  # exits when all runs are complete
    logger.info("===== {} successfully processed! =====".format(run_name))
    return kive_runs


def main():
    global logger
    processed_runs = set()
    failure_message = None
    # main loop
    # Process runs flagged for processing not already processed by this version of the pipeline
    while True:
        if failure_message is not None:
            # We're still logging to a run folder that failed.
            logging.shutdown()
            logger = None

        if logger is None:
            logger = init_logging(settings.home + '/MISEQ_MONITOR_OUTPUT.log')

        if failure_message is not None:
            # Log it here again, now logging to home folder's log file.
            logger.error(failure_message)
            failure_message = None

        runs_needing_processing = find_runs(processed_runs)

        if not runs_needing_processing:
            logger.info('No runs need processing')
            if not settings.production:
                # development mode - exit main loop
                break
            time.sleep(settings.delay)
            continue

        # Process most recently generated run and work backwards
        runs_needing_processing.sort(reverse=True)
        curr_run = runs_needing_processing[0]
        root = curr_run.replace(settings.NEEDS_PROCESSING, '')
        run_name = format_run_name(root)
        run_folder = os.path.join(settings.home, run_name)

        # Make folder on the cluster for intermediary files and outputs
        if not os.path.exists(run_folder):
            os.mkdir(run_folder)

        # Determine if intermediary files should persist past end of run
        all_runs = map(lambda x: x.split('/')[-1], glob('%s*%s*' % (settings.home, settings.instrument_number)))
        all_runs.sort()
        runs_to_keep = all_runs[-settings.nruns_to_store:]  # X most recent runs
        do_cleanup = (run_name not in runs_to_keep)  # true/false
        if do_cleanup:
            logger.info('Set clean mode for run %s', run_name)

        # Record standard input / output of monitor
        logger.info('Starting run %s', root)
        logging.shutdown()
        log_file = run_folder + '/MISEQ_MONITOR_OUTPUT.log'
        logger = init_logging(log_file)
        logger.info('===== Processing {} with pipeline version {} ====='.format(root, settings.pipeline_version))

        fastqs, failure_message, quality_input = upload_data(root, run_folder)
        if not (fastqs and quality_input):
            continue

        try:
            kive_runs = launch_runs(fastqs, quality_input, run_name)
        except Exception as e:
            failure_message = mark_run_as_disabled(
                root,
                "Failed to launch runs: '{}'".format(e),
                exc_info=True)
            continue

        try:
            if not settings.production:
                results_folder = os.path.join(run_folder, 'results')
            else:
                results_parent = os.path.join(root, 'Results')
                if not os.path.exists(results_parent):
                    os.mkdir(results_parent)
                results_folder = os.path.join(results_parent,
                                              'version_' + settings.pipeline_version)
            if not os.path.exists(results_folder):
                os.mkdir(results_folder)
            download_results(kive_runs, results_folder, run_folder)
            if settings.production:
                update_qai.process_folder(results_folder, logger)
            mark_run_as_done(results_folder)
            logger.info("===== %s file transfer completed =====", run_name)
        except:
            failure_message = mark_run_as_disabled(
                root,
                "Failed to post pipeline results",
                exc_info=True)
        if not settings.production:
            processed_runs.add(curr_run)

        # Reset logger so it will switch back to the root directory log file.
        logging.shutdown()
        logger = None
main()
