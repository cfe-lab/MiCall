#! /usr/bin/env python

"""
MISEQ_MONITOR.py
1) For runs flagged 'needsprocessing' that have not yet been processed, copy fastqs to local disk
2) Process the run by calling MISEQ_MONITOR.py
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

from micall.core import miseq_logging
from micall.monitor import qai_helper
from micall.utils.sample_sheet_parser import sample_sheet_parser
import micall.settings as settings
from micall.monitor import update_qai
from micall.utils.collate import collate_labeled_files
import itertools
import operator


from kiveapi import KiveAPI
import sched
import time



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
            pass # Leave the file empty

def execute_command(command):
    logger.info(" ".join(command))
    stdout, stderr = subprocess.Popen(command, stdout = subprocess.PIPE, stderr = subprocess.PIPE).communicate()
    if stderr != "": logging.warn(stderr)
    return stdout

def post_files(files, destination):
    for f in files:
        execute_command(['rsync', '-a', f, '{}/{}'.format(destination, os.path.basename(f))])


def download_quality(run_info_path, destination, read_lengths):
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
                expected_cycle += 16
    return qcRunId


def collate_results(run_folder, results_folder, logger):
    """
    Collate pipeline logs and outputs into a small set of files.
    :param run_folder: Path to write collated files
    :param logger: miseq_logging object
    :return:
    """
    logger.info("Collating *.coverage.log files")
    miseq_logging.collate_logs(run_folder, "coverage.log", "coverage.log")

    logger.info("Collating *.mapping.log files")
    miseq_logging.collate_logs(run_folder, "mapping.log", "mapping.log")

    logger.info("Collating *.sam2aln.*.log files")
    miseq_logging.collate_logs(run_folder, "sam2aln.log", "sam2aln.log")

    logger.info("Collating *.g2p.log files")
    miseq_logging.collate_logs(run_folder, "g2p.log", "g2p.log")

    logger.info("Collating aln2counts.log files")
    miseq_logging.collate_logs(run_folder,
                               "aln2counts.log",
                               "aln2counts.log")

    logger.info("Collating aln2nuc.log files")
    miseq_logging.collate_logs(run_folder, "aln2nuc.log", "aln2nuc.log")

    files_to_collate = (('amino_frequencies.csv', '*.amino.csv'),
                        ('collated_conseqs.csv', '*.conseq.csv'),
                        ('coverage_scores.csv', None),
                        ('failed_align.csv', None),
                        ('failed_read.csv', None),
                        ('g2p.csv', None),
                        ('conseq_ins.csv', None),
                        ('coord_ins.csv', None),
                        ('nucleotide_frequencies.csv', '*.nuc.csv'),
                        ('nuc_variants.csv', None),
                        ('collated_counts.csv', '*.remap_counts.csv'))

    for target_file, pattern in files_to_collate:
        logger.info("Collating {}".format(target_file))
        if pattern is None:
            pattern = '*.' + target_file
        collate_labeled_files(os.path.join(run_folder, pattern),
                              os.path.join(results_folder, target_file))



#####################################################################################

KiveAPI.SERVER_URL = settings.kive_server_url
kive = KiveAPI(settings.kive_user, settings.kive_password, verify=False)

# retrieve Pipeline object based on version
pipeline = kive.get_pipeline(settings.pipeline_version_kive_id)

# retrieve quality.csv compound data type
quality_cdt = kive.get_cdt(settings.quality_cdt_kive_id)

processed_runs = set()
logger = None
failure_message = None


def check_kive(scheduler, runs):
    """
    :param scheduler: instance of Python sched
    :param runs: a List of Kive pipeline run objects
    :return:
    """
    for run in runs:
        logger.info('{} {}%'.format(run.get_status(), run.get_progress_percent()))
    if not all(map(lambda x: x.is_complete(), runs)):
        # reschedule a check on runs at subsequent time
        scheduler.enter(5, 1, check_kive, (scheduler, runs,))


## main loop
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
        
    # flag indicates that Illumina MiseqReporter has completed pre-processing, files available on NAS
    runs = glob(settings.rawdata_mount + 'MiSeq/runs/*/{}'.format(settings.NEEDS_PROCESSING))
    #runs = glob(rawdata_mount + 'MiSeq/runs/131119_M01841_0041_000000000-A5EPY/{}'.format(NEEDS_PROCESSING))

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
        else:
            # running in development mode - do all runs even if already processed
            # note that results will not be uploaded
            if run in processed_runs:
                continue

        runs_needing_processing.append(run)

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
    run_name = root.split('/')[-2]
    run_folder = settings.home+run_name
    results_folder = os.path.join(run_folder, 'results')

    # Make folder on the cluster for intermediary files and outputs
    if not os.path.exists(run_folder):
        os.mkdir(run_folder)
    if not os.path.exists(results_folder):
        os.mkdir(results_folder)

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


    # transfer SampleSheet.csv
    remote_file = curr_run.replace(settings.NEEDS_PROCESSING, 'SampleSheet.csv')
    local_file = run_folder + '/SampleSheet.csv'
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


    gz_files = glob(root + 'Data/Intensities/BaseCalls/*R?_001.fastq.gz')
    if not gz_files:
        failure_message = mark_run_as_disabled(root, "No data files found")
        continue

    # use Kive-API to transfer fastq.gz files
    fastqs = {}
    R1_files = filter(lambda x: '_R1_' in os.path.basename(x), gz_files)
    for R1_file in R1_files:
        R2_file = R1_file.replace('_R1_', '_R2_')
        if R2_file not in gz_files:
            logger.info("ERROR: Detected an unpaired R1 FASTQ file: {}".format(R1_file))
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
            logger.info("{} reads failed to demultiplex in {} (removing file)".format(failed_demultiplexing, filename))
            continue

        logger.info(filename)

        R1_obj = kive.add_dataset(name=filename,
                                  description='',
                                  handle=open(R1_file, 'rb'),
                                  cdt=None,
                                  users=None,
                                  groups=['Everyone'])

        R2_obj = kive.add_dataset(name=os.path.basename(R2_file),
                                  description='',
                                  handle=open(R2_file, 'rb'),
                                  cdt=None,
                                  users=None,
                                  groups=['Everyone'])
        fastqs.update({(sample, snum): (R1_obj, R2_obj)})
        break  #FIXME: for debugging - append only one sample

    # generate quality.csv
    try:
        quality_csv = os.path.join(settings.home, run_name, 'quality.csv')
        download_quality(run_info_path=os.path.join(root, 'RunInfo.xml'),
                         destination=quality_csv,
                         read_lengths=read_lengths)
    except StandardError as e:
        failure_message = mark_run_as_disabled(root,
                                               "Quality could not be downloaded.",
                                               exc_info=True)
        continue

    # transfer quality.csv with Kive API
    quality_input = kive.add_dataset(name='quality.csv',
                                     description='phiX174 quality scores per tile and cycle for run %s' % (run_name,),
                                     handle=open(quality_csv, 'rU'),
                                     cdt=quality_cdt,
                                     users=None,
                                     groups=['Everyone'])

    # Standard out/error concatenates to the log
    logger.info("Launching pipeline for %s%s", settings.home, run_name)
    kive_runs = []
    try:
        monitor_path = os.path.abspath(os.path.dirname(__file__))
        environment = dict(os.environ)
        old_path = environment.get('PYTHONPATH', '')
        environment['PYTHONPATH'] = os.pathsep.join((old_path, monitor_path))
        #subprocess.check_call([os.path.join(monitor_path, 'micall/monitor/run_processor.py'),
        #                       home+run_name,
        #                       '--clean' if do_cleanup else ''],
        #                      env=environment)

        # push all samples into the queue
        for key, (fastq1, fastq2) in fastqs.iteritems():
            status = kive.run_pipeline(pipeline=pipeline, inputs=[fastq1, fastq2, quality_input])
            kive_runs.append(status)

        # initialize progress monitoring
        sc = sched.scheduler(time.time, time.sleep)
        sc.enter(5, 1, check_kive, (sc, kive_runs))
        sc.run()  # exits when all runs are complete

        logger.info("===== {} successfully processed! =====".format(run_name))
    except Exception as e:
        failure_message = mark_run_as_disabled(
            root,
            "MISEQ_MONITOR.py failed: '{}'".format(e),
            exc_info=True)
        continue


    # Retrieve pipeline output files from Kive
    for kive_run in kive_runs:
        for dataset in kive_run.get_results():
            with open(os.path.join(run_folder, dataset.filename), 'wb') as handle:
                dataset.download(handle)

    # Collate pipeline logs and outputs to location where following code expects
    collate_results(run_folder, results_folder, logger)


    if not settings.production:
        processed_runs.add(curr_run)
    else:
        try:
            # Determine output paths
            result_path = curr_run.replace(settings.NEEDS_PROCESSING, 'Results')
            result_path_final = '{}/version_{}'.format(result_path, settings.pipeline_version)
    
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
                full_pattern = os.path.join(settings.home, run_name, pattern)
                if not os.path.isdir(destination_path):
                    os.mkdir(destination_path)
                post_files(glob(full_pattern), destination_path)


            untar_path = os.path.join(result_path_final, 'untar')
            coverage_source_path = os.path.join(untar_path, 'coverage_maps')
            coverage_dest_path = os.path.join(result_path_final, 'coverage_maps')
            os.mkdir(untar_path)
            os.mkdir(coverage_source_path)
            os.mkdir(coverage_dest_path)

            for tar_path in glob(settings.home + run_name + '/*.coverage_maps.tar'):
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

    logging.shutdown()  #FIXME: should this not be outside the loop?
    logger = None
