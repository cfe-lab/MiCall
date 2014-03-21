"""
MISEQ_MONITOR.py
1) For runs flagged 'needsprocessing' (And no error flag), unzip fastqs to local disk
2) Process the run by clling MISEQ_PIPELINE.py
3) Upload results back to macdatafile
"""

pipeline_version = '5.3'

import logging, miseq_logging, miseqUtils, os, subprocess, sys, time
from settings import *
from glob import glob

if sys.version_info[:2] != (2, 7):
    raise Exception("Python 2.7 not detected")

def mark_run_as_disabled(curr_run):
    open(curr_run.replace(NEEDS_PROCESSING, ERROR_PROCESSING), 'w').close()

def execute_command(command):
    logger.info(" ".join(command))
    stdout, stderr = subprocess.Popen(command, stdout = subprocess.PIPE, stderr = subprocess.PIPE).communicate()
    if stderr != "": logging.warn(stderr)
    return stdout

def post_files(files, destination):
    for file in files:
        execute_command(['rsync', '-a', file, '{}/{}'.format(destination, os.path.basename(file))])


# Process runs flagged for processing not already processed by this version of the pipeline
while True:

    runs = glob(macdatafile_mount + 'MiSeq/runs/*/{}'.format(NEEDS_PROCESSING))

	# OVERRIDE
    runs = glob(macdatafile_mount + 'MiSeq/runs/140212_M01841_0056_000000000-A64G4/{}'.format(NEEDS_PROCESSING))

    runs_needing_processing = []
    for run in runs:
        result_path = '{}/version_{}'.format(run.replace(NEEDS_PROCESSING, 'Results'), pipeline_version)

        if os.path.exists(run.replace(NEEDS_PROCESSING, ERROR_PROCESSING)):
            continue

        if not os.path.exists(result_path):
            runs_needing_processing.append(run)

    if not runs_needing_processing:
        logging.info('No runs need processing')
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
    try:
        logger = miseq_logging.init_logging(log_file, file_log_level=logging.DEBUG, console_log_level=logging.INFO)
        logger.info('===== Processing {} with pipeline version {} ====='.format(root, pipeline_version))
    except Exception as e:
        raise Exception("Couldn't setup logging (init_logging() threw exception '{}') - HALTING NOW!".format(str(e)))

    # SampleSheet.csv needed to determine the mode and sample T-primer
    remote_file = curr_run.replace(NEEDS_PROCESSING, 'SampleSheet.csv')
    local_file = home + run_name + '/SampleSheet.csv'
    execute_command(['rsync', '-a', remote_file, local_file])

    try:
        with open(local_file, 'rU') as sample_sheet:
            run_info = miseqUtils.sampleSheetParser(sample_sheet)
            mode = run_info['Description']
    except Exception as e:
        logger.error("Exception thrown while parsing sample sheet: '{}' - skipping run and marking with flag '{}'".format(str(e),ERROR_PROCESSING))
        mark_run_as_disabled(curr_run)
        continue

    if mode not in ['Nextera', 'Amplicon']:
        logger.error("{} not a valid mode: skipping run and marking with {}".format(mode, ERROR_PROCESSING))
        mark_run_as_disabled(curr_run)
        continue

    # Copy fastq.gz files to the cluster and unzip them
    for gz_file in glob(root+'Data/Intensities/BaseCalls/*.fastq.gz'):
        filename = os.path.basename(gz_file)

        # Report number of reads failing to demultiplex to the log
        if filename.startswith('Undetermined'):
            #if filename.endswith('_L001_R2_001.fastq.gz'):
                # the second file will have exactly the same number of lines
            #    continue
            # do word count directly on stream redirected from gunzip
            #p1 = subprocess.Popen(['gunzip', '-c', gz_file], stdout=subprocess.PIPE)
            #p2 = subprocess.Popen(['wc', '-l'], stdin=p1.stdout, stdout=subprocess.PIPE)
            #output = p2.communicate()[0]
            #failed_demultiplexing = output.strip(' \n')
            #logger.info("{} reads failed to demultiplex in {} (removing file)".format(failed_demultiplexing, filename))
            continue

        local_file = home + run_name + '/' + filename
        # if local copy of unzipped fastq exists, skip
        if os.path.exists(local_file.replace('.gz', '')):
            continue

        execute_command(['rsync', '-a', gz_file, local_file])
        execute_command(['gunzip', '-f', local_file])

    # Store output of MISEQ_PIPELINE.py in a log + poll continuously to display output to console
    pipeline_log_path = home + run_name + '/MISEQ_PIPELINE_OUTPUT.log'
    with open(pipeline_log_path, "wb") as PIPELINE_log:

        # Standard out/error concatenates to the log
        command = ['python2.7', '-u', 'MISEQ_PIPELINE.py', home+run_name]
        p = subprocess.Popen(command, stdout = PIPELINE_log, stderr = PIPELINE_log)
        logger.info(" ".join(command))

        # Poll the log of MISEQ_PIPELINE.py and display output to console as it appears
        with open(pipeline_log_path, 'rb') as cursor:
            while True:
                if p.poll() == 0:
                    PIPELINE_log.flush()
                    sys.stdout.write(cursor.read())
                    break
                time.sleep(1)
                PIPELINE_log.flush()
                sys.stdout.write(cursor.read())

    # Determine output paths
    result_path = curr_run.replace(NEEDS_PROCESSING, 'Results')
    result_path_final = '{}/version_{}'.format(result_path, pipeline_version)
    log_path = '{}/logs'.format(result_path_final)
    counts_path = '{}/counts'.format(log_path)
    coverage_maps_path = '{}/coverage_maps'.format(result_path_final)

    # Create sub-folders if needed
    for path in [result_path, result_path_final, log_path, counts_path, coverage_maps_path]:
        if not os.path.exists(path): os.mkdir(path)

    # Post files to each appropriate sub-folder
    logging.info("Posting results to {}".format(result_path))
    if mode == 'Amplicon':
        v3_path = '{}/v3_tropism'.format(result_path_final)
        if not os.path.exists(v3_path): os.mkdir(v3_path)
        post_files(glob(home + run_name + '/*.v3prot'), v3_path)
        post_files(glob(home + run_name + '/v3_tropism_summary.txt'), v3_path)

    #post_files(glob(home + run_name + '/*.counts'), counts_path)
    post_files(glob(home + run_name + '/*.log'), log_path)
    post_files(glob(home + run_name + '/collated_counts.csv'), counts_path)
    post_files(glob(home + run_name + '/collated_conseqs.csv'), result_path_final)
    post_files(glob(home + run_name + '/amino_frequencies.csv'), result_path_final)
    post_files(glob(home + run_name + '/nucleotide_frequencies.csv'), result_path_final)
    post_files(glob(home + run_name + '/coverage_maps/*.png'), coverage_maps_path)

    # Close the log and copy it to macdatafile
    logger.info("===== {} successfully completed! =====".format(run_name))
    logging.shutdown()
    execute_command(['rsync', '-a', log_file, '{}/{}'.format(log_path, os.path.basename(log_file))])
