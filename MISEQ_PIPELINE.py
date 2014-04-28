from settings import *

import logging, miseq_logging, miseq_modules, miseqUtils, os, subprocess, sys, time
sys.path.append(path_to_fifo_scheduler)
from fifo_scheduler import Factory
from glob import glob

# Mapping factory suitable for multi-thread jobs (4 processes * 8 threads / job = 32 cores allocated)
mapping_factory = Factory(mapping_factory_resources)
single_thread_factory = Factory(single_thread_resources)

root = sys.argv[1]			# MONITOR parameter: Location of fastq files to process

# Logging parameters
log_file = "{}/pipeline_output.log".format(root)
logger = miseq_logging.init_logging(log_file, file_log_level=logging.DEBUG, console_log_level=logging.INFO)

def factory_barrier(my_factory):
    """Wait until factory completes queued jobs"""
    while True:
        processes_spawned = my_factory.supervise()
        if processes_spawned:
            for p, command in processes_spawned:
                logger.info("pID {}: {}".format(p.pid, command))
        if my_factory.completely_idle(): break
        time.sleep(1)
    return

if len(sys.argv) == 3:
    mode = sys.argv[2]
    run_info = None
else:
    with open(root+'/SampleSheet.csv', 'rU') as sample_sheet:
        logger.debug("sampleSheetParser({})".format(sample_sheet))
        run_info = miseqUtils.sampleSheetParser(sample_sheet)
        mode = run_info['Description']

### Begin Mapping
fastq_files = glob(root + '/*R1*.fastq')
fastq_files = [f for f in fastq_files if not f.endswith('.Tcontaminants.fastq')]
for fastq in fastq_files:
    fastq_filename = os.path.basename(fastq)
    sample_name, sample_number = fastq_filename.split('_')[:2]
    key = sample_name + '_' + sample_number
    if run_info and not run_info['Data'].has_key(key):
        logger.error('{} not in SampleSheet.csv - cannot initiate mapping for this sample'.format(key))
        continue
    is_t_primer = run_info['Data'][key]['is_T_primer'] if run_info else '0'
    command = "python2.7 STEP_1_MAPPING.py {} {} {} {} {} {} {} {}".format(mapping_ref_path,
            fastq, consensus_q_cutoff, mode, is_t_primer, min_mapping_efficiency, max_remaps, bowtie_threads)
    log_path = "{}.mapping.log".format(fastq)
    queue_request = mapping_factory.queue_work(command, log_path, log_path)
    if queue_request:
        p, command = queue_request
        logger.info("pID {}: {}".format(p.pid, command))
factory_barrier(mapping_factory)
logger.info("Collating *.mapping.log files")
miseq_logging.collate_logs(root, "mapping.log", "mapping.log")

### Begin sam2csf
for file in glob(root + '/*.remap.sam'):
    filename = file.split('/')[-1]
    # Generate csf with different q cutoff censoring rules
    for qcut in sam2csf_q_cutoffs:
        # CSFs are sorted by read prevalence for Amplicon and left-offset for Nextera
        command = "python2.7 STEP_2_SAM2CSF.py {} {} {} {} {}".format(file, qcut, read_mapping_cutoff, mode, max_prop_N)
        log_path = "{}.sam2csf.{}.log".format(file, qcut)
        queue_request = single_thread_factory.queue_work(command, log_path, log_path)
        if queue_request:
            p, command = queue_request
            logger.info("pID {}: {}".format(p.pid, command))
factory_barrier(single_thread_factory)
logger.info("Collating *.sam2csf.*.log files")
miseq_logging.collate_logs(root, "sam2csf.*.log", "sam2csf.log")



### Begin csf2counts
for csf_file in glob(root + '/*.csf'):
    if csf_file.endswith('.clean.csf') or csf_file.endswith('.contam.csf'):
        # in case this has been re-run
        continue
    # Determine nucleotide/amino counts, along with the consensus, in HXB2/H77 space
    mixture_cutoffs = ",".join(map(str,conseq_mixture_cutoffs))
    command = "python2.7 STEP_4_CSF2COUNTS.py {} {} {} {}".format(csf_file, mode, mixture_cutoffs,
                                                                  final_alignment_ref_path)
    log_path = "{}.csf2counts.log".format(csf_file)
    queue_request = single_thread_factory.queue_work(command, log_path, log_path)
    if queue_request:
        p, command = queue_request
        logger.info("pID {}: {}".format(p.pid, command))

factory_barrier(single_thread_factory)


### Begin cross-contamination filter
logger.info('Filtering for cross-contamination')

for qcut in sam2csf_q_cutoffs:
    command = 'python2.7 FILTER_CONTAMINANTS.py %s %d %d' % (root, qcut, bowtie_threads)
    log_path = root + '/filtering.log'
    queue_request = single_thread_factory.queue_work(command, log_path, log_path)
    if queue_request:
        p, command = queue_request
        logger.info('pID {}: {}'.format(p.pid, command))

factory_barrier(single_thread_factory)



### Begin g2p (For Amplicon)
if mode == 'Amplicon':

    # Compute g2p V3 tropism scores from HIV1B-env csf files and store in v3prot files
    for env_csf_file in glob(root + '/*.HIV1B-env.*.clean.csf'):
        command = "python2.7 STEP_3_G2P.py {} {}".format(env_csf_file, g2p_alignment_cutoff)
        log_path = "{}.g2p.log".format(env_csf_file)
        queue_request = single_thread_factory.queue_work(command, log_path, log_path)
        if queue_request:
            p, command = queue_request
            logger.info("pID {}: {}".format(p.pid, command))
    factory_barrier(single_thread_factory)
    logger.info("Collating *.g2p.log files")
    miseq_logging.collate_logs(root, "g2p.log", "g2p.log")

    # Summarize the v3prot files in v3_tropism_summary.txt
    with open("{}/v3_tropism_summary.txt".format(root), 'wb') as summary_file:
        summary_file.write("sample,q_cutoff,fpr_cutoff,min_count,total_x4,total_reads,proportion_x4\n")
        for file in glob(root + '/*.v3prot'):
            prefix, gene, sam2csf_q_cutoff = file.split('.')[:3]

            # Determine proportion x4 under different FPR cutoffs / min counts
            for fpr_cutoff in g2p_fpr_cutoffs:
                for mincount in v3_mincounts:
                    try:
                        sample = prefix.split('/')[-1]
                        proportion_x4, total_x4_count, total_count = miseqUtils.prop_x4(file, fpr_cutoff, mincount)
                        summary_file.write("{},{},{},{},{},{},{:.3}\n".format(sample, sam2csf_q_cutoff,
                                fpr_cutoff, mincount, total_x4_count, total_count, proportion_x4))
                    except Exception as e:
                        logger.warn("miseqUtils.prop_x4() threw exception '{}'".format(str(e)))

### Repeat csf2counts
for csf_file in glob(root + '/*.clean.csf'):
    # Determine nucleotide/amino counts, along with the consensus, in HXB2/H77 space
    mixture_cutoffs = ",".join(map(str,conseq_mixture_cutoffs))
    command = "python2.7 STEP_4_CSF2COUNTS.py {} {} {} {}".format(csf_file,mode,mixture_cutoffs,final_alignment_ref_path)
    log_path = "{}.csf2counts.log".format(csf_file)
    queue_request = single_thread_factory.queue_work(command, log_path, log_path)
    if queue_request:
        p, command = queue_request
        logger.info("pID {}: {}".format(p.pid, command))


factory_barrier(single_thread_factory)
logger.info("Collating csf2counts.log files")
miseq_logging.collate_logs(root, "csf2counts.log", "csf2counts.log")


# Collate lines from csf2counts-derived frequency files
collated_amino_freqs_path = "{}/amino_frequencies.csv".format(root)
logger.info("collate_frequencies({},{},{})".format(root,collated_amino_freqs_path,"amino"))
miseq_modules.collate_frequencies(root,collated_amino_freqs_path, "amino")

collated_nuc_freqs_path = "{}/nucleotide_frequencies.csv".format(root)
logger.info("collate_frequencies({},{},{})".format(root,collated_nuc_freqs_path,"nuc"))
miseq_modules.collate_frequencies(root,collated_nuc_freqs_path, "nuc")

# Collate consensus sequences from csf2counts-derived conseq files
collated_conseq_path = "{}/collated_conseqs.csv".format(root)
logger.info("collate_conseqs({},{})".format(root,collated_conseq_path))
miseq_modules.collate_conseqs(root,collated_conseq_path)

# Collate read mapping efficiency from counts files
collated_counts_path = "{}/collated_counts.csv".format(root)
logger.info("collate_counts({},{})".format(root,collated_counts_path))
miseq_modules.collate_counts(root,collated_counts_path)

# Generate coverage maps (For clean and unclean data...)
command = ["python2.7","generate_coverage_plots.py",collated_amino_freqs_path,"{}/coverage_maps".format(root)]
logger.info(" ".join(command))
subprocess.call(command)

collated_clean_amino_freqs_path = "{}/amino_cleaned_frequencies.csv".format(root)
command = ["python2.7","generate_coverage_plots.py",collated_clean_amino_freqs_path,"{}/coverage_maps".format(root)]
logger.info(" ".join(command))
subprocess.call(command)

# Delete local files on the cluster that shouldn't be stored
logger.info("Deleting intermediary files")
for extension in file_extensions_to_delete:
    for file in glob("{}/*.{}".format(root, extension)):
        logging.debug("os.remove({})".format(file))
        os.remove(file)

logging.shutdown()
