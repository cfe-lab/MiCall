#!/usr/bin/env python

from glob import glob
import logging, os, sys

from collate import collate_frequencies, collate_conseqs, collate_counts
import miseq_logging
from sample_sheet_parser import sampleSheetParser
from settings import file_extensions_to_delete, file_extensions_to_keep, \
    final_alignment_ref_path, final_nuc_align_ref_path, \
    mapping_factory_resources, mapping_ref_path, \
    path_to_fifo_scheduler, production, single_thread_resources
    
sys.path.append(path_to_fifo_scheduler)
from fifo_scheduler import Factory, Job

def launch_callback(command):
    logger.info("Launching {!r}".format(command))

# Mapping factory suitable for multi-thread jobs (4 processes * 8 threads / job = 32 cores allocated)
mapping_factory = Factory(mapping_factory_resources,
                          launch_callback=launch_callback)
single_thread_factory = Factory(single_thread_resources,
                                launch_callback=launch_callback)
#TODO: Do we need to set the working directory for the factories?

root = sys.argv[1]			# MONITOR parameter: Location of fastq files to process

# Logging parameters
log_file = "{}/pipeline_output.log".format(root)
logger = miseq_logging.init_logging(log_file, file_log_level=logging.DEBUG, console_log_level=logging.INFO)

if len(sys.argv) == 3:
    mode = sys.argv[2]
    run_info = None
else:
    with open(root+'/SampleSheet.csv', 'rU') as sample_sheet:
        logger.debug("sampleSheetParser({})".format(sample_sheet))
        run_info = sampleSheetParser(sample_sheet)
        mode = run_info['Description']


class SampleInfo(object):
    def __init__(self, fastq1):
        self.fastq1 = fastq1
        basename = os.path.basename(fastq1)
        dirname = os.path.dirname(fastq1)
        self.fastq2 = os.path.join(dirname, basename.replace('_R1_', '_R2_'))
        sample_name, sample_number = basename.split('_')[:2]
        self.key = sample_name + '_' + sample_number
        self.output_root = os.path.join(dirname, self.key)
        
fastq_samples = []
fastq_files = glob(root + '/*_R1_001.fastq')
for fastq in fastq_files:
    sample_info = SampleInfo(fastq)

    # verify this sample is in SampleSheet.csv
    if run_info and not run_info['Data'].has_key(sample_info.key):
        logger.error(
            '{} not in SampleSheet.csv - cannot map this sample'.format(
                sample_info.key))
        continue

    fastq_samples.append(sample_info)
    
# Preliminary map
for sample_info in fastq_samples:
    log_path = "{}.mapping.log".format(sample_info.fastq1)
    mapping_factory.queue_job(Job(script='prelim_map.py',
                                  helpers=('settings.py',
                                           mapping_ref_path + '.fasta'),
                                  args=(sample_info.fastq1,
                                        sample_info.fastq2,
                                        sample_info.output_root + '.prelim.sam'),
                                  stdout=log_path,
                                  stderr=log_path))

mapping_factory.wait()

# Iterative re-mapping
for sample_info in fastq_samples:
    log_path = "{}.mapping.log".format(sample_info.fastq1)
    mapping_factory.queue_job(Job(script='remap.py',
                                  helpers=('settings.py',
                                           mapping_ref_path + '.fasta'),
                                  args=(sample_info.fastq1, 
                                        sample_info.fastq2,
                                        sample_info.output_root + '.prelim.sam',
                                        sample_info.output_root + '.output.csv',
                                        sample_info.output_root + '.counts',
                                        sample_info.output_root + '.conseq.csv'),
                                  stdout=log_path,
                                  stderr=log_path))

mapping_factory.wait()
logger.info("Collating *.mapping.log files")
miseq_logging.collate_logs(root, "mapping.log", "mapping.log")


########################
### Begin sam2csf

logger.info('Removing old CSF files')
old_csf_files = glob(root+'/*.csf')
for f in old_csf_files:
    os.remove(f)

for sample_info in fastq_samples:
    log_path = "{}.sam2csf.log".format(sample_info.fastq1)
    single_thread_factory.queue_job(Job(script='sam2csf.py',
                                        helpers=('settings.py', ),
                                        args=(sample_info.output_root + '.prelim.sam',
                                              sample_info.output_root + '.csf.csv'),
                                        stdout=log_path,
                                        stderr=log_path))
single_thread_factory.wait()
logger.info("Collating *.sam2csf.*.log files")
miseq_logging.collate_logs(root, "sam2csf.log", "sam2csf.log")

###############################
### Begin g2p (For Amplicon covering HIV-1 env only!)
if mode == 'Amplicon':
    # Compute g2p V3 tropism scores from HIV1B-env csf files and store in v3prot files
    for sample_info in fastq_samples:
        log_path = "{}.g2p.log".format(sample_info.fastq1)
        single_thread_factory.queue_job(Job(script='fasta_to_g2p.rb',
                                            args=(sample_info.output_root + '.csf.csv',
                                                  sample_info.output_root + '.g2p.csv'),
                                            stdout=log_path,
                                            stderr=log_path))

    single_thread_factory.wait()
    logger.info("Collating *.g2p.log files")
    miseq_logging.collate_logs(root, "g2p.log", "g2p.log")

#######################
### Begin csf2counts

for sample_info in fastq_samples:
    log_path = "{}.csf2counts.log".format(sample_info.fastq1)
    single_thread_factory.queue_job(Job(script='csf2counts.py',
                                        helpers=(final_alignment_ref_path, ),
                                        args=(sample_info.output_root + '.csf.csv',
                                              sample_info.output_root + '.conseq.csv',
                                              sample_info.output_root + '.nuc.freqs',
                                              sample_info.output_root + '.amino.freqs',
                                              sample_info.output_root + '.output_indels.csv',
                                              sample_info.output_root + '.conseq'),
                                        stdout=log_path,
                                        stderr=log_path))
    
# No dependency between csf2counts and csf2nuc, so no need to wait for factory
for sample_info in fastq_samples:
    log_path = "{}.csf2nuc.log".format(sample_info.fastq1)
    single_thread_factory.queue_job(Job(script='csf2nuc.py',
                                        helpers=(final_nuc_align_ref_path, ),
                                        args=(sample_info.output_root + '.csf.csv',
                                              sample_info.output_root + '.nuc2.csv'), #TODO: extension?
                                        stdout=log_path,
                                        stderr=log_path))

single_thread_factory.wait()

#########################
### Collate results files
logger.info("Collating csf2counts.log files")
miseq_logging.collate_logs(root, "csf2counts.log", "csf2counts.log")

collated_amino_freqs_path = "{}/amino_frequencies.csv".format(root)
logger.info("collate_frequencies({},{},{})".format(root, collated_amino_freqs_path, "amino"))
collate_frequencies(root, collated_amino_freqs_path, "amino")

collated_nuc_freqs_path = "{}/nucleotide_frequencies.csv".format(root)
logger.info("collate_frequencies({},{},{})".format(root, collated_nuc_freqs_path, "nuc"))
collate_frequencies(root, collated_nuc_freqs_path, "nuc")

# TODO: collate output from csf2nuc.py?

collated_conseq_path = "{}/collated_conseqs.csv".format(root)
logger.info("collate_conseqs({},{})".format(root, collated_conseq_path))
collate_conseqs(root, collated_conseq_path)

collated_counts_path = "{}/collated_counts.csv".format(root)
logger.info("collate_counts({},{})".format(root, collated_counts_path))
collate_counts(root, collated_counts_path)

coverage_map_path = "{}/coverage_maps.tar".format(root)
coverage_score_path = "{}/coverage_scores.csv".format(root)


# generate coverage plots
single_thread_factory.queue_job(Job(script='coverage_plots.R',
                                    helpers=('key_positions.csv', ),
                                    args=(collated_amino_freqs_path,
                                          coverage_map_path,
                                          coverage_score_path)))
single_thread_factory.wait()

###############################
### (optional) filter for cross-contamination
#TODO: Should we put this back in or remove it completely?
# if filter_cross_contaminants:
#     logger.info('Filtering for cross-contamination')
# 
#     for qcut in sam2csf_q_cutoffs:
#         command = 'python2.7 FILTER_CONTAMINANTS.py %s %d %d' % (root, qcut, bowtie_threads)
#         log_path = root + '/filtering.log'
#         single_thread_factory.queue_work(command, log_path, log_path)
#         # yields *.clean.csf and *.contam.csf
# 
#     single_thread_factory.wait()
# 
#     # re-do g2p scoring
#     for env_csf_file in glob(root + '/*.HIV1B-env.*.clean.csf'):
#         command = "python2.7 STEP_3_G2P.py {} {}".format(env_csf_file, g2p_alignment_cutoff)
#         log_path = "{}.g2p.log".format(env_csf_file)
#         single_thread_factory.queue_work(command, log_path, log_path)
#     single_thread_factory.wait()
#     logger.info("Collating *.g2p.log files")
#     miseq_logging.collate_logs(root, "g2p.log", "g2p.log")
# 
#     with open("{}/v3_tropism_summary_clean.csv".format(root), 'wb') as summary_file:
#         summary_file.write("sample,q_cutoff,fpr_cutoff,min_count,total_x4,total_reads,proportion_x4\n")
#         for f in glob(root + '/*.clean.v3prot'):
#             prefix, gene, sam2csf_q_cutoff = f.split('.')[:3]
# 
#             # Determine proportion x4 under different FPR cutoffs / min counts
#             for fpr_cutoff in g2p_fpr_cutoffs:
#                 for mincount in v3_mincounts:
#                     try:
#                         sample = prefix.split('/')[-1]
#                         proportion_x4, total_x4_count, total_count = prop_x4(f, fpr_cutoff, mincount)
#                         summary_file.write("{},{},{},{},{},{},{:.3}\n".format(
#                             sample, sam2csf_q_cutoff, fpr_cutoff, mincount, total_x4_count,
#                             total_count, proportion_x4
#                         ))
#                     except Exception as e:
#                         logger.warn("miseqUtils.prop_x4() threw exception '{}'".format(str(e)))
# 
#     # regenerate counts
#     for csf_file in glob(root + '/*.clean.csf'):
#         mixture_cutoffs = ",".join(map(str,conseq_mixture_cutoffs))
#         command = "python2.7 STEP_4_CSF2COUNTS.py {} {} {} {}".format(
#             csf_file, mode, mixture_cutoffs, final_alignment_ref_path
#         )
#         log_path = "{}.csf2counts.log".format(csf_file)
#         single_thread_factory.queue_work(command, log_path, log_path)
# 
#     single_thread_factory.wait()
# 
#     collated_clean_amino_freqs_path = "{}/amino_cleaned_frequencies.csv".format(root)
#     command = ["python2.7","generate_coverage_plots.py",collated_clean_amino_freqs_path,"{}/coverage_maps".format(root)]
#     logger.info(" ".join(command))
#     subprocess.call(command)


###########################
# Delete local files on the cluster that shouldn't be stored
if production:
    logger.info("Deleting intermediary files")
    # TODO: a bunch of these files are created in the temporary folders that now 
    # get deleted. Do we still need this?
    for extension in file_extensions_to_delete:
        for f in glob("{}/*.{}".format(root, extension)):
            if any([f.endswith(ext2) for ext2 in file_extensions_to_keep]):
                continue
            logging.debug("os.remove({})".format(f))
            os.remove(f)

logging.shutdown()
