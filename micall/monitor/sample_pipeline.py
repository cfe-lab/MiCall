#!/usr/bin/env python

import argparse
from glob import glob
import logging, os

from micall.utils.collate import collate_labeled_files
from fifo_scheduler import Job, Worker
from micall.core import miseq_logging
from micall.utils.sample_sheet_parser import sample_sheet_parser
from micall.settings import are_temp_folders_deleted
from mpi4py import MPI

def build_path(path):
    base_path = os.path.abspath(os.path.dirname(os.path.dirname(__file__)))
    return os.path.join(base_path, path)

def parseOptions(comm_world):
    parser = argparse.ArgumentParser(
        description='Process all the samples in a single run folder.')
    
    parser.add_argument('run_folder',
                        help='Path to sample fastq files and SampleSheet.csv')
    parser.add_argument('mode',
                        help='Amplicon or Nextera, default from sample sheet',
                        nargs=argparse.OPTIONAL)
    parser.add_argument('--phase', '-p',
                        help='Phase to execute: mapping, counting, plotting, or all',
                        default='all')
    
    # Rank 0 process parses the arguments, complains about errors, 
    # then broadcasts the parsed arguments (or None) to all other ranks.
    args = None
    try:
        if comm_world.Get_rank() == 0:
            args = parser.parse_args()
    finally:
        args = comm_world.bcast(args, root=0)
    
    if args is None:
        exit(0)
    return args

class SampleInfo(object):
    def __init__(self, fastq1):
        self.fastq1 = fastq1
        basename = os.path.basename(fastq1)
        dirname = os.path.dirname(fastq1)
        self.fastq2 = os.path.join(dirname, basename.replace('_R1_', '_R2_'))
        sample_name, sample_number = basename.split('_')[:2]
        self.key = sample_name + '_' + sample_number
        self.output_root = os.path.join(dirname, self.key)
        
def filter_quality(run_folder, worker):
    log_path = os.path.join(run_folder, 'quality.log')
    worker.run_job(Job(script=build_path('core/filter_quality.py'),
                       helpers=(build_path('settings.py'),
                                build_path('core/miseq_logging.py'),
                                build_path('projects.json')),
                       args=(os.path.join(run_folder, 'quality.csv'),
                             os.path.join(run_folder, 'bad_cycles.csv')),
                       stdout=log_path,
                       stderr=log_path))
    
def map_samples(run_folder, fastq_samples, worker):
    for sample_info in fastq_samples:
        log_path = "{}.censor.log".format(sample_info.output_root)
        worker.run_job(Job(script=build_path('core/censor_fastq.py'),
                           args=(sample_info.fastq1,
                                 os.path.join(run_folder, 'bad_cycles.csv'),
                                 sample_info.output_root + '.censored1.fastq'),
                           stdout=log_path,
                           stderr=log_path))
        worker.run_job(Job(script=build_path('core/censor_fastq.py'),
                           args=(sample_info.fastq2,
                                 os.path.join(run_folder, 'bad_cycles.csv'),
                                 sample_info.output_root + '.censored2.fastq'),
                           stdout=log_path,
                           stderr=log_path))
        log_path = "{}.mapping.log".format(sample_info.output_root)
        # Preliminary map
        worker.run_job(Job(script=build_path('core/prelim_map.py'),
                           helpers=(build_path('settings.py'),
                                    build_path('core/miseq_logging.py'),
                                    build_path('core/project_config.py'),
                                    build_path('projects.json')),
                           args=(sample_info.output_root + '.censored1.fastq',
                                 sample_info.output_root + '.censored2.fastq',
                                 sample_info.output_root + '.prelim.csv'),
                           stdout=log_path,
                           stderr=log_path))
    
        # Iterative re-mapping
        worker.run_job(Job(script=build_path('core/remap.py'),
                           helpers=(build_path('settings.py'),
                                    build_path('core/miseq_logging.py'),
                                    build_path('core/project_config.py'),
                                    build_path('projects.json')),
                           args=(sample_info.output_root + '.censored1.fastq', 
                                 sample_info.output_root + '.censored2.fastq',
                                 sample_info.output_root + '.prelim.csv',
                                 sample_info.output_root + '.remap.csv',
                                 sample_info.output_root + '.remap_counts.csv',
                                 sample_info.output_root + '.remap_conseq.csv',
                                 sample_info.output_root + '.unmapped1.fastq',
                                 sample_info.output_root + '.unmapped2.fastq'),
                           stdout=log_path,
                           stderr=log_path))

def count_samples(fastq_samples, worker, args):
    ########################
    ### Begin sam2aln
    
    
    for sample_info in fastq_samples:
        log_path = "{}.sam2aln.log".format(sample_info.output_root)
        worker.run_job(Job(script=build_path('core/sam2aln.py'),
                           helpers=(build_path('settings.py'), ),
                           args=(sample_info.output_root + '.remap.csv',
                                 sample_info.output_root + '.aligned.csv',
                                 sample_info.output_root + '.conseq_ins.csv',
                                 sample_info.output_root + '.failed_read.csv'),
                           stdout=log_path,
                           stderr=log_path))
    
    #######################
    ### Begin aln2counts
    
    for sample_info in fastq_samples:
        log_path = "{}.aln2counts.log".format(sample_info.output_root)
        worker.run_job(Job(script=build_path('core/aln2counts.py'),
                           helpers=(build_path('settings.py'),
                                    build_path('core/project_config.py'),
                                    build_path('projects.json'),
                                    build_path('core/miseq_logging.py')),
                           args=(sample_info.output_root + '.aligned.csv',
                                 sample_info.output_root + '.nuc.csv',
                                 sample_info.output_root + '.amino.csv',
                                 sample_info.output_root + '.coord_ins.csv',
                                 sample_info.output_root + '.conseq.csv',
                                 sample_info.output_root + '.failed_align.csv',
                                 sample_info.output_root + '.nuc_variants.csv'),
                           stdout=log_path,
                           stderr=log_path))
        
    for sample_info in fastq_samples:
        log_path = "{}.g2p.log".format(sample_info.output_root)
        worker.run_job(Job(script=build_path('g2p/sam_g2p.py'),
                           helpers=(build_path('core/sam2aln.py'),
                                    build_path('utils/translation.py'),
                                    build_path('g2p/g2p.matrix'),
                                    build_path('g2p/g2p_fpr.txt')),
                           args=(sample_info.output_root + '.remap.csv',
                                 sample_info.output_root + '.nuc.csv',
                                 sample_info.output_root + '.g2p.csv'),
                           stdout=log_path,
                           stderr=log_path))

    for sample_info in fastq_samples:
        log_path = "{}.coverage.log".format(sample_info.output_root)
        worker.run_job(Job(script=build_path('monitor/coverage_plots.R'),
                           helpers=(build_path('projects.json'), ),
                           args=(sample_info.output_root + '.amino.csv',
                                 sample_info.output_root + '.coverage_maps.tar',
                                 sample_info.output_root + '.coverage_scores.csv'),
                           stdout=log_path,
                           stderr=log_path))

            
def collate_results(fastq_samples, worker, args, logger):

    results_folder = os.path.join(args.run_folder, 'results')
    if not os.path.isdir(results_folder):
        os.mkdir(results_folder)
    
    logger.info("Collating *.coverage.log files")
    miseq_logging.collate_logs(args.run_folder, "coverage.log", "coverage.log")
    
    logger.info("Collating *.mapping.log files")
    miseq_logging.collate_logs(args.run_folder, "mapping.log", "mapping.log")
    
    logger.info("Collating *.sam2aln.*.log files")
    miseq_logging.collate_logs(args.run_folder, "sam2aln.log", "sam2aln.log")
    
    logger.info("Collating *.g2p.log files")
    miseq_logging.collate_logs(args.run_folder, "g2p.log", "g2p.log")
    
    logger.info("Collating aln2counts.log files")
    miseq_logging.collate_logs(args.run_folder,
                               "aln2counts.log",
                               "aln2counts.log")
    
    logger.info("Collating aln2nuc.log files")
    miseq_logging.collate_logs(args.run_folder, "aln2nuc.log", "aln2nuc.log")
    
    logger.info("Collating censor.log files")
    miseq_logging.collate_logs(args.run_folder, "censor.log", "censor.log")
    
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
        collate_labeled_files(os.path.join(args.run_folder, pattern),
                              os.path.join(results_folder, target_file))



def main():
    comm = MPI.COMM_WORLD  # @UndefinedVariable
    process_rank = comm.Get_rank()
    process_count = comm.Get_size()
    
    args = parseOptions(comm)
    log_file = "{}/pipeline{}.log".format(args.run_folder, process_rank)
    logger = miseq_logging.init_logging(log_file,
                                        file_log_level=logging.DEBUG,
                                        console_log_level=logging.INFO)
    logger.info('Start processing run %s, rank %d',
                args.run_folder,
                process_rank)
    
    if args.mode is not None:
        run_info = None
    else:
        with open(args.run_folder+'/SampleSheet.csv', 'rU') as sample_sheet:
            logger.debug("sample_sheet_parser({})".format(sample_sheet))
            run_info = sample_sheet_parser(sample_sheet)
            args.mode = run_info['Description']

    fastq_samples = []
    fastq_files = glob(args.run_folder + '/*_R1_001.fastq')
    for i, fastq in enumerate(fastq_files):
        if i % process_count != process_rank:
            # skip samples that are assigned to other worker processes
            continue
        
        sample_info = SampleInfo(fastq)
    
        # verify this sample is in SampleSheet.csv
        if run_info and not run_info['Data'].has_key(sample_info.key):
            logger.error(
                '{} not in SampleSheet.csv - cannot map this sample'.format(
                    sample_info.key))
            continue
    
        fastq_samples.append(sample_info)
    
    def launch_callback(command):
        logger.info("Launching {!r}".format(command))

    worker = Worker(launch_callback=launch_callback,
                    working_path=args.run_folder,
                    are_temp_folders_deleted=are_temp_folders_deleted,
                    logger=logger)
    
    if args.phase in ('filter', 'all'):
        filter_quality(args.run_folder, worker)
    
    if args.phase in ('mapping', 'all'):
        map_samples(args.run_folder, fastq_samples, worker)
    
    if args.phase in ('counting', 'all'):
        count_samples(fastq_samples, worker, args)
    
    if args.phase in ('summarizing', 'all') and process_rank == 0:
        collate_results(fastq_samples, worker, args, logger)

    # FIXME: this log message gets sent before workers start
    logger.info('Finish processing run %s, rank %d',
                args.run_folder,
                process_rank)

    

###############################
### (optional) filter for cross-contamination
#TODO: Once the rest of the pipeline stabilizes, try filtering again
# if filter_cross_contaminants:
#     logger.info('Filtering for cross-contamination')
# 
#     for qcut in sam2aln_q_cutoffs:
#         command = 'python2.7 FILTER_CONTAMINANTS.py %s %d %d' % (root, qcut, bowtie_threads)
#         log_path = root + '/filtering.log'
#         single_thread_factory.queue_work(command, log_path, log_path)
#         # yields *.clean.aln and *.contam.aln
# 
#     single_thread_factory.wait()
# 
#     # re-do g2p scoring
#     for env_aln_file in glob(root + '/*.HIV1B-env.*.clean.aln'):
#         command = "python2.7 STEP_3_G2P.py {} {}".format(env_aln_file, g2p_alignment_cutoff)
#         log_path = "{}.g2p.log".format(env_aln_file)
#         single_thread_factory.queue_work(command, log_path, log_path)
#     single_thread_factory.wait()
#     logger.info("Collating *.g2p.log files")
#     miseq_logging.collate_logs(root, "g2p.log", "g2p.log")
# 
#     with open("{}/v3_tropism_summary_clean.csv".format(root), 'wb') as summary_file:
#         summary_file.write("sample,q_cutoff,fpr_cutoff,min_count,total_x4,total_reads,proportion_x4\n")
#         for f in glob(root + '/*.clean.v3prot'):
#             prefix, gene, sam2aln_q_cutoff = f.split('.')[:3]
# 
#             # Determine proportion x4 under different FPR cutoffs / min counts
#             for fpr_cutoff in g2p_fpr_cutoffs:
#                 for mincount in v3_mincounts:
#                     try:
#                         sample = prefix.split('/')[-1]
#                         proportion_x4, total_x4_count, total_count = prop_x4(f, fpr_cutoff, mincount)
#                         summary_file.write("{},{},{},{},{},{},{:.3}\n".format(
#                             sample, sam2aln_q_cutoff, fpr_cutoff, mincount, total_x4_count,
#                             total_count, proportion_x4
#                         ))
#                     except Exception as e:
#                         logger.warn("miseqUtils.prop_x4() threw exception '{}'".format(str(e)))
# 
#     # regenerate counts
#     for aln_file in glob(root + '/*.clean.aln'):
#         mixture_cutoffs = ",".join(map(str,conseq_mixture_cutoffs))
#         command = "python2.7 STEP_4_CSF2COUNTS.py {} {} {} {}".format(
#             aln_file, mode, mixture_cutoffs, final_alignment_ref_path
#         )
#         log_path = "{}.aln2counts.log".format(aln_file)
#         single_thread_factory.queue_work(command, log_path, log_path)
# 
#     single_thread_factory.wait()
# 
#     collated_clean_amino_freqs_path = "{}/amino_cleaned_frequencies.csv".format(root)
#     command = ["python2.7","generate_coverage_plots.py",collated_clean_amino_freqs_path,"{}/coverage_maps".format(root)]
#     logger.info(" ".join(command))
#     subprocess.call(command)

if __name__ == '__main__':
    main()
