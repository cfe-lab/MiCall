#!/usr/bin/env python

import argparse
from glob import glob
import logging, os

from collate import collate_conseqs, collate_labeled_files
from fifo_scheduler import Job, Worker
import miseq_logging
from sample_sheet_parser import sample_sheet_parser
from settings import final_alignment_ref_path, final_nuc_align_ref_path, \
    are_temp_folders_deleted, mapping_ref_path, base_path
from mpi4py import MPI

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
        
def map_samples(fastq_samples, worker):
    # Preliminary map
    for sample_info in fastq_samples:
        log_path = "{}.mapping.log".format(sample_info.output_root)
        worker.run_job(Job(script=base_path + 'prelim_map.py',
                           helpers=(base_path + 'settings.py',
                                    base_path + 'miseq_logging.py',
                                    mapping_ref_path + '.fasta'),
                           args=(sample_info.fastq1,
                                 sample_info.fastq2,
                                 sample_info.output_root + '.prelim.csv'),
                           stdout=log_path,
                           stderr=log_path))
    
    # Iterative re-mapping
    for sample_info in fastq_samples:
        log_path = "{}.mapping.log".format(sample_info.output_root)
        worker.run_job(Job(script=base_path + 'remap.py',
                           helpers=(base_path + 'settings.py',
                                    base_path + 'miseq_logging.py',
                                    mapping_ref_path + '.fasta'),
                           args=(sample_info.fastq1, 
                                 sample_info.fastq2,
                                 sample_info.output_root + '.prelim.csv',
                                 sample_info.output_root + '.remap.csv',
                                 sample_info.output_root + '.remap_counts.csv',
                                 sample_info.output_root + '.remap_conseq.csv'),
                           stdout=log_path,
                           stderr=log_path))

def count_samples(fastq_samples, worker, args):
    ########################
    ### Begin sam2csf
    
    
    for sample_info in fastq_samples:
        log_path = "{}.sam2csf.log".format(sample_info.output_root)
        worker.run_job(Job(script=base_path + 'sam2csf.py',
                           helpers=(base_path + 'settings.py', ),
                           args=(sample_info.output_root + '.remap.csv',
                                 sample_info.output_root + '.aligned.csv',
                                 sample_info.output_root + '.failed.csv'),
                           stdout=log_path,
                           stderr=log_path))
    
    #######################
    ### Begin csf2counts
    
    for sample_info in fastq_samples:
        log_path = "{}.csf2counts.log".format(sample_info.output_root)
        worker.run_job(Job(script=base_path + 'csf2counts.py',
                           helpers=(base_path + 'settings.py',
                                    final_alignment_ref_path,
                                    base_path + 'miseq_logging.py',
                                    base_path + 'hyphyAlign.py'),
                           args=(sample_info.output_root + '.aligned.csv',
                                 sample_info.output_root + '.remap_conseq.csv',
                                 sample_info.output_root + '.nuc.csv',
                                 sample_info.output_root + '.amino.csv',
                                 sample_info.output_root + '.indels.csv',
                                 sample_info.output_root + '.conseq.csv'),
                           stdout=log_path,
                           stderr=log_path))
        
    for sample_info in fastq_samples:
        log_path = "{}.csf2nuc.log".format(sample_info.output_root)
        worker.run_job(Job(script=base_path + 'csf2nuc.py',
                           helpers=(base_path + 'settings.py',
                                    final_nuc_align_ref_path,
                                    base_path + 'hyphyAlign.py'),
                           args=(sample_info.output_root + '.aligned.csv',
                                 sample_info.output_root + '.nuc_variants.csv'),
                           stdout=log_path,
                           stderr=log_path))
        
    for sample_info in fastq_samples:
        log_path = "{}.coverage.log".format(sample_info.output_root)
        worker.run_job(Job(script=base_path + 'coverage_plots.R',
                           helpers=(base_path + 'key_positions.csv', ),
                           args=(sample_info.output_root + '.amino.csv',
                                 sample_info.output_root + '.coverage_maps.tar',
                                 sample_info.output_root + '.coverage_scores.csv'),
                           stdout=log_path,
                           stderr=log_path))

def collate_results(fastq_samples, worker, args, logger):
    
    #TODO: Move fasta_to_g2p back to count_samples().
    # We have to make it thread safe before it can go there.
    # Running it single threaded adds roughly an hour to a typical run.
    # Remove fastq_samples and worker parameters when it moves.
    
    ###############################
    ### Begin g2p (For Amplicon covering HIV-1 env only!)
    if args.mode == 'Amplicon':
        # Compute g2p V3 tropism scores from HIV1B-env csf files and store in v3prot files
        for sample_info in fastq_samples:
            log_path = "{}.g2p.log".format(sample_info.output_root)
            worker.run_job(Job(script=base_path + 'sam_g2p.sh',
                               helpers=(base_path + 'sam_g2p.rb',
                                        base_path + 'pssm_lib.rb',
                                        base_path + 'alignment.so',
                                        base_path + 'g2p.matrix',
                                        base_path + 'g2p_fpr.txt'),
                               args=(sample_info.output_root + '.remap.csv',
                                     sample_info.output_root + '.g2p.csv'),
                               stdout=log_path,
                               stderr=log_path))
    
    logger.info("Collating *.mapping.log files")
    miseq_logging.collate_logs(args.run_folder, "mapping.log", "mapping.log")
    
    logger.info("Collating *.sam2csf.*.log files")
    miseq_logging.collate_logs(args.run_folder, "sam2csf.log", "sam2csf.log")
    
    logger.info("Collating *.g2p.log files")
    miseq_logging.collate_logs(args.run_folder, "g2p.log", "g2p.log")
    
    logger.info("Collating csf2counts.log files")
    miseq_logging.collate_logs(args.run_folder,
                               "csf2counts.log",
                               "csf2counts.log")
    
    logger.info("Collating csf2nuc.log files")
    miseq_logging.collate_logs(args.run_folder, "csf2nuc.log", "csf2nuc.log")
    
    collated_conseq_path = "{}/collated_conseqs.csv".format(args.run_folder)
    logger.info("collate_conseqs({},{})".format(args.run_folder,
                                                collated_conseq_path))
    collate_conseqs(args.run_folder, collated_conseq_path)
    
    files_to_collate = (('coverage_scores.csv', None),
                        ('collated_counts.csv', '*.remap_counts.csv'),
                        ('amino_frequencies.csv', '*.amino.csv'),
                        ('nucleotide_frequencies.csv', '*.nuc.csv'))
    
    for target_file, pattern in files_to_collate:
        logger.info("Collating {}".format(target_file))
        if pattern is None:
            pattern = '*.' + target_file
        collate_labeled_files(os.path.join(args.run_folder, pattern),
                              os.path.join(args.run_folder, target_file))

def main():
    comm = MPI.COMM_WORLD
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
                    are_temp_folders_deleted=are_temp_folders_deleted)
    
    if args.phase in ('mapping', 'all'):
        map_samples(fastq_samples, worker)
    
    if args.phase in ('counting', 'all'):
        count_samples(fastq_samples, worker, args)
    
    if args.phase in ('summarizing', 'all') and process_rank == 0:
        collate_results(fastq_samples, worker, args, logger)
    
    logger.info('Finish processing run %s, rank %d',
                args.run_folder,
                process_rank)
    

###############################
### (optional) filter for cross-contamination
#TODO: Once the rest of the pipeline stabilizes, try filtering again
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

if __name__ == '__main__':
    main()
