from argparse import ArgumentParser
import csv
import errno
import fnmatch
import json
import logging
from operator import itemgetter
import os
import shutil
import subprocess

from micall.core.aln2counts import aln2counts
from micall.core.censor_fastq import censor
from micall.core.filter_quality import report_bad_cycles
from micall.core.remap import remap
from micall.core.prelim_map import prelim_map
from micall.core.sam2aln import sam2aln
from micall.monitor import phix_parser
from micall.g2p.sam_g2p import sam_g2p
from micall.g2p.pssm_lib import Pssm

logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s[%(levelname)s]%(name)s.%(funcName)s(): %(message)s')
logger = logging.getLogger('micall')


def parse_args():
    parser = ArgumentParser(description='Map FASTQ files to references.')
    parser.add_argument('data_path',
                        default='/data',
                        nargs='?',
                        help='data folder filled in by BaseSpace')
    return parser.parse_args()


def parse_json(json_file):
    """ Load JSON from an open file, and pull out the arguments for this run.

    :param json_file: an open file that contains JSON in the BaseSpace
    AppSession format.
    :return: an object with an attribute for each argument
    """
    class Args(object):
        pass
    args = Args()

    raw_args = json.load(json_file)
    arg_map = {item['Name']: item
               for item in raw_args['Properties']['Items']}

    args.name = arg_map['Input.app-session-name']['Content']
    run_content = arg_map['Input.run-id']['Content']
    args.run_id = run_content['Id']
    args.read_length1 = run_content['SequencingStats']['NumCyclesRead1']
    args.read_length2 = run_content['SequencingStats']['NumCyclesRead2']
    args.index_length1 = run_content['SequencingStats']['NumCyclesIndex1']
    args.index_length2 = run_content['SequencingStats']['NumCyclesIndex2']
    args.samples = sorted(arg_map['Input.sample-ids']['Items'],
                          key=itemgetter('Name'))
    args.project_id = arg_map['Input.project-id']['Content']['Id']

    return args


def censor_sample(filename, bad_cycles_path, censored_name):
    with open(filename, 'rb') as fastq, \
            open(bad_cycles_path, 'rU') as bad_cycles, \
            open(censored_name, 'w') as dest:
        censor(fastq, csv.DictReader(bad_cycles), dest)
        return dest.name


def process_sample(sample_info, project_id, data_path, pssm):
    scratch_path = os.path.join(data_path, 'scratch')
    sample_id = sample_info['Id']
    sample_name = sample_info['Name']
    sample_dir = os.path.join(data_path,
                              'input',
                              'samples',
                              sample_id,
                              'Data',
                              'Intensities',
                              'BaseCalls')
    if not os.path.exists(sample_dir):
        sample_dir = os.path.join(data_path,
                                  'input',
                                  'samples',
                                  sample_id)
    sample_path = None
    for root, _dirs, files in os.walk(sample_dir):
        sample_paths = fnmatch.filter(files, '*_R1_*')
        if sample_paths:
            sample_path = os.path.join(root, sample_paths[0])
            break
    if sample_path is None:
        raise RuntimeError('No R1 file found for sample id {}.'.format(sample_id))
    sample_path2 = sample_path.replace('_R1_', '_R2_')
    if not os.path.exists(sample_path2):
        raise RuntimeError('R2 file missing for sample id {}: {!r}.'.format(
            sample_id,
            sample_path2))
    logger.info('Processing sample %s: %s (%s).', sample_id, sample_name, sample_path)

    sample_out_path = os.path.join(data_path,
                                   'output',
                                   'appresults',
                                   project_id,
                                   sample_name)
    makedirs(sample_out_path)

    sample_scratch_path = os.path.join(scratch_path, sample_name)
    makedirs(sample_scratch_path)

    censored_path1 = os.path.join(sample_scratch_path, 'censored1.fastq')
    censor_sample(sample_path,
                  os.path.join(scratch_path, 'bad_cycles.csv'),
                  censored_path1)
    censored_path2 = os.path.join(sample_scratch_path, 'censored2.fastq')
    censor_sample(sample_path2,
                  os.path.join(scratch_path, 'bad_cycles.csv'),
                  censored_path2)

    logger.info('Running prelim_map.')
    with open(os.path.join(sample_scratch_path, 'prelim.csv'), 'wb') as prelim_csv:
        prelim_map(censored_path1,
                   censored_path2,
                   prelim_csv)

    logger.info('Running remap.')
    with open(os.path.join(sample_scratch_path, 'prelim.csv'), 'rU') as prelim_csv, \
            open(os.path.join(sample_scratch_path, 'remap.csv'), 'wb') as remap_csv, \
            open(os.path.join(sample_out_path, 'remap_counts.csv'), 'wb') as counts_csv, \
            open(os.path.join(sample_out_path, 'remap_conseq.csv'), 'wb') as conseq_csv, \
            open(os.path.join(sample_out_path, 'unmapped1.fastq'), 'w') as unmapped1, \
            open(os.path.join(sample_out_path, 'unmapped2.fastq'), 'w') as unmapped2:

        remap(censored_path1,
              censored_path2,
              prelim_csv,
              remap_csv,
              counts_csv,
              conseq_csv,
              unmapped1,
              unmapped2,
              scratch_path)

    logger.info("Running sam2aln.")
    with open(os.path.join(sample_scratch_path, 'remap.csv'), 'rU') as remap_csv, \
            open(os.path.join(sample_scratch_path, 'aligned.csv'), 'wb') as aligned_csv, \
            open(os.path.join(sample_out_path, 'conseq_ins.csv'), 'wb') as insert_csv, \
            open(os.path.join(sample_out_path, 'failed_read.csv'), 'wb') as failed_csv:

        sam2aln(remap_csv, aligned_csv, insert_csv, failed_csv)

    logger.info("Running aln2counts.")
    with open(os.path.join(sample_scratch_path, 'aligned.csv'), 'rU') as aligned_csv, \
            open(os.path.join(sample_out_path, 'nuc.csv'), 'wb') as nuc_csv, \
            open(os.path.join(sample_out_path, 'amino.csv'), 'wb') as amino_csv, \
            open(os.path.join(sample_out_path, 'coord_ins.csv'), 'wb') as coord_ins_csv, \
            open(os.path.join(sample_out_path, 'conseq.csv'), 'wb') as conseq_csv, \
            open(os.path.join(sample_out_path, 'failed_align.csv'), 'wb') as failed_align_csv, \
            open(os.path.join(sample_out_path, 'nuc_variants.csv'), 'wb') as nuc_variants_csv:

        aln2counts(aligned_csv,
                   nuc_csv,
                   amino_csv,
                   coord_ins_csv,
                   conseq_csv,
                   failed_align_csv,
                   nuc_variants_csv)

    logger.info("Running sam_g2p.")
    with open(os.path.join(sample_scratch_path, 'remap.csv'), 'rU') as remap_csv, \
            open(os.path.join(sample_out_path, 'nuc.csv'), 'rU') as nuc_csv, \
            open(os.path.join(sample_out_path, 'g2p.csv'), 'wb') as g2p_csv, \
            open(os.path.join(sample_out_path, 'g2p_summary.csv'), 'wb') as g2p_summary_csv:

        sam_g2p(pssm=pssm,
                remap_csv=remap_csv,
                nuc_csv=nuc_csv,
                g2p_csv=g2p_csv,
                g2p_summary_csv=g2p_summary_csv)


def parse_phix(args, json):
    read_lengths = [json.read_length1,
                    json.index_length1,
                    json.index_length2,
                    json.read_length2]

    phix_path = os.path.join(args.data_path,
                             'input',
                             'runs',
                             json.run_id,
                             'InterOp',
                             'ErrorMetricsOut.bin')
    quality_path = os.path.join(args.data_path,
                                'scratch',
                                'quality.csv')
    bad_cycles_path = os.path.join(args.data_path,
                                   'scratch',
                                   'bad_cycles.csv')
    with open(phix_path, 'rb') as phix, open(quality_path, 'w') as quality:
        records = phix_parser.read_phix(phix)
        phix_parser.write_phix_csv(quality, records, read_lengths)
    with open(quality_path, 'rU') as quality, open(bad_cycles_path, 'w') as bad_cycles:
        report_bad_cycles(quality, bad_cycles)


def makedirs(path):
    try:
        os.makedirs(path)
    except OSError as exc:
        if exc.errno != errno.EEXIST or not os.path.isdir(path):
            raise


def main():
    logger.info('Starting.')
    args = parse_args()
    json_path = os.path.join(args.data_path, 'input', 'AppSession.json')
    with open(json_path, 'rU') as json_file:
        json = parse_json(json_file)
    pssm = Pssm()

    scratch_path = os.path.join(args.data_path, 'scratch')
    makedirs(scratch_path)
    for filename in os.listdir(scratch_path):
        filepath = os.path.join(scratch_path, filename)
        if os.path.isdir(filepath):
            shutil.rmtree(filepath)
        else:
            os.remove(filepath)

    logger.info('Processing error rates.')
    if json.run_id is not None:
        parse_phix(args, json)

    for sample_info in json.samples:
        process_sample(sample_info, json.project_id, args.data_path, pssm)

    listing_path = os.path.join(args.data_path,
                                'output',
                                'appresults',
                                json.project_id,
                                json.samples[0]['Name'],
                                'listing.txt')
    with open(listing_path, 'w') as listing:
        listing.write(subprocess.check_output(['ls', '-R', args.data_path]))
    logger.info('Done.')

if __name__ == '__main__':
    main()
elif __name__ == '__live_coding__':
    from cStringIO import StringIO
    json_file = StringIO("""\
{"Properties": {"Items": [{"Name": "Input.app-session-name",
                           "Content": "MiCall 04/05/2016 3:14:23"},
                          {"Name": "Input.sample-ids",
                           "Items": [{"Id": "11111",
                                      "Name": "Example-SampleB"},
                                     {"Id": "22222",
                                      "Name": "Example-SampleA"}]},
                          {"Name": "Input.project-id",
                           "Content": {"Id": "33333",
                                       "Name": "Example-Project"}},
                          {"Name": "Input.run-id",
                           "Content": {"Id": "44444",
                                       "Name": "160115_Project",
                                       "SequencingStats": {
                                            "NumCyclesIndex1": 6,
                                            "NumCyclesIndex2": 0,
                                            "NumCyclesRead1": 251,
                                            "NumCyclesRead2": 251}}}]}}
""")
    parse_json(json_file)
