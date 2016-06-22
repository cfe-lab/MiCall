from argparse import ArgumentParser
import csv
import errno
import fnmatch
from glob import glob
import json
import logging
import multiprocessing
from operator import itemgetter
import os
import shutil
import socket
import subprocess
from xml.etree import ElementTree

from micall.core.aln2counts import aln2counts
from micall.core.censor_fastq import censor
from micall.core.filter_quality import report_bad_cycles
from micall.core.remap import remap
from micall.core.prelim_map import prelim_map
from micall.core.sam2aln import sam2aln
from micall.monitor import error_metrics_parser, quality_metrics_parser
from micall.g2p.sam_g2p import sam_g2p
from micall.g2p.pssm_lib import Pssm
from micall.monitor.tile_metrics_parser import summarize_tiles
from micall.utils.coverage_plots import coverage_plot

logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s[%(levelname)s]%(name)s.%(funcName)s(): %(message)s')
logger = logging.getLogger('micall')

MIN_G2P_COUNT = 3


def parse_args():
    parser = ArgumentParser(description='Map FASTQ files to references.')
    parser.add_argument('data_path',
                        default='/data',
                        nargs='?',
                        help='data folder filled in by BaseSpace')
    parser.add_argument('--link_run',
                        '-l',
                        help='Run folder to link into the data folder')
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
    args.href_app_session = raw_args['Href']
    run = arg_map.get('Input.run-id')
    if run is None:
        args.run_id = None
    else:
        run_content = run['Content']
        args.run_id = run_content['Id']
        args.read_length1 = run_content['SequencingStats']['NumCyclesRead1']
        args.read_length2 = run_content['SequencingStats']['NumCyclesRead2']
        args.index_length1 = run_content['SequencingStats']['NumCyclesIndex1']
        args.index_length2 = run_content['SequencingStats']['NumCyclesIndex2']
    args.samples = sorted(arg_map['Input.sample-ids']['Items'],
                          key=itemgetter('Name'))
    args.project_id = arg_map['Input.project-id']['Content']['Id']

    return args


def link_json(run_path, data_path):
    """ Load the data from a run folder into the BaseSpace layout. """
    class Args(object):
        pass
    args = Args()

    shutil.rmtree(data_path, ignore_errors=True)
    makedirs(data_path)

    args.run_id = os.path.basename(run_path)
    runs_path = os.path.join(data_path, 'input', 'runs')
    makedirs(runs_path)
    new_run_path = os.path.join(runs_path, args.run_id)
    os.symlink(run_path, new_run_path)
    run_info_path = os.path.join(new_run_path, 'RunInfo.xml')
    run_info = ElementTree.parse(run_info_path).getroot()
    read1 = run_info.find('.//Read[@Number="1"][@IsIndexedRead="N"]')
    args.read_length1 = int(read1.attrib['NumCycles'])
    read2 = run_info.find('.//Read[@IsIndexedRead="N"][last()]')
    args.read_length2 = int(read2.attrib['NumCycles'])
    index1 = run_info.find('.//Read[@Number="2"][@IsIndexedRead="Y"]')
    args.index_length1 = int(index1.attrib['NumCycles'])
    index2 = run_info.find('.//Read[@Number="3"][@IsIndexedRead="Y"]')
    if index2 is None:
        args.index_length2 = 0
    else:
        args.index_length2 = int(index2.attrib['NumCycles'])
    args.project_id = args.href_app_session = '1'

    args.samples = []
    samples_path = os.path.join(data_path, 'input', 'samples')
    fastq_files = glob(os.path.join(run_path,
                                    'Data',
                                    'Intensities',
                                    'BaseCalls',
                                    '*_R1_*'))
    fastq_files.sort()
    for i, fastq_file in enumerate(fastq_files, 1):
        sample_file = os.path.basename(fastq_file)
        if not sample_file.startswith('Undetermined'):
            sample_id = str(i)
            sample_path = os.path.join(samples_path, sample_id)
            makedirs(sample_path)
            os.symlink(fastq_file, os.path.join(sample_path, sample_file))
            fastq_file = fastq_file.replace('_R1_', '_R2_')
            sample_file = os.path.basename(fastq_file)
            os.symlink(fastq_file, os.path.join(sample_path, sample_file))
            sample_name = '_'.join(sample_file.split('_')[:2])
            args.samples.append(dict(Id=sample_id,
                                     Href="v1pre3/samples/" + sample_id,
                                     Name=sample_name))
    return args


def censor_sample(filename, bad_cycles_path, censored_name, read_summary_name):
    if not os.path.exists(bad_cycles_path):
        bad_cycles = []
    else:
        with open(bad_cycles_path, 'rU') as bad_cycles:
            bad_cycles = list(csv.DictReader(bad_cycles))
    with open(filename, 'rb') as fastq_src,\
            open(censored_name, 'w') as fastq_dest,\
            open(read_summary_name, 'w') as read_summary:
        censor(fastq_src, bad_cycles, fastq_dest, summary_file=read_summary)


def build_app_result_path(data_path,
                          run_info,
                          sample_info,
                          suffix=None):
    dir_name = sample_info['Name']
    if suffix is not None:
        dir_name += suffix
    sample_out_path = os.path.join(data_path,
                                   'output',
                                   'appresults',
                                   run_info.project_id,
                                   dir_name)
    return sample_out_path


def create_app_result(data_path,
                      run_info,
                      sample_info,
                      description='',
                      suffix=None):
    sample_out_path = build_app_result_path(data_path,
                                            run_info,
                                            sample_info,
                                            suffix)
    makedirs(sample_out_path)
    metadata = dict(Name=os.path.basename(sample_out_path),
                    Description=description,
                    HrefAppSession=run_info.href_app_session,
                    Properties=[dict(Type='sample',
                                     Name='Input.Samples',
                                     Content=sample_info['Href'])])
    with open(os.path.join(sample_out_path, '_metadata.json'), 'w') as json_file:
        json.dump(metadata, json_file, indent=4)
    return sample_out_path


def process_sample(sample_index, run_info, data_path, pssm):
    """ Process a single sample.

    :param sample_index: which sample to process from the session JSON
    :param run_info: run parameters loaded from the session JSON
    :param str data_path: the root folder for all BaseSpace data
    :param pssm: the pssm library for running G2P analysis
    """
    scratch_path = os.path.join(data_path, 'scratch')
    sample_info = run_info.samples[sample_index]
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
    logger.info('Processing sample %s (%d of %d): %s (%s).',
                sample_id,
                sample_index+1,
                len(run_info.samples),
                sample_name,
                sample_path)

    sample_out_path = create_app_result(data_path,
                                        run_info,
                                        sample_info,
                                        description='Mapping results',
                                        suffix='_QC')

    sample_scratch_path = os.path.join(scratch_path, sample_name)
    makedirs(sample_scratch_path)

    censored_path1 = os.path.join(sample_scratch_path, 'censored1.fastq')
    read_summary_path1 = os.path.join(sample_scratch_path, 'read1_summary.csv')
    censor_sample(sample_path,
                  os.path.join(scratch_path, 'bad_cycles.csv'),
                  censored_path1,
                  read_summary_path1)
    censored_path2 = os.path.join(sample_scratch_path, 'censored2.fastq')
    read_summary_path2 = os.path.join(sample_scratch_path, 'read2_summary.csv')
    censor_sample(sample_path2,
                  os.path.join(scratch_path, 'bad_cycles.csv'),
                  censored_path2,
                  read_summary_path2)

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
            open(os.path.join(sample_out_path, 'nuc_variants.csv'), 'wb') as nuc_variants_csv, \
            open(os.path.join(sample_scratch_path, 'coverage_summary.csv'), 'wb') as coverage_summary_csv:

        aln2counts(aligned_csv,
                   nuc_csv,
                   amino_csv,
                   coord_ins_csv,
                   conseq_csv,
                   failed_align_csv,
                   nuc_variants_csv,
                   coverage_summary_csv=coverage_summary_csv)

    logger.info("Running coverage_plots.")
    coverage_path = os.path.join(sample_out_path, 'coverage')
    with open(os.path.join(sample_out_path, 'amino.csv'), 'rU') as amino_csv, \
            open(os.path.join(sample_out_path, 'coverage_scores.csv'), 'w') as coverage_scores_csv:
        coverage_plot(amino_csv, coverage_scores_csv, path_prefix=coverage_path)

    with open(os.path.join(sample_out_path, 'coverage_scores.csv'), 'rU') as coverage_scores_csv:
        reader = csv.DictReader(coverage_scores_csv)
        is_v3loop_good = False
        for row in reader:
            if row['region'] == 'V3LOOP':
                is_v3loop_good = row['on.score'] == '4'
                break

    if is_v3loop_good:
        logger.info("Running sam_g2p.")
        g2p_path = create_app_result(data_path,
                                     run_info,
                                     sample_info,
                                     description='Geno To Pheno results',
                                     suffix='_G2P')
        with open(os.path.join(sample_scratch_path, 'remap.csv'), 'rU') as remap_csv, \
                open(os.path.join(sample_out_path, 'nuc.csv'), 'rU') as nuc_csv, \
                open(os.path.join(g2p_path, 'g2p.csv'), 'wb') as g2p_csv, \
                open(os.path.join(g2p_path, 'g2p_summary.csv'), 'wb') as g2p_summary_csv:

            sam_g2p(pssm=pssm,
                    remap_csv=remap_csv,
                    nuc_csv=nuc_csv,
                    g2p_csv=g2p_csv,
                    g2p_summary_csv=g2p_summary_csv,
                    min_count=MIN_G2P_COUNT)


def summarize_run(args, json):
    """ Summarize the run data from the InterOp folder.

    Writes some summary files.
    :return: a dictionary with summary values.
    """
    read_lengths = [json.read_length1,
                    json.index_length1,
                    json.index_length2,
                    json.read_length2]
    summary = {}

    interop_path = os.path.join(args.data_path,
                                'input',
                                'runs',
                                json.run_id,
                                'InterOp')
    phix_path = os.path.join(interop_path, 'ErrorMetricsOut.bin')
    quality_path = os.path.join(args.data_path, 'scratch', 'quality.csv')
    bad_cycles_path = os.path.join(args.data_path, 'scratch', 'bad_cycles.csv')
    summary_path = build_app_result_path(args.data_path,
                                         json,
                                         json.samples[0],
                                         suffix='_QC')
    makedirs(summary_path)
    bad_tiles_path = os.path.join(summary_path, 'bad_tiles.csv')
    with open(phix_path, 'rb') as phix, open(quality_path, 'w') as quality:
        records = error_metrics_parser.read_errors(phix)
        error_metrics_parser.write_phix_csv(quality,
                                            records,
                                            read_lengths,
                                            summary)
    with open(quality_path, 'rU') as quality, \
            open(bad_cycles_path, 'w') as bad_cycles, \
            open(bad_tiles_path, 'w') as bad_tiles:
        report_bad_cycles(quality, bad_cycles, bad_tiles)

    quality_metrics_path = os.path.join(interop_path, 'QMetricsOut.bin')
    quality_metrics_parser.summarize_quality(quality_metrics_path,
                                             summary,
                                             read_lengths)

    tile_metrics_path = os.path.join(interop_path, 'TileMetricsOut.bin')
    summarize_tiles(tile_metrics_path, summary)
    return summary


def summarize_samples(args, json, run_summary):
    summary_path = build_app_result_path(args.data_path,
                                         json,
                                         json.samples[0],
                                         suffix='_QC')

    score_sum = 0.0
    base_count = 0
    coverage_sum = 0.0
    coverage_count = 0
    for sample in json.samples:
        sample_scratch_path = os.path.join(args.data_path,
                                           'scratch',
                                           sample['Name'])
        for filename in ('read1_summary.csv', 'read2_summary.csv'):
            read_summary_path = os.path.join(sample_scratch_path, filename)
            with open(read_summary_path, 'rU') as read_summary:
                reader = csv.DictReader(read_summary)
                row = reader.next()
                sample_base_count = int(row['base_count'])
                if sample_base_count:
                    score_sum += float(row['avg_quality']) * sample_base_count
                    base_count += sample_base_count
        coverage_summary_path = os.path.join(sample_scratch_path,
                                             'coverage_summary.csv')
        with open(coverage_summary_path, 'rU') as coverage_summary:
            reader = csv.DictReader(coverage_summary)
            for row in reader:
                region_width = int(row['region_width'])
                coverage_sum += float(row['avg_coverage']) * region_width
                coverage_count += region_width

    if base_count > 0:
        run_summary['avg_quality'] = score_sum / base_count
    if coverage_count > 0:
        run_summary['avg_coverage'] = coverage_sum / coverage_count

    run_quality_path = os.path.join(summary_path, 'run_quality.csv')
    with open(run_quality_path, 'w') as run_quality:
        writer = csv.DictWriter(run_quality,
                                ['q30_fwd',
                                 'q30_rev',
                                 'cluster_density',
                                 'pass_rate',
                                 'error_rate_fwd',
                                 'error_rate_rev',
                                 'avg_quality',
                                 'avg_coverage'],
                                lineterminator=os.linesep)
        writer.writeheader()
        writer.writerow(run_summary)

    listing_path = os.path.join(summary_path, 'listing.txt')
    with open(listing_path, 'w') as listing:
        listing.write(subprocess.check_output(['ls', '-R', args.data_path]))


def makedirs(path):
    try:
        os.makedirs(path)
    except OSError as exc:
        if exc.errno != errno.EEXIST or not os.path.isdir(path):
            raise


def main():
    logger.info("Starting on %s with %d CPU's.",
                socket.gethostname(),
                multiprocessing.cpu_count())
    args = parse_args()
    if args.link_run is not None:
        json = link_json(args.link_run, args.data_path)
    else:
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

    if json.run_id is not None:
        logger.info('Summarizing run.')
        run_summary = summarize_run(args, json)

    for sample_index in range(len(json.samples)):
        process_sample(sample_index, json, args.data_path, pssm)

    if json.run_id is not None:
        summarize_samples(args, json, run_summary)
    logger.info('Done.')

if __name__ == '__main__':
    main()
elif __name__ == '__live_coding__':
    from cStringIO import StringIO
    json_file = StringIO("""\
{"Href": "v1pre3/appsessions/1234",
 "Properties": {"Items": [{"Name": "Input.app-session-name",
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
