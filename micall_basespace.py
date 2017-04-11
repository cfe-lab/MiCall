from argparse import ArgumentParser
import csv
import errno
import fnmatch
import functools
from glob import glob
import json
import logging
import multiprocessing
from multiprocessing.pool import Pool
from operator import itemgetter
import os
import shutil
import socket
import requests
from xml.etree import ElementTree

from micall.core.aln2counts import aln2counts
from micall.core.trim_fastqs import trim
from micall.core.filter_quality import report_bad_cycles
from micall.core.remap import remap
from micall.core.prelim_map import prelim_map
from micall.core.sam2aln import sam2aln
from micall.monitor import error_metrics_parser, quality_metrics_parser
from micall.g2p.fastq_g2p import fastq_g2p, DEFAULT_MIN_COUNT
from micall.g2p.pssm_lib import Pssm
from micall.monitor.tile_metrics_parser import summarize_tiles
from micall.core.coverage_plots import coverage_plot

EXCLUDED_SEEDS = ['HLA-B-seed']  # Not ready yet.
EXCLUDED_PROJECTS = ['HCV-NS5a',
                     'HIV-GP41',
                     'HIV-PRI',
                     'INT',
                     'MidHCV',
                     'MiniHCV',
                     'MiniRT',
                     'PR-RT',
                     'RT',
                     'V3LOOP',
                     'wg1HCV']  # Avoid useless duplicates for BaseSpace version.
DOWNLOAD_BATCH_SIZE = 1000
logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s[%(levelname)s]%(name)s.%(funcName)s(): %(message)s')
logger = logging.getLogger('micall')


class BSrequest:
    """A class that wraps some of the Basespace RESTful API"""

    def __init__(self,
                 base_url="https://api.basespace.illumina.com/v1pre3",
                 access_token=None):
        self._burl = base_url
        if access_token is None:
            envtoken = os.getenv("AccessToken", None)
            if envtoken is None:
                logger.error("Cannot find environment variable called 'AccessToken'")
            self._access_token = envtoken
        else:
            self._access_token = access_token

    @staticmethod
    def _responsestatus(jobj):
        """ Return the response code of jobj (expected to be a response from a API query).
        This consists of an Basespace errorcode and a message (two strings).

        The dictionary ResponseStatus can be empty if the request was successful.
        Return empty strings in this case.

        NOTE: There is an inconsistency in the API response codes:
        If the API is given an unrecognised URL path, it returns an
        empty ErrorCode, and a nonempty Message. For consistency, it should
        return a nonempty ErrorCode in this case.
        This would confuse routines calling this one, because they check for errors
        like this:
        err_code, err_msg = bs.responsestatus(jobj)
        if err_code:
           print "Houston, we have a problem"
        else:
           print "Everything is hunky dory"

        """
        resp = jobj["ResponseStatus"]
        err_code = resp.get("ErrorCode", "")
        err_msg = resp.get("Message", "")
        # see docstring for the reasoning behind this logic
        if err_msg and not err_code:
            err_code = "API ERROR"
        return err_code, err_msg

    def _raw_get_file(self, locstr, paramdct=None):
        """ Perform a RESTful basespace API request using
        locstr: the API request without the base URL, and
        optional parameters in the paramdct.
        These must be of the form defined the ine basespace API docs.

        Return the requests object. Note that if this expected to contain a json
        object, the this must be extracted by obj.json().
        """
        if self._access_token is None:
            logger.error("skipping _raw_get_file '%s': no access token" % locstr)
            raise RuntimeError("no access_token for basespace API")
        if paramdct is None:
            pdct = {}
        else:
            pdct = paramdct.copy()
        assert "access_token" not in pdct, " access_token may not be in paramdct"
        pdct["access_token"] = self._access_token
        reqstr = "/".join([self._burl, locstr])
        logger.info("API: '%s'" % reqstr)
        return requests.get(reqstr, params=pdct)

    def _runssamples(self, runid, limit, offset):
        """ Return the json result of the runs/id/samples API request:
        https://developer.basespace.illumina.com/docs/content/documentation/rest-api/api-reference#Samples
        This returns a json object describing samples  belonging to a specified run, or
        an object with an error message.

        Limit: limit the number of responses to Limit number of items.
        Offset: the offset of the first file to read in the whole collection
        """
        jobj = self._raw_get_file("/".join(["runs", runid, "samples"]),
                                  paramdct={"Limit": limit, "Offset": offset}).json()
        err_code, _err_msg = self._responsestatus(jobj)
        if err_code:
            raise RuntimeError("runsamples API error")
        return jobj

    def _get_all_sample_ids_from_run_id(self, runid):
        """Retrieve a set of the sampleids of all the samples in the specified run.
        Raise a Runtimeerror iff an API error occurs.

        NOTE: This routine will download sample_id in batches of DOWNLOAD_BATCH_SIZE.
        """
        numgot = 0
        retset = set()
        needmore = True
        while needmore:
            jobj = self._runssamples(runid, DOWNLOAD_BATCH_SIZE, numgot)
            response = jobj["Response"]
            retset |= set([s["Id"] for s in response["Items"]])
            numneed = response["TotalCount"]
            numgot += response["DisplayedCount"]
            needmore = (numgot < numneed)
        return retset

    def check_run_sample_ids(self, run_id_lst, sample_id_lst):
        """"Make sure that all samples in the sample_id_lst
        belong to runs specified in the run_id_lst.

        Return a set of those sample_ids in sample_id_lst that are found in the
        runs specified in run_id_lst.
        """
        # a: determine the set of all sampleids in all runs
        sample_set = set().union(*[self._get_all_sample_ids_from_run_id(runid) for runid in run_id_lst])
        return set(sample_id_lst) & sample_set


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


class Args(object):
    pass


def parse_json(json_file):
    """ Load JSON from an open file, and pull out the arguments for this run.

    :param json_file: an open file that contains JSON in the BaseSpace
    AppSession format.
    :return: an object with an attribute for each argument
    """
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
    args = Args()
    args.project_id = args.href_app_session = '1'

    shutil.rmtree(data_path, ignore_errors=True)
    makedirs(data_path)

    results_path = os.path.join(run_path, 'Results', 'basespace')
    makedirs(results_path)
    appresults_path = os.path.join(data_path, 'output', 'appresults')
    makedirs(appresults_path)
    os.symlink(results_path, os.path.join(appresults_path, args.project_id))

    args.run_id = os.path.basename(run_path)
    input_runs_path = os.path.join(data_path, 'input', 'runs')
    makedirs(input_runs_path)
    new_run_path = os.path.join(input_runs_path, args.run_id)
    os.symlink(run_path, new_run_path)
    run_info_path = os.path.join(new_run_path, 'RunInfo.xml')
    if not os.path.exists(run_info_path):
        args.run_id = None
    else:
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

    args.samples = []
    samples_path = os.path.join(data_path, 'input', 'samples')
    fastq_files = (glob(os.path.join(run_path,
                                     'Data',
                                     'Intensities',
                                     'BaseCalls',
                                     '*_R1_*')) or
                   glob(os.path.join(run_path,
                                     '*_R1_*')))
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


def build_app_result_path(data_path,
                          run_info,
                          sample_info=None,
                          suffix=None):
    dir_name = sample_info['Name'] if sample_info else ''
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
                      sample_info=None,
                      description='',
                      suffix=None):
    sample_out_path = build_app_result_path(data_path,
                                            run_info,
                                            sample_info,
                                            suffix)
    makedirs(sample_out_path)
    if sample_info is None:
        sample_properties = dict(Type='sample[]',
                                 Name='Input.Samples',
                                 Items=list(map(itemgetter('Href'),
                                                run_info.samples)))
    else:
        sample_properties = dict(Type='sample',
                                 Name='Input.Samples',
                                 Content=sample_info['Href'])
    metadata = dict(Name=os.path.basename(sample_out_path),
                    Description=description,
                    HrefAppSession=run_info.href_app_session,
                    Properties=sample_properties)
    with open(os.path.join(sample_out_path, '_metadata.json'), 'w') as json_file:
        json.dump(metadata, json_file, indent=4)
    return sample_out_path


def try_sample(sample_index, run_info, args, pssm):
    """ Try processing a single sample.

    Tracebacks and some errors can't be pickled across process boundaries, so
    log detailed error before raising a RuntimeError.
    """
    try:
        process_sample(sample_index, run_info, args, pssm)
    except:
        message = 'Failed to process sample {}.'.format(sample_index+1)
        logger.error(message, exc_info=True)
        raise RuntimeError(message)


def process_sample(sample_index, run_info, args, pssm):
    """ Process a single sample.

    :param sample_index: which sample to process from the session JSON
    :param run_info: run parameters loaded from the session JSON
    :param args: the command-line arguments
    :param pssm: the pssm library for running G2P analysis
    """
    scratch_path = os.path.join(args.data_path, 'scratch')
    sample_info = run_info.samples[sample_index]
    sample_id = sample_info['Id']
    sample_name = sample_info['Name']
    sample_dir = os.path.join(args.data_path,
                              'input',
                              'samples',
                              sample_id,
                              'Data',
                              'Intensities',
                              'BaseCalls')
    if not os.path.exists(sample_dir):
        sample_dir = os.path.join(args.data_path,
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

    sample_qc_path = os.path.join(args.qc_path, sample_name)
    makedirs(sample_qc_path)
    sample_scratch_path = os.path.join(scratch_path, sample_name)
    makedirs(sample_scratch_path)

    bad_cycles_path = os.path.join(scratch_path, 'bad_cycles.csv')
    trimmed_path1 = os.path.join(sample_scratch_path, 'trimmed1.fastq')
    read_summary_path = os.path.join(sample_scratch_path, 'read_summary.csv')
    trimmed_path2 = os.path.join(sample_scratch_path, 'trimmed2.fastq')
    with open(read_summary_path, 'w') as read_summary:
        trim((sample_path, sample_path2),
             bad_cycles_path,
             (trimmed_path1, trimmed_path2),
             summary_file=read_summary,
             use_gzip=sample_path.endswith('.gz'))

    logger.info('Running prelim_map (%d of %d).', sample_index+1, len(run_info.samples))
    with open(os.path.join(sample_scratch_path, 'prelim.csv'), 'w') as prelim_csv:
        prelim_map(trimmed_path1,
                   trimmed_path2,
                   prelim_csv,
                   work_path=sample_scratch_path,
                   excluded_seeds=EXCLUDED_SEEDS)

    logger.info('Running remap (%d of %d).', sample_index+1, len(run_info.samples))
    debug_file_prefix = None  # os.path.join(sample_scratch_path, 'debug')
    with open(os.path.join(sample_scratch_path, 'prelim.csv'), 'rU') as prelim_csv, \
            open(os.path.join(sample_scratch_path, 'remap.csv'), 'w') as remap_csv, \
            open(os.path.join(sample_scratch_path, 'remap_counts.csv'), 'w') as counts_csv, \
            open(os.path.join(sample_scratch_path, 'remap_conseq.csv'), 'w') as conseq_csv, \
            open(os.path.join(sample_qc_path, 'unmapped1.fastq'), 'w') as unmapped1, \
            open(os.path.join(sample_qc_path, 'unmapped2.fastq'), 'w') as unmapped2:

        remap(trimmed_path1,
              trimmed_path2,
              prelim_csv,
              remap_csv,
              counts_csv,
              conseq_csv,
              unmapped1,
              unmapped2,
              sample_scratch_path,
              debug_file_prefix=debug_file_prefix,
              nthreads=1)

    logger.info('Running sam2aln (%d of %d).', sample_index+1, len(run_info.samples))
    with open(os.path.join(sample_scratch_path, 'remap.csv'), 'rU') as remap_csv, \
            open(os.path.join(sample_scratch_path, 'aligned.csv'), 'w') as aligned_csv, \
            open(os.path.join(sample_scratch_path, 'conseq_ins.csv'), 'w') as conseq_ins_csv, \
            open(os.path.join(sample_scratch_path, 'failed_read.csv'), 'w') as failed_csv, \
            open(os.path.join(sample_scratch_path, 'clipping.csv'), 'w') as clipping_csv:

        sam2aln(remap_csv,
                aligned_csv,
                conseq_ins_csv,
                failed_csv,
                clipping_csv=clipping_csv)

    logger.info('Running aln2counts (%d of %d).', sample_index+1, len(run_info.samples))
    with open(os.path.join(sample_scratch_path, 'aligned.csv'), 'rU') as aligned_csv, \
            open(os.path.join(sample_scratch_path, 'clipping.csv'), 'rU') as clipping_csv, \
            open(os.path.join(sample_scratch_path, 'conseq_ins.csv'), 'rU') as conseq_ins_csv, \
            open(os.path.join(sample_scratch_path, 'nuc.csv'), 'w') as nuc_csv, \
            open(os.path.join(sample_scratch_path, 'amino.csv'), 'w') as amino_csv, \
            open(os.path.join(sample_scratch_path, 'coord_ins.csv'), 'w') as coord_ins_csv, \
            open(os.path.join(sample_scratch_path, 'conseq.csv'), 'w') as conseq_csv, \
            open(os.path.join(sample_scratch_path, 'failed_align.csv'), 'w') as failed_align_csv, \
            open(os.path.join(sample_scratch_path, 'nuc_variants.csv'), 'w') as nuc_variants_csv, \
            open(os.path.join(sample_scratch_path, 'coverage_summary.csv'), 'w') as coverage_summary_csv:

        aln2counts(aligned_csv,
                   nuc_csv,
                   amino_csv,
                   coord_ins_csv,
                   conseq_csv,
                   failed_align_csv,
                   nuc_variants_csv,
                   coverage_summary_csv=coverage_summary_csv,
                   clipping_csv=clipping_csv,
                   conseq_ins_csv=conseq_ins_csv)

    logger.info('Running coverage_plots (%d of %d).', sample_index+1, len(run_info.samples))
    coverage_maps_path = os.path.join(args.qc_path, 'coverage_maps')
    makedirs(coverage_maps_path)
    with open(os.path.join(sample_scratch_path, 'amino.csv'), 'rU') as amino_csv, \
            open(os.path.join(sample_scratch_path, 'coverage_scores.csv'), 'w') as coverage_scores_csv:
        coverage_plot(amino_csv,
                      coverage_scores_csv,
                      coverage_maps_path=coverage_maps_path,
                      coverage_maps_prefix=sample_name,
                      excluded_projects=EXCLUDED_PROJECTS)

    logger.info('Running fastq_g2p (%d of %d).', sample_index+1, len(run_info.samples))
    with open(os.path.join(sample_scratch_path, 'trimmed1.fastq'), 'rU') as fastq1, \
            open(os.path.join(sample_scratch_path, 'trimmed2.fastq'), 'rU') as fastq2, \
            open(os.path.join(sample_scratch_path, 'g2p.csv'), 'w') as g2p_csv, \
            open(os.path.join(sample_scratch_path, 'g2p_summary.csv'), 'w') as g2p_summary_csv:

        fastq_g2p(pssm=pssm,
                  fastq1=fastq1,
                  fastq2=fastq2,
                  g2p_csv=g2p_csv,
                  g2p_summary_csv=g2p_summary_csv,
                  min_count=DEFAULT_MIN_COUNT)
    logger.info('Finished sample (%d of %d).', sample_index+1, len(run_info.samples))


def summarize_run(args, run_json):
    """ Summarize the run data from the InterOp folder.

    Writes some summary files.
    :return: a dictionary with summary values.
    """
    read_lengths = [run_json.read_length1,
                    run_json.index_length1,
                    run_json.index_length2,
                    run_json.read_length2]
    summary = {}

    has_error_metrics = run_json.has_runinfo

    if has_error_metrics:
        interop_path = os.path.join(args.data_path,
                                    'input',
                                    'runs',
                                    run_json.run_id,
                                    'InterOp')
        phix_path = os.path.join(interop_path, 'ErrorMetricsOut.bin')
        quality_path = os.path.join(args.data_path, 'scratch', 'quality.csv')
        bad_cycles_path = os.path.join(args.data_path, 'scratch', 'bad_cycles.csv')
        bad_tiles_path = os.path.join(args.qc_path, 'bad_tiles.csv')
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


def summarize_samples(args, run_json, run_summary):
    score_sum = 0.0
    base_count = 0
    coverage_sum = 0.0
    coverage_count = 0
    for sample in run_json.samples:
        sample_scratch_path = os.path.join(args.data_path,
                                           'scratch',
                                           sample['Name'])
        read_summary_path = os.path.join(sample_scratch_path, 'read_summary.csv')
        with open(read_summary_path, 'rU') as read_summary:
            reader = csv.DictReader(read_summary)
            for row in reader:
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

    run_quality_path = os.path.join(args.qc_path, 'run_quality.csv')
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


def collate_samples(args, run_info):
    scratch_path = os.path.join(args.data_path, 'scratch')
    filenames = ['remap_counts.csv',
                 'remap_conseq.csv',
                 'conseq_ins.csv',
                 'failed_read.csv',
                 'nuc.csv',
                 'amino.csv',
                 'coord_ins.csv',
                 'conseq.csv',
                 'failed_align.csv',
                 'nuc_variants.csv',
                 'coverage_scores.csv',
                 'g2p.csv',
                 'g2p_summary.csv']
    for filename in filenames:
        out_path = args.qc_path if filename == 'coverage_scores.csv' else args.g2p_path
        with open(os.path.join(out_path, filename), 'w') as fout:
            writer = csv.writer(fout, lineterminator=os.linesep)
            is_header_written = False
            for sample_info in run_info.samples:
                sample_name = sample_info['Name']
                sample_scratch_path = os.path.join(scratch_path, sample_name)
                srcfile = os.path.join(sample_scratch_path, filename)
                try:
                    with open(srcfile, 'rU') as fin:
                        reader = csv.reader(fin)
                        for i, row in enumerate(reader):
                            if i == 0:
                                if not is_header_written:
                                    row.insert(0, 'sample')
                                    writer.writerow(row)
                                    is_header_written = True
                            else:
                                row.insert(0, sample_name)
                                writer.writerow(row)
                except IOError as ex:
                    if ex.errno != errno.ENOENT:
                        raise


def makedirs(path):
    try:
        os.makedirs(path)
    except OSError as exc:
        if exc.errno != errno.EEXIST or not os.path.isdir(path):
            raise


# noinspection PyTypeChecker,PyUnresolvedReferences
def main():
    logger.info("Starting on %s with %d CPU's.",
                socket.gethostname(),
                multiprocessing.cpu_count())
    args = parse_args()
    if args.link_run is not None:
        run_json = link_json(args.link_run, args.data_path)
        run_json.has_runinfo = True
    else:
        json_path = os.path.join(args.data_path, 'input', 'AppSession.json')
        try:
            with open(json_path, 'rU') as json_file:
                run_json = parse_json(json_file)
        except:
            if os.path.exists(json_path):
                # copy the input file to the output dir for postmortem analysis
                logger.error("Error occurred while parsing '%s'" % json_path)
                with open(json_path, 'rU') as json_file:
                    file_cont = json_file.read()
                out_path = os.path.join(args.data_path, 'logs', 'AppSession.json')
                with open(out_path, 'w') as json_file:
                    json_file.write(file_cont)
            else:
                logger.error("Error: no such file as '%s'" % json_path)
            raise
        # Do we have run_ids for all sample_ids ?
        if run_json.run_id is None:
            run_json.has_runinfo = False
        else:
            bs = BSrequest()
            sample_id_set = bs.check_run_sample_ids([run_json.run_id],
                                                    [s["Id"] for s in run_json.samples])
            run_json.has_runinfo = (len(sample_id_set) == len(run_json.samples))
        logger.info("setting json.has_run_info to %s" % run_json.has_runinfo)
    pssm = Pssm()

    scratch_path = os.path.join(args.data_path, 'scratch')
    makedirs(scratch_path)
    for filename in os.listdir(scratch_path):
        filepath = os.path.join(scratch_path, filename)
        if os.path.isdir(filepath):
            shutil.rmtree(filepath)
        else:
            os.remove(filepath)
    args.g2p_path = args.qc_path = create_app_result(args.data_path,
                                                     run_json, suffix='results')
    if run_json.run_id is None:
        run_summary = None
    else:
        logger.info('Summarizing run.')
        run_summary = summarize_run(args, run_json)

    pool = Pool()
    pool.map(functools.partial(try_sample,
                               run_info=run_json,
                               args=args,
                               pssm=pssm),
             range(len(run_json.samples)))

    pool.close()
    pool.join()
    collate_samples(args, run_json)
    if run_json.run_id is not None:
        summarize_samples(args, run_json, run_summary)
    logger.info('Done.')

if __name__ == '__main__':
    main()
