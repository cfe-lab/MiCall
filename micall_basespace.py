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
from micall.core.censor_fastq import censor
from micall.core.filter_quality import report_bad_cycles
from micall.core.remap import remap
from micall.core.prelim_map import prelim_map
from micall.core.sam2aln import sam2aln
from micall.monitor import error_metrics_parser, quality_metrics_parser
from micall.g2p.sam_g2p import sam_g2p, DEFAULT_MIN_COUNT
from micall.g2p.pssm_lib import Pssm
from micall.monitor.tile_metrics_parser import summarize_tiles
from micall.core.coverage_plots import coverage_plot

logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s[%(levelname)s]%(name)s.%(funcName)s(): %(message)s')
logger = logging.getLogger('micall')


class BSrequest:
    """A class that wraps some of the Basespace API"""

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

    def _responsestatus(self, jobj):
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
        logger.info("REQ '%s'" % reqstr)
        return requests.get(reqstr, params=pdct)

    def _runssamples(self, runid, Limit, Offset):
        """ Return the json result of the runs/id/samples API request:
        https://developer.basespace.illumina.com/docs/content/documentation/rest-api/api-reference#Samples
        This returns a json object describing samples  belonging to a specified run, or
        an object with an error message.

        Limit: limit the number of responses to Limit number of items.
        Offset: the offset of the first file to read in the whole collection
        """
        return self._raw_get_file("/".join(["runs", runid, "samples"]),
                                  paramdct={"Limit": Limit, "Offset": Offset}).json()

    def _get_all_sample_ids_from_run_id(self, runid):
        """Retrieve a set of the sampleids of all the samples in the specified run.
        Raise an Runtimeerror iff an API error occurs.

        NOTE: This routine will download sample_id in batches of NUM_BATCH defined in
        the routine itself.
        """
        NUM_BATCH = 1000
        jobj = self._runssamples(runid, NUM_BATCH, 0)
        err_code, err_msg = self._responsestatus(jobj)
        if err_code:
            raise RuntimeError("runsamples API error")
        response = jobj["Response"]
        retset = set([s["Id"] for s in response["Items"]])
        numgot = response["DisplayedCount"]
        numneed = response["TotalCount"]

        while (numgot < numneed):
            jobj = self._runssamples(runid, NUM_BATCH, numgot)
            err_code, err_msg = self._responsestatus(jobj)
            if err_code:
                raise RuntimeError("runsamples API error")
            response = jobj["Response"]
            numgot += response["DisplayedCount"]
            retset |= set([s["Id"] for s in response["Items"]])
        return retset

    def Check_Run_Sample_IDs(self, run_id_lst, sample_id_lst):
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
                                 Items=map(itemgetter('Href'), run_info.samples))
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

    logger.info('Running prelim_map (%d of %d).', sample_index+1, len(run_info.samples))
    with open(os.path.join(sample_scratch_path, 'prelim.csv'), 'wb') as prelim_csv:
        prelim_map(censored_path1,
                   censored_path2,
                   prelim_csv,
                   work_path=sample_scratch_path)

    logger.info('Running remap (%d of %d).', sample_index+1, len(run_info.samples))
    with open(os.path.join(sample_scratch_path, 'prelim.csv'), 'rU') as prelim_csv, \
            open(os.path.join(sample_scratch_path, 'remap.csv'), 'wb') as remap_csv, \
            open(os.path.join(sample_scratch_path, 'remap_counts.csv'), 'wb') as counts_csv, \
            open(os.path.join(sample_scratch_path, 'remap_conseq.csv'), 'wb') as conseq_csv, \
            open(os.path.join(sample_qc_path, 'unmapped1.fastq'), 'w') as unmapped1, \
            open(os.path.join(sample_qc_path, 'unmapped2.fastq'), 'w') as unmapped2:

        remap(censored_path1,
              censored_path2,
              prelim_csv,
              remap_csv,
              counts_csv,
              conseq_csv,
              unmapped1,
              unmapped2,
              sample_scratch_path,
              nthreads=1)

    logger.info('Running sam2aln (%d of %d).', sample_index+1, len(run_info.samples))
    with open(os.path.join(sample_scratch_path, 'remap.csv'), 'rU') as remap_csv, \
            open(os.path.join(sample_scratch_path, 'aligned.csv'), 'wb') as aligned_csv, \
            open(os.path.join(sample_scratch_path, 'conseq_ins.csv'), 'wb') as conseq_ins_csv, \
            open(os.path.join(sample_scratch_path, 'failed_read.csv'), 'wb') as failed_csv, \
            open(os.path.join(sample_scratch_path, 'clipping.csv'), 'wb') as clipping_csv:

        sam2aln(remap_csv,
                aligned_csv,
                conseq_ins_csv,
                failed_csv,
                clipping_csv=clipping_csv)

    logger.info('Running aln2counts (%d of %d).', sample_index+1, len(run_info.samples))
    with open(os.path.join(sample_scratch_path, 'aligned.csv'), 'rU') as aligned_csv, \
            open(os.path.join(sample_scratch_path, 'clipping.csv'), 'rU') as clipping_csv, \
            open(os.path.join(sample_scratch_path, 'conseq_ins.csv'), 'rU') as conseq_ins_csv, \
            open(os.path.join(sample_scratch_path, 'nuc.csv'), 'wb') as nuc_csv, \
            open(os.path.join(sample_scratch_path, 'amino.csv'), 'wb') as amino_csv, \
            open(os.path.join(sample_scratch_path, 'coord_ins.csv'), 'wb') as coord_ins_csv, \
            open(os.path.join(sample_scratch_path, 'conseq.csv'), 'wb') as conseq_csv, \
            open(os.path.join(sample_scratch_path, 'failed_align.csv'), 'wb') as failed_align_csv, \
            open(os.path.join(sample_scratch_path, 'nuc_variants.csv'), 'wb') as nuc_variants_csv, \
            open(os.path.join(sample_scratch_path, 'coverage_summary.csv'), 'wb') as coverage_summary_csv:

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
                      coverage_maps_prefix=sample_name)

    with open(os.path.join(sample_scratch_path, 'coverage_scores.csv'), 'rU') as coverage_scores_csv:
        reader = csv.DictReader(coverage_scores_csv)
        is_v3loop_good = False
        for row in reader:
            if row['region'] == 'V3LOOP':
                is_v3loop_good = row['on.score'] == '4'
                break

    if is_v3loop_good:
        logger.info('Running sam_g2p (%d of %d).', sample_index+1, len(run_info.samples))
        with open(os.path.join(sample_scratch_path, 'remap.csv'), 'rU') as remap_csv, \
                open(os.path.join(sample_scratch_path, 'nuc.csv'), 'rU') as nuc_csv, \
                open(os.path.join(sample_scratch_path, 'g2p.csv'), 'wb') as g2p_csv, \
                open(os.path.join(sample_scratch_path, 'g2p_summary.csv'), 'wb') as g2p_summary_csv:

            sam_g2p(pssm=pssm,
                    remap_csv=remap_csv,
                    nuc_csv=nuc_csv,
                    g2p_csv=g2p_csv,
                    g2p_summary_csv=g2p_summary_csv,
                    min_count=DEFAULT_MIN_COUNT)
    logger.info('Finished sample (%d of %d).', sample_index+1, len(run_info.samples))


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


def summarize_samples(args, json, run_summary):
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


def main():
    logger.info("Starting on %s with %d CPU's.",
                socket.gethostname(),
                multiprocessing.cpu_count())
    args = parse_args()
    if args.link_run is not None:
        json = link_json(args.link_run, args.data_path)
    else:
        json_path = os.path.join(args.data_path, 'input', 'AppSession.json')
        try:
            with open(json_path, 'rU') as json_file:
                json = parse_json(json_file)
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
        if json.run_id is None:
            json.has_runinfo = False
        else:
            bs = BSrequest()
            sample_id_lst = [s["Id"] for s in json.samples]
            sample_id_set = bs.Check_Run_Sample_IDs([json.run_id], sample_id_lst)
            json.has_runinfo = (len(sample_id_set) == len(json.samples))
        logger.info("setting json.has_run_info to %s" % json.has_runinfo)
    pssm = Pssm()

    scratch_path = os.path.join(args.data_path, 'scratch')
    makedirs(scratch_path)
    for filename in os.listdir(scratch_path):
        filepath = os.path.join(scratch_path, filename)
        if os.path.isdir(filepath):
            shutil.rmtree(filepath)
        else:
            os.remove(filepath)
    args.qc_path = create_app_result(args.data_path,
                                     json,
                                     suffix='quality_control')
    args.g2p_path = create_app_result(args.data_path,
                                      json,
                                      suffix='geno_to_pheno')

    if json.run_id is not None:
        logger.info('Summarizing run.')
        run_summary = summarize_run(args, json)

    pool = Pool()
    pool.map(functools.partial(try_sample,
                               run_info=json,
                               args=args,
                               pssm=pssm),
             range(len(json.samples)))

    pool.close()
    pool.join()
    collate_samples(args, json)
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
