from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
import csv
import errno
import fnmatch
import functools
from glob import glob
import json
import logging
import multiprocessing
from multiprocessing.pool import Pool
from operator import attrgetter
import os
import shutil
import socket
from zipfile import ZipFile, ZIP_DEFLATED
import tarfile

from micall.core.filter_quality import report_bad_cycles
from micall.drivers.run_info import RunInfo, ReadSizes, parse_read_sizes
from micall.drivers.sample import Sample
from micall.drivers.sample_group import SampleGroup
from micall.monitor.find_groups import find_groups
from micall.monitor import error_metrics_parser, quality_metrics_parser
from micall.g2p.pssm_lib import Pssm
from micall.monitor.tile_metrics_parser import summarize_tiles

EXCLUDED_SEEDS = ['HLA-B-seed']  # Not ready yet.
EXCLUDED_PROJECTS = ['HCV-NS5a',
                     'HIV-GP41',
                     'HIV',
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
        import requests

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
        logger.info("API: %r" % reqstr)
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
        err_code, err_msg = self._responsestatus(jobj)
        if err_code:
            raise RuntimeError("runsamples API error {}: {}".format(err_code, err_msg))
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
    parser = ArgumentParser(description='Map FASTQ files to references.',
                            formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument(
        "--denovo",
        action="store_true",
        help="Use de novo assembly instead of mapping to reference sequences.",
    )
    parser.add_argument(
        "--input_dir",
        help="Directory to look for input files in (if Dockerized, this is the "
             "path *inside* the container; point a bind mount here).  If input "
             "file paths are not specified as absolute paths, they will be taken "
             "as being relative to this path.",
        default=".",
    )
    parser.add_argument(
        "--data_path",
        help="Directory to write outputs to (if Dockerized, this is the "
             "path *inside* the container; point a bind mount here).  If output "
             "file paths are not specified as absolute paths, they will be taken "
             "as being relative to this path.",
        default=".",
    )

    subparsers = parser.add_subparsers()

    # ####
    # The "process a full BaseSpace-esque directory" subcommand.
    # ####
    full_run_parser = subparsers.add_parser(
        "folder",
        help="Process an entire run directory (e.g. for BaseSpace)",
    )
    full_run_parser.add_argument(
        '--link_run',
        '-l',
        help='Run folder to link into the data folder'
    )
    full_run_parser.add_argument(
        '--all_projects',
        '-a',
        action='store_true',
        help="Don't exclude any projects or seeds."
    )
    full_run_parser.add_argument(
        '--debug_remap',
        '-d',
        action='store_true',
        help="Write debug files for remapping steps."
    )
    full_run_parser.add_argument(
        '--max_active',
        '-m',
        type=int,
        help="Maximum number of samples to process at once, "
             "if not the number of CPU's."
    )
    full_run_parser.set_defaults(func=full_run)

    # ####
    # The "process a single sample" subcommand.
    # ####
    # First, the inputs.
    single_sample_parser = subparsers.add_parser(
        "sample",
        help="Process a single sample"
    )
    single_sample_parser.add_argument(
        "fastq1",
        help="FASTQ file containing forward reads (either an absolute path or relative"
             " to the bind-mounted input directory)"
    )
    single_sample_parser.add_argument(
        "fastq2",
        help="FASTQ file containing reverse reads (either an absolute path or relative"
             " to the bind-mounted input directory)"
    )
    single_sample_parser.add_argument(
        "--bad_cycles_csv",
        help="list of tiles and cycles rejected for poor quality (either an absolute "
             "path or relative to the bind-mounted input directory)"
    )

    # Next, the outputs.
    single_sample_parser.add_argument(
        "--g2p_csv",
        help="CSV containing g2p predictions.",
        default="g2p.csv",
    )
    single_sample_parser.add_argument(
        "--g2p_summary_csv",
        help="CSV containing overall call for the sample.",
        default="g2p_summary.csv",

    )
    single_sample_parser.add_argument(
        "--remap_counts_csv",
        help="CSV containing numbers of mapped reads",
        default="remap_counts.csv",
    )
    single_sample_parser.add_argument(
        "--remap_conseq_csv",
        help="CSV containing mapping consensus sequences",
        default="remap_conseq.csv",
    )
    single_sample_parser.add_argument(
        "--unmapped1_fastq",
        help="FASTQ R1 of reads that failed to map to any region",
        default="unmapped1.fastq",
    )
    single_sample_parser.add_argument(
        "--unmapped2_fastq",
        help="FASTQ R2 of reads that failed to map to any region",
        default="unmapped2.fastq",
    )
    single_sample_parser.add_argument(
        "--conseq_ins_csv",
        help="CSV containing insertions relative to sample consensus",
        default="conseq_ins.csv",
    )
    single_sample_parser.add_argument(
        "--failed_csv",
        help="CSV containing reads that failed to merge",
        default="failed.csv",
    )
    single_sample_parser.add_argument(
        "--cascade_csv",
        help="count of reads at each step",
        default="cascade.csv",
    )
    single_sample_parser.add_argument(
        "--nuc_csv",
        help="CSV containing nucleotide frequencies",
        default="nuc.csv",
    )
    single_sample_parser.add_argument(
        "--amino_csv",
        help="CSV containing amino frequencies",
        default="amino.csv",
    )
    single_sample_parser.add_argument(
        "--coord_ins_csv",
        help="CSV containing insertions relative to coordinate reference",
        default="coord_ins.csv",
    )
    single_sample_parser.add_argument(
        "--conseq_csv",
        help="CSV containing consensus sequences",
        default="conseq.csv",
    )
    single_sample_parser.add_argument(
        "--conseq_region_csv",
        help="CSV containing consensus sequences, split by region",
        default="conseq_region.csv",
    )
    single_sample_parser.add_argument(
        "--failed_align_csv",
        help="CSV containing any consensus that failed to align",
        default="failed_align.csv",
    )
    single_sample_parser.add_argument(
        "--coverage_scores_csv",
        help="CSV coverage scores.",
        default="coverage_scores.csv",
    )
    single_sample_parser.add_argument(
        "--coverage_maps_tar",
        help="tar file of coverage maps.",
        default="coverage_maps.tar",
    )
    single_sample_parser.add_argument(
        "--aligned_csv",
        help="CSV containing individual reads aligned to consensus",
        default="aligned.csv",
    )
    single_sample_parser.add_argument(
        "--g2p_aligned_csv",
        help="CSV containing individual reads aligned to V3LOOP",
        default="g2p_aligned.csv",
    )
    single_sample_parser.add_argument(
        "--genome_coverage_csv",
        help="CSV of coverage levels in full-genome coordinates",
        default="genome_coverage.csv",
    )
    single_sample_parser.add_argument(
        "--genome_coverage_svg",
        help="SVG diagram of coverage in full-genome coordinates",
        default="genome_coverage.svg",
    )
    single_sample_parser.add_argument(
        "--contigs_csv",
        help="CSV containing contigs built by de novo assembly",
        default="contigs.csv",
    )
    single_sample_parser.add_argument(
        "--read_entropy_csv",
        help="CSV containing read pair length counts",
        default="read_entropy.csv",
    )
    single_sample_parser.set_defaults(func=single_sample)

    # ####
    # The "process a single HCV sample" subcommand.
    # ####
    # First, the inputs.
    hcv_sample_parser = subparsers.add_parser(
        "hcv_sample",
        help="Process a single HCV sample"
    )
    hcv_sample_parser.add_argument(
        "fastq1",
        help="FASTQ file containing forward reads (either an absolute path or relative"
             " to the bind-mounted input directory)"
    )
    hcv_sample_parser.add_argument(
        "fastq2",
        help="FASTQ file containing reverse reads (either an absolute path or relative"
             " to the bind-mounted input directory)"
    )
    hcv_sample_parser.add_argument(
        "midi1",
        help="FASTQ file containing HCV MIDI forward reads (either an absolute path "
             "or relative to the bind-mounted input directory)"
    )
    hcv_sample_parser.add_argument(
        "midi2",
        help="FASTQ file containing HCV MIDI reverse reads (either an absolute path "
             "or relative to the bind-mounted input directory)"
    )
    hcv_sample_parser.add_argument(
        "--bad_cycles_csv",
        help="list of tiles and cycles rejected for poor quality (either an absolute "
             "path or relative to the bind-mounted input directory)"
    )
    hcv_sample_parser.add_argument(
        "--midi_bad_cycles_csv",
        help="list of tiles and cycles rejected for poor quality (either an absolute "
             "path or relative to the bind-mounted input directory)"
    )

    # Next, the outputs pertaining to the main samples.
    hcv_sample_parser.add_argument(
        "--g2p_csv",
        help="CSV containing g2p predictions.",
        default="g2p.csv",
    )
    hcv_sample_parser.add_argument(
        "--g2p_summary_csv",
        help="CSV containing overall call for the sample.",
        default="g2p_summary.csv",

    )
    hcv_sample_parser.add_argument(
        "--remap_counts_csv",
        help="CSV containing numbers of mapped reads",
        default="remap_counts.csv",
    )
    hcv_sample_parser.add_argument(
        "--remap_conseq_csv",
        help="CSV containing mapping consensus sequences",
        default="remap_conseq.csv",
    )
    hcv_sample_parser.add_argument(
        "--unmapped1_fastq",
        help="FASTQ R1 of reads that failed to map to any region",
        default="unmapped1.fastq",
    )
    hcv_sample_parser.add_argument(
        "--unmapped2_fastq",
        help="FASTQ R2 of reads that failed to map to any region",
        default="unmapped2.fastq",
    )
    hcv_sample_parser.add_argument(
        "--conseq_ins_csv",
        help="CSV containing insertions relative to sample consensus",
        default="conseq_ins.csv",
    )
    hcv_sample_parser.add_argument(
        "--failed_csv",
        help="CSV containing reads that failed to merge",
        default="failed.csv",
    )
    hcv_sample_parser.add_argument(
        "--cascade_csv",
        help="count of reads at each step",
        default="cascade.csv",
    )
    hcv_sample_parser.add_argument(
        "--nuc_csv",
        help="CSV containing nucleotide frequencies",
        default="nuc.csv",
    )
    hcv_sample_parser.add_argument(
        "--amino_csv",
        help="CSV containing amino frequencies",
        default="amino.csv",
    )
    hcv_sample_parser.add_argument(
        "--coord_ins_csv",
        help="CSV containing insertions relative to coordinate reference",
        default="coord_ins.csv",
    )
    hcv_sample_parser.add_argument(
        "--conseq_csv",
        help="CSV containing consensus sequences",
        default="conseq.csv",
    )
    hcv_sample_parser.add_argument(
        "--conseq_region_csv",
        help="CSV containing consensus sequences, split by region",
        default="conseq_region.csv",
    )
    hcv_sample_parser.add_argument(
        "--failed_align_csv",
        help="CSV containing any consensus that failed to align",
        default="failed_align.csv",
    )
    hcv_sample_parser.add_argument(
        "--coverage_scores_csv",
        help="CSV coverage scores.",
        default="coverage_scores.csv",
    )
    hcv_sample_parser.add_argument(
        "--coverage_maps_tar",
        help="tar file of coverage maps.",
        default="coverage_maps.tar",
    )
    hcv_sample_parser.add_argument(
        "--aligned_csv",
        help="CSV containing individual reads aligned to consensus",
        default="aligned.csv",
    )
    hcv_sample_parser.add_argument(
        "--g2p_aligned_csv",
        help="CSV containing individual reads aligned to V3LOOP",
        default="g2p_aligned.csv",
    )
    hcv_sample_parser.add_argument(
        "--genome_coverage_csv",
        help="CSV of coverage levels in full-genome coordinates",
        default="genome_coverage.csv",
    )
    hcv_sample_parser.add_argument(
        "--genome_coverage_svg",
        help="SVG diagram of coverage in full-genome coordinates",
        default="genome_coverage.svg",
    )
    hcv_sample_parser.add_argument(
        "--contigs_csv",
        help="CSV containing contigs built by de novo assembly",
        default="contigs.csv",
    )
    hcv_sample_parser.add_argument(
        "--read_entropy_csv",
        help="CSV containing read pair length counts",
        default="read_entropy.csv",
    )

    # Next, outputs pertaining to the MIDI samples.
    hcv_sample_parser.add_argument(
        "--midi_g2p_csv",
        help="CSV containing g2p predictions (MIDI).",
        default="midi_g2p.csv",
    )
    hcv_sample_parser.add_argument(
        "--midi_g2p_summary_csv",
        help="CSV containing overall call for the sample (MIDI).",
        default="midi_g2p_summary.csv",

    )
    hcv_sample_parser.add_argument(
        "--midi_remap_counts_csv",
        help="CSV containing numbers of mapped reads (MIDI)",
        default="midi_remap_counts.csv",
    )
    hcv_sample_parser.add_argument(
        "--midi_remap_conseq_csv",
        help="CSV containing mapping consensus sequences (MIDI)",
        default="midi_remap_conseq.csv",
    )
    hcv_sample_parser.add_argument(
        "--midi_unmapped1_fastq",
        help="FASTQ R1 of reads that failed to map to any region (MIDI)",
        default="midi_unmapped1.fastq",
    )
    hcv_sample_parser.add_argument(
        "--midi_unmapped2_fastq",
        help="FASTQ R2 of reads that failed to map to any region (MIDI)",
        default="midi_unmapped2.fastq",
    )
    hcv_sample_parser.add_argument(
        "--midi_conseq_ins_csv",
        help="CSV containing insertions relative to sample consensus (MIDI)",
        default="midi_conseq_ins.csv",
    )
    hcv_sample_parser.add_argument(
        "--midi_failed_csv",
        help="CSV containing reads that failed to merge (MIDI)",
        default="midi_failed.csv",
    )
    hcv_sample_parser.add_argument(
        "--midi_cascade_csv",
        help="count of reads at each step (MIDI)",
        default="midi_cascade.csv",
    )
    hcv_sample_parser.add_argument(
        "--midi_nuc_csv",
        help="CSV containing nucleotide frequencies (MIDI)",
        default="midi_nuc.csv",
    )
    hcv_sample_parser.add_argument(
        "--midi_amino_csv",
        help="CSV containing amino frequencies (MIDI)",
        default="midi_amino.csv",
    )
    hcv_sample_parser.add_argument(
        "--midi_coord_ins_csv",
        help="CSV containing insertions relative to coordinate reference (MIDI)",
        default="midi_coord_ins.csv",
    )
    hcv_sample_parser.add_argument(
        "--midi_conseq_csv",
        help="CSV containing consensus sequences (MIDI)",
        default="midi_conseq.csv",
    )
    hcv_sample_parser.add_argument(
        "--midi_conseq_region_csv",
        help="CSV containing consensus sequences, split by region (MIDI)",
        default="midi_conseq_region.csv",
    )
    hcv_sample_parser.add_argument(
        "--midi_failed_align_csv",
        help="CSV containing any consensus that failed to align (MIDI)",
        default="midi_failed_align.csv",
    )
    hcv_sample_parser.add_argument(
        "--midi_coverage_scores_csv",
        help="CSV coverage scores (MIDI).",
        default="midi_coverage_scores.csv",
    )
    hcv_sample_parser.add_argument(
        "--midi_coverage_maps_tar",
        help="tar file of coverage maps (MIDI).",
        default="midi_coverage_maps.tar",
    )
    hcv_sample_parser.add_argument(
        "--midi_aligned_csv",
        help="CSV containing individual reads aligned to consensus (MIDI)",
        default="midi_aligned.csv",
    )
    hcv_sample_parser.add_argument(
        "--midi_g2p_aligned_csv",
        help="CSV containing individual reads aligned to V3LOOP (MIDI)",
        default="midi_g2p_aligned.csv",
    )
    hcv_sample_parser.add_argument(
        "--midi_genome_coverage_csv",
        help="CSV of coverage levels in full-genome coordinates (MIDI)",
        default="midi_genome_coverage.csv",
    )
    hcv_sample_parser.add_argument(
        "--midi_genome_coverage_svg",
        help="SVG diagram of coverage in full-genome coordinates (MIDI)",
        default="midi_genome_coverage.svg",
    )
    hcv_sample_parser.add_argument(
        "--midi_contigs_csv",
        help="CSV containing contigs built by de novo assembly (MIDI)",
        default="midi_contigs.csv",
    )
    hcv_sample_parser.add_argument(
        "--midi_read_entropy_csv",
        help="CSV containing read pair length counts (MIDI)",
        default="midi_read_entropy.csv",
    )

    hcv_sample_parser.add_argument(
        "--resistance_csv",
        help="CSV containing resistance calls.",
        default="resistance.csv",
    )
    hcv_sample_parser.add_argument(
        "--mutations_csv",
        help="CSV containing resistance mutations.",
        default="mutations.csv",
    )
    hcv_sample_parser.add_argument(
        "--resistance_fail_csv",
        help="CSV containing failure reasons.",
        default="resistance_fail.csv",
    )
    hcv_sample_parser.add_argument(
        "--resistance_pdf",
        help="resistance report",
        default="resistance.pdf",
    )
    hcv_sample_parser.add_argument(
        "--resistance_consensus_csv",
        help="CSV with amino consensus used for resistance.",
        default="resistance_consensus.csv",
    )

    hcv_sample_parser.set_defaults(func=hcv_sample)

    return parser.parse_args()


class MiCallArgs:
    """
    Wrapper that performs some path translation on our inputs and outputs.
    """
    INPUTS = [
        "fastq1",
        "fastq2",
        "bad_cycles_csv",
        "midi1",
        "midi2",
        "midi_bad_cycles_csv"
    ]

    def __init__(self, args):
        self.original_args = args

    def __getattr__(self, arg_name):
        if arg_name.startswith('__'):
            raise AttributeError(arg_name)
        resolved_path = self.original_args.get(arg_name)
        if not os.path.isabs(resolved_path):
            io_prefix = (self.original_args.input_dir if arg_name in self.INPUTS
                         else self.original_args.data_path)
            resolved_path = os.path.join(io_prefix, resolved_path)
        return resolved_path


def full_run(args):
    resolved_args = MiCallArgs(args)
    if args.link_run is not None:
        run_info = link_samples(resolved_args.link_run, resolved_args.data_path, args.denovo)
    elif os.path.exists(resolved_args.data_path, "input", "AppSession.json"):
        run_info = load_samples(resolved_args.data_path)
    elif os.path.exists(os.path.join(args.input_dir, "SampleSheet.csv")):
        run_info = link_samples(args.input_dir, resolved_args.data_path, args.denovo)
    else:
        raise ValueError("Could not find the required inputs to process.")
    pssm = Pssm()

    for filename in os.listdir(run_info.scratch_path):
        filepath = os.path.join(run_info.scratch_path, filename)
        if os.path.isdir(filepath):
            shutil.rmtree(filepath)
        else:
            os.remove(filepath)

    if run_info.interop_path is None:
        run_summary = None
    else:
        logger.info('Summarizing run.')
        run_summary = summarize_run(run_info)

    pool = Pool(processes=args.max_active)
    pool.map(functools.partial(process_sample,
                               args=args,
                               pssm=pssm,
                               use_denovo=run_info.is_denovo),
             run_info.get_all_samples())

    pool.close()
    pool.join()
    pool = Pool()
    pool.map(functools.partial(process_resistance,
                               run_info=run_info),
             run_info.sample_groups)

    pool.close()
    pool.join()
    collate_samples(run_info)
    if run_summary is not None:
        summarize_samples(run_info, run_summary)
    logger.info('Done.')


def single_sample(args):
    resolved_args = MiCallArgs(args)
    scratch_path = os.path.join(os.path.dirname(resolved_args.cascade_csv), "scratch")
    shutil.rmtree(scratch_path, ignore_errors=True)

    sample = Sample(
        fastq1=resolved_args.fastq1,
        fastq2=resolved_args.fastq2,
        bad_cycles_csv=resolved_args.bad_cycles_csv,
        g2p_csv=resolved_args.g2p_csv,
        g2p_summary_csv=resolved_args.g2p_summary_csv,
        remap_counts_csv=resolved_args.remap_counts_csv,
        remap_conseq_csv=resolved_args.remap_conseq_csv,
        unmapped1_fastq=resolved_args.unmapped1_fastq,
        unmapped2_fastq=resolved_args.unmapped2_fastq,
        conseq_ins_csv=resolved_args.conseq_ins_csv,
        failed_csv=resolved_args.failed_csv,
        cascade_csv=resolved_args.cascade_csv,
        nuc_csv=resolved_args.nuc_csv,
        amino_csv=resolved_args.amino_csv,
        coord_ins_csv=resolved_args.coord_ins_csv,
        conseq_csv=resolved_args.conseq_csv,
        conseq_region_csv=resolved_args.conseq_region_csv,
        failed_align_csv=resolved_args.failed_align_csv,
        coverage_scores_csv=resolved_args.coverage_scores_csv,
        aligned_csv=resolved_args.aligned_csv,
        g2p_aligned_csv=resolved_args.g2p_aligned_csv,
        contigs_csv=resolved_args.contigs_csv,
        genome_coverage_csv=resolved_args.genome_coverage_csv,
        genome_coverage_svg=resolved_args.genome_coverage_svg,
        read_entropy_csv=resolved_args.read_entropy_csv,
        scratch_path=scratch_path
    )

    pssm = Pssm()
    sample.process(pssm,
                   force_gzip=True,  # dataset files change .gz to .raw
                   use_denovo=args.denovo)

    with tarfile.open(resolved_args.coverage_maps_tar, mode='w') as tar:
        for image_name in os.listdir(sample.coverage_maps):
            image_path = os.path.join(sample.coverage_maps, image_name)
            archive_path = os.path.join('coverage_maps', image_name)
            tar.add(image_path, archive_path)

    return sample


def hcv_sample(args):
    # First, process the main samples.
    single_sample(args)

    # Do the same for the MIDI samples.
    resolved_args = MiCallArgs(args)
    scratch_path = os.path.join(os.path.dirname(resolved_args.midi_cascade_csv), "scratch")
    shutil.rmtree(scratch_path, ignore_errors=True)

    # First process the main samples.
    midi_sample = Sample(
        fastq1=resolved_args.midi_fastq1,
        fastq2=resolved_args.midi_fastq2,
        bad_cycles_csv=resolved_args.midi_bad_cycles_csv,
        g2p_csv=resolved_args.midi_g2p_csv,
        g2p_summary_csv=resolved_args.midi_g2p_summary_csv,
        remap_counts_csv=resolved_args.midi_remap_counts_csv,
        remap_conseq_csv=resolved_args.midi_remap_conseq_csv,
        unmapped1_fastq=resolved_args.midi_unmapped1_fastq,
        unmapped2_fastq=resolved_args.midi_unmapped2_fastq,
        conseq_ins_csv=resolved_args.midi_conseq_ins_csv,
        failed_csv=resolved_args.midi_failed_csv,
        cascade_csv=resolved_args.midi_cascade_csv,
        nuc_csv=resolved_args.midi_nuc_csv,
        amino_csv=resolved_args.midi_amino_csv,
        coord_ins_csv=resolved_args.midi_coord_ins_csv,
        conseq_csv=resolved_args.midi_conseq_csv,
        conseq_region_csv=resolved_args.midi_conseq_region_csv,
        failed_align_csv=resolved_args.midi_failed_align_csv,
        coverage_scores_csv=resolved_args.midi_coverage_scores_csv,
        aligned_csv=resolved_args.midi_aligned_csv,
        g2p_aligned_csv=resolved_args.midi_g2p_aligned_csv,
        contigs_csv=resolved_args.midi_contigs_csv,
        genome_coverage_csv=resolved_args.midi_genome_coverage_csv,
        genome_coverage_svg=resolved_args.midi_genome_coverage_svg,
        read_entropy_csv=resolved_args.midi_read_entropy_csv,
        scratch_path=scratch_path
    )

    pssm = Pssm()
    midi_sample.process(
        pssm,
        force_gzip=True,  # dataset files change .gz to .raw
        use_denovo=args.denovo,
    )

    with tarfile.open(resolved_args.midi_coverage_maps_tar, mode='w') as tar:
        for image_name in os.listdir(midi_sample.coverage_maps):
            image_path = os.path.join(midi_sample.coverage_maps, image_name)
            archive_path = os.path.join('coverage_maps', image_name)
            tar.add(image_path, archive_path)

    # Now, analyze the two samples together for resistance.
    # scratch_path = os.path.join(os.path.dirname(resolved_args.main_amino_csv), 'scratch')
    # shutil.rmtree(scratch_path, ignore_errors=True)
    sample1 = Sample(
        amino_csv=resolved_args.amino_csv,
        resistance_csv=resolved_args.resistance_csv,
        mutations_csv=resolved_args.mutations_csv,
        resistance_fail_csv=resolved_args.resistance_fail_csv,
        resistance_pdf=resolved_args.resistance_pdf,
        resistance_consensus_csv=resolved_args.resistance_consensus_csv,
        scratch_path=scratch_path
    )
    sample2 = Sample(amino_csv=resolved_args.midi_amino_csv)
    main_and_midi = SampleGroup(sample1, sample2)

    main_and_midi.process_resistance(RunInfo([main_and_midi]))
    return main_and_midi


class Args(object):
    pass


def load_samples(data_path):
    """ Load JSON file from the data path, and pull out the arguments for this run.

    :param str data_path: folder that contains a JSON file in the BaseSpace
    AppSession format.
    :return RunInfo: details about the run and samples
    """
    json_path = os.path.join(data_path, 'input', 'AppSession.json')
    try:
        with open(json_path, 'r') as json_file:
            raw_args = json.load(json_file)

        arg_map = {item['Name']: item
                   for item in raw_args['Properties']['Items']}

        href_app_session = raw_args['Href']
        run = arg_map.get('Input.run-id')
        if run is None:
            run_id = interop_path = read_sizes = None
        else:
            run_content = run['Content']
            run_id = run_content['Id']
            interop_path = os.path.join(data_path,
                                        'input',
                                        'runs',
                                        run_id,
                                        'InterOp')
            read_sizes = ReadSizes(
                run_content['SequencingStats']['NumCyclesRead1'],
                run_content['SequencingStats']['NumCyclesRead2'],
                run_content['SequencingStats']['NumCyclesIndex1'],
                run_content['SequencingStats']['NumCyclesIndex2'])
        project_id = arg_map['Input.project-id']['Content']['Id']
        output_path = os.path.join(data_path,
                                   'output',
                                   'appresults',
                                   project_id,
                                   'results')
        makedirs(output_path)
        reports = arg_map['Input.reports']['Items']
        builder_node = arg_map.get('Input.builder')
        if builder_node is None:
            is_denovo = False
        else:
            is_denovo = builder_node['Content'] == 'denovo'

        scratch_path = os.path.join(data_path, 'scratch')
        sample_groups = []
        run_info = RunInfo(sample_groups,
                           reports,
                           interop_path,
                           scratch_path,
                           output_path,
                           read_sizes,
                           href_app_session,
                           is_denovo)
        main_samples = arg_map['Input.sample-ids.main']['Items']
        midi_samples = arg_map['Input.sample-ids.midi']['Items']
        for main_sample_json, midi_sample_json in zip(main_samples, midi_samples):
            sample_group = SampleGroup(load_sample(main_sample_json,
                                                   data_path,
                                                   scratch_path),
                                       load_sample(midi_sample_json,
                                                   data_path,
                                                   scratch_path))
            sample_groups.append(sample_group)

        # Do we have run_ids for all sample_ids ?
        if run_id is not None:
            bs = BSrequest()
            all_ids = {s.basespace_id for s in run_info.get_all_samples()}
            sample_id_set = bs.check_run_sample_ids(
                [run_id],
                all_ids)
            if len(sample_id_set) != len(all_ids):
                for s in run_info.get_all_samples():
                    if s.basespace_id not in sample_id_set:
                        logger.warning(
                            'Run info not found for %s, skipping error rate data.',
                            s)
                run_info.read_sizes = run_info.interop_path = None

        create_app_result(run_info)
    except IOError:
        if os.path.exists(json_path):
            # copy the input file to the output dir for postmortem analysis
            logger.error("Error occurred while parsing %r.", json_path)
            with open(json_path, 'r') as json_file:
                file_cont = json_file.read()
            out_path = os.path.join(data_path, 'logs', 'AppSession.json')
            with open(out_path, 'w') as json_file:
                json_file.write(file_cont)
        else:
            logger.error("Error: no such file as %r.", json_path)
        raise

    return run_info


def load_sample(sample_json, data_path, scratch_path):
    if sample_json is None:
        return None
    sample_dir = os.path.join(data_path,
                              'input',
                              'samples',
                              sample_json['Id'],
                              'Data',
                              'Intensities',
                              'BaseCalls')
    if not os.path.exists(sample_dir):
        sample_dir = os.path.join(data_path,
                                  'input',
                                  'samples',
                                  sample_json['Id'])
    fastq1 = None
    for root, _dirs, files in os.walk(sample_dir):
        fastq_files = fnmatch.filter(files, '*_R1_*')
        if fastq_files:
            fastq1 = os.path.join(root, fastq_files[0])
            break
    if fastq1 is None:
        raise RuntimeError(
            'No R1 file found for sample id {}.'.format(sample_json['Id']))
    sample = Sample(fastq1=fastq1,
                    basespace_id=sample_json['Id'],
                    basespace_href=sample_json['Href'])
    sample.scratch_path = os.path.join(scratch_path, sample.name)
    return sample


def link_samples(run_path: str, data_path: str, is_denovo: bool):
    """ Load the data from a run folder into the BaseSpace layout. """

    shutil.rmtree(data_path, ignore_errors=True)
    makedirs(data_path)

    results_dir = os.path.join(run_path, "Results")
    basespace_dir = os.path.join(results_dir, "basespace")
    output_path = os.path.join(data_path, 'output')

    results_dir_exists = os.path.exists(results_dir)
    basespace_dir_exists = os.path.exists(basespace_dir)
    makedirs(basespace_dir)
    try:
        os.link(basespace_dir, output_path)
    except OSError as e:
        if e.errno == errno.EXDEV:
            logger.info(
                "Unable to link %s to %s; writing results directly to the latter.",
                basespace_dir,
                output_path,
            )
            if not results_dir_exists:
                shutil.rmtree(results_dir)
            elif not basespace_dir_exists:
                shutil.rmtree(basespace_dir)
            makedirs(output_path)
        else:
            raise
    scratch_path = os.path.join(data_path, 'scratch')
    makedirs(scratch_path)

    sample_groups = []
    run_info_path = os.path.join(run_path, 'RunInfo.xml')
    interop_path = os.path.join(run_path, 'InterOp')
    if not (os.path.exists(run_info_path) and os.path.exists(interop_path)):
        read_sizes = None
    else:
        read_sizes = parse_read_sizes(run_info_path)
    run_info = RunInfo(sample_groups,
                       reports=['PR_RT', 'IN', 'NS3', 'NS5a', 'NS5b'],
                       interop_path=interop_path,
                       scratch_path=scratch_path,
                       output_path=output_path,
                       read_sizes=read_sizes,
                       is_denovo=is_denovo)

    fastq_files = list(glob(os.path.join(run_path,
                                         'Data',
                                         'Intensities',
                                         'BaseCalls',
                                         '*_R1_*')) or
                       glob(os.path.join(run_path,
                                         '*_R1_*')))
    source_folder = fastq_files and os.path.dirname(fastq_files[0])
    file_names = [os.path.basename(fastq_file) for fastq_file in fastq_files]
    groups = find_groups(file_names,
                         os.path.join(run_path, 'SampleSheet.csv'))
    for group in groups:
        main_file, midi_file = group.names
        if main_file.startswith('Undetermined'):
            continue
        main_sample = Sample(fastq1=os.path.join(source_folder, main_file))
        if midi_file is None:
            midi_sample = None
        else:
            midi_sample = Sample(fastq1=os.path.join(source_folder, midi_file))
        sample_groups.append(SampleGroup(main_sample, midi_sample))

    sample_count = sum(1 for _ in run_info.get_all_samples())
    for i, sample in enumerate(run_info.get_all_samples(), 1):
        sample.rank = '{} of {}'.format(i, sample_count)
        sample.bad_cycles_csv = run_info.bad_cycles_csv
        sample.scratch_path = os.path.join(scratch_path, sample.name)

    return run_info


def create_app_result(run_info,
                      description=''):
    """ Create the metadata for BaseSpace output files.

    :param RunInfo run_info: details of the run
    :param str description: description text for the run
    """
    sample_properties = dict(Type='sample[]',
                             Name='Input.Samples',
                             Items=list(map(attrgetter('basespace_href'),
                                            run_info.get_all_samples())))
    metadata = dict(Name=os.path.basename(run_info.output_path),
                    Description=description,
                    HrefAppSession=run_info.href_app_session,
                    Properties=sample_properties)
    metadata_path = os.path.join(run_info.output_path, '_metadata.json')
    with open(metadata_path, 'w') as json_file:
        json.dump(metadata, json_file, indent=4)


def process_sample(sample, args, pssm, use_denovo=False):
    """ Process a single sample.

    :param Sample sample: the sample to process
    :param args: the command-line arguments
    :param pssm: the pssm library for running G2P analysis
    :param use_denovo: use denovo assembly instead of mapping to references
    """
    sample.debug_remap = args.debug_remap
    try:
        excluded_seeds = [] if args.all_projects else EXCLUDED_SEEDS
        excluded_projects = [] if args.all_projects else EXCLUDED_PROJECTS
        sample.process(pssm,
                       excluded_seeds,
                       excluded_projects,
                       use_denovo=use_denovo)
    except Exception:
        message = 'Failed to process {}.'.format(sample)
        logger.error(message, exc_info=True)
        raise RuntimeError(message)


def process_resistance(sample_group, run_info):
    """ Process the resistance steps for a sample group.

    :param SampleGroup sample_group: the group to process
    :param RunInfo run_info: the details about the whole run
    """
    try:
        sample_group.process_resistance(run_info)
    except Exception:
        message = 'Failed to process resistance of {}.'.format(
            sample_group.main_sample)
        logger.error(message, exc_info=True)
        raise RuntimeError(message)


def summarize_run(run_info):
    """ Summarize the run data from the InterOp folder.

    Writes some summary files.
    :param RunInfo run_info: details of the run
    :return: a dictionary with summary values.
    """
    summary = {}

    if run_info.read_sizes is not None:
        read_lengths = [run_info.read_sizes.read1,
                        run_info.read_sizes.index1,
                        run_info.read_sizes.index2,
                        run_info.read_sizes.read2]
        phix_path = os.path.join(run_info.interop_path, 'ErrorMetricsOut.bin')
        with open(phix_path, 'rb') as phix, \
                open(run_info.quality_csv, 'w') as quality:
            records = error_metrics_parser.read_errors(phix)
            error_metrics_parser.write_phix_csv(quality,
                                                records,
                                                read_lengths,
                                                summary)
        with open(run_info.quality_csv) as quality, \
                open(run_info.bad_cycles_csv, 'w') as bad_cycles, \
                open(run_info.bad_tiles_csv, 'w') as bad_tiles:
            report_bad_cycles(quality, bad_cycles, bad_tiles)

        quality_metrics_path = os.path.join(run_info.interop_path,
                                            'QMetricsOut.bin')
        quality_metrics_parser.summarize_quality(quality_metrics_path,
                                                 summary,
                                                 read_lengths)

        tile_metrics_path = os.path.join(run_info.interop_path,
                                         'TileMetricsOut.bin')
        summarize_tiles(tile_metrics_path, summary)
    return summary


def summarize_samples(run_info, run_summary):
    score_sum = 0.0
    base_count = 0
    coverage_sum = 0.0
    coverage_count = 0
    for sample in run_info.get_all_samples():
        with open(sample.read_summary_csv, 'r') as read_summary:
            reader = csv.DictReader(read_summary)
            for row in reader:
                sample_base_count = int(row['base_count'])
                if sample_base_count:
                    score_sum += float(row['avg_quality']) * sample_base_count
                    base_count += sample_base_count
        with open(sample.coverage_summary_csv, 'r') as coverage_summary:
            reader = csv.DictReader(coverage_summary)
            for row in reader:
                region_width = int(row['region_width'])
                coverage_sum += float(row['avg_coverage']) * region_width
                coverage_count += region_width

    if base_count > 0:
        run_summary['avg_quality'] = score_sum / base_count
    if coverage_count > 0:
        run_summary['avg_coverage'] = coverage_sum / coverage_count

    run_quality_path = os.path.join(run_info.output_path, 'run_quality.csv')
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


def zip_folder(parent_path, folder_name):
    source_path = os.path.join(parent_path, folder_name)
    zip_filename = os.path.join(parent_path, folder_name + '.zip')
    with ZipFile(zip_filename, mode='w', compression=ZIP_DEFLATED) as zip_file:
        for dirpath, dirnames, filenames in os.walk(source_path):
            relpath = os.path.relpath(dirpath, source_path)
            for filename in filenames:
                zip_file.write(os.path.join(dirpath, filename),
                               os.path.join(folder_name, relpath, filename))


def collate_samples(run_info):
    """ Combine all the sample files into run files.

    :param RunInfo run_info: details of the run and samples
    """
    filenames = ['remap_counts.csv',
                 'remap_conseq.csv',
                 'conseq_ins.csv',
                 'failed_read.csv',
                 'nuc.csv',
                 'amino.csv',
                 'coord_ins.csv',
                 'conseq.csv',
                 'conseq_region.csv',
                 'failed_align.csv',
                 'coverage_scores.csv',
                 'g2p.csv',
                 'g2p_summary.csv',
                 'resistance.csv',
                 'mutations.csv',
                 'resistance_fail.csv',
                 'resistance_consensus.csv',
                 'cascade.csv',
                 'merge_lengths.csv']
    for filename in filenames:
        out_path = run_info.output_path
        with open(os.path.join(out_path, filename), 'w') as fout:
            writer = csv.writer(fout, lineterminator=os.linesep)
            is_header_written = False
            for sample_info in run_info.get_all_samples():
                sample_name = sample_info.name
                sample_scratch_path = sample_info.scratch_path
                srcfile = os.path.join(sample_scratch_path, filename)
                try:
                    with open(srcfile, 'r') as fin:
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
    resistance_reports_path = os.path.join(run_info.output_path,
                                           'resistance_reports')
    makedirs(resistance_reports_path)
    coverage_maps_path = os.path.join(run_info.output_path, 'coverage_maps')
    genome_coverage_path = os.path.join(coverage_maps_path, 'genome')
    makedirs(genome_coverage_path)
    merge_lengths_path = os.path.join(run_info.output_path, 'merge_lengths')
    makedirs(merge_lengths_path)
    for sample_info in run_info.get_all_samples():
        if os.path.exists(sample_info.coverage_maps):
            for map_file in os.listdir(sample_info.coverage_maps):
                safe_file_move(os.path.join(sample_info.coverage_maps, map_file),
                               os.path.join(coverage_maps_path, map_file))
        if os.path.exists(sample_info.contigs_svg):
            safe_file_move(sample_info.contigs_svg,
                           os.path.join(coverage_maps_path,
                                        sample_info.name + '_contigs.svg'))
        if os.path.exists(sample_info.genome_coverage_svg):
            safe_file_move(sample_info.genome_coverage_svg,
                           os.path.join(genome_coverage_path,
                                        sample_info.name + '_genome_coverage.svg'))
        if os.path.exists(sample_info.merge_lengths_svg):
            safe_file_move(sample_info.merge_lengths_svg,
                           os.path.join(merge_lengths_path,
                                        sample_info.name + '_merge_lengths.svg'))
        if os.path.exists(sample_info.resistance_pdf):
            safe_file_move(sample_info.resistance_pdf,
                           os.path.join(resistance_reports_path,
                                        sample_info.name + '_resistance.pdf'))
    try:
        # Remove directory, if it's empty.
        os.rmdir(genome_coverage_path)
    except OSError:
        # Guess it wasn't empty.
        pass
    zip_folder(run_info.output_path, 'resistance_reports')
    zip_folder(run_info.output_path, 'coverage_maps')


def safe_file_move(src, dst):
    """
    Helper that attempts to move a file from src to dst.

    Because os.rename may fail on certain platforms, we fall back to
    the boneheaded way of copying-and-deleting.
    :param src:
    :param dst:
    :return:
    """
    try:
        os.rename(src, dst)
    except OSError as e:
        if e.errno == errno.EXDEV:
            logger.debug(
                "Failed to rename %s to %s; copying and deleting the original.",
                src,
                dst,
            )
            shutil.copy(src, dst)
            os.remove(src)
        else:
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
    args.func(args)


if __name__ == '__main__':
    main()
