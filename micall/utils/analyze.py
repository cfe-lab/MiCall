#! /usr/bin/env python

"""
Entry script that serves as an entry point of MiCall's Docker image.
"""

from argparse import ArgumentParser
import csv
import errno
import functools
from concurrent.futures.process import ProcessPoolExecutor
from glob import glob
import json
import logging
from operator import attrgetter
import os
import shutil
import socket
from zipfile import ZipFile, ZIP_DEFLATED
import typing
from pathlib import Path

from micall.core.filter_quality import report_bad_cycles
from micall.core.trim_fastqs import TrimSteps
from micall.drivers.run_info import RunInfo, ReadSizes, parse_read_sizes
from micall.drivers.sample import Sample
from micall.drivers.sample_group import SampleGroup, load_git_version
from micall.monitor.find_groups import find_groups
from micall.utils.interop_wrappers import summarize_quality, summarize_tiles
from micall.g2p.pssm_lib import Pssm
from micall.utils.list_fastq_files import list_fastq_files
from micall.utils.driver_utils import MiCallFormatter, safe_file_move, makedirs, \
    MiCallArgs
from miseqinteropreader.interop_reader import InterOpReader, MetricFile
from miseqinteropreader.error_metrics_parser import write_phix_csv
from miseqinteropreader import ReadLengths4

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
# logging.getLogger('micall.core.trim_fastqs').setLevel(logging.DEBUG)


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


MAIN_DESCRIPTION = """\
Map FASTQ files to references.

This driver can run in several different modes; for more specific help, 
see the help for the different subcommands.

This will typically be Dockerized so that the container itself may be
invoked as if it is the driver: you can invoke it as

    docker run [image name] [subcommand] [arguments]

All input and output paths will be handled as if they are inside the container,
so you will have to set up bind mounts for your inputs and outputs.
"""

BASESPACE_DESCRIPTION = """\
Map FASTQ files to references for a BaseSpace run.

This will look for MiSeq output in the folder /data, and a BaseSpace-generated
file at "/data/input/AppSession.json"; output will be written to
"/data/output/".  This should always be run in a Docker container with the
appropriate bind mounts to "/data", as done by BaseSpace. 

For example: if you built your container as "micall:v0.1.2-3" and your data is
on your host machine as "/path/on/host/", you would invoke it as follows:

    docker run \\
        --mount type=bind,source=/path/on/host/,destination=/data \\
        micall:v0.1.2-3 basespace
"""

FOLDER_DESCRIPTION = """\
Map FASTQ files to references for a MiSeq run folder.

This will look for MiSeq output in the path specified by the optional 
"run_folder" argument (by default, this is ".").  The output will be 
written to the path specified by the option "results_folder" argument (by 
default, "micall-results"), which itself will be taken as relative to the path 
specified by "run_folder" if it is not an absolute path.

As this will typically be run in a Dockerized setting, you should make sure
your bind mounts are appropriately set; also note that any input/output paths
will be handled as if they are *inside* the container.  The working directory
in the container will be "/data", so that the default "run_folder" argument
will resolve as "/data"; this is typically where you will mount your folder
containing MiSeq output.

For example: if you built your container as "micall:v0.1.2-3", your MiSeq data is
on your host machine as "/path/on/host/", and you want the outputs to be  
written to "/path/on/host/Results", you would invoke it as follows:

    docker run \\ 
        --mount type=bind,source=/path/on/host/,destination=/data \\ 
        micall:v0.1.2-3 folder
        
If you wanted the outputs to be written to a path that is *not* a subfolder of
the run directory, you can achieve this with an additional bind mount, and 
specifying both of the aforementioned optional parameters.  For example,
if you want to write the outputs to "/path/for/output/on/host":

    docker run \\ 
        --mount type=bind,source=/path/on/host/,destination=/data \\ 
        --mount type=bind,source=/output/path/on/host,destination=/results \\ 
        micall:v0.1.2-3 folder /data /results
"""

SAMPLE_DESCRIPTION = """\
Map FASTQ files to references for a single sample.

This will process the sample specified by the command-line arguments.  Input 
paths specified by these arguments will be taken as relative to the "base" 
path specified by "--run_folder" (default "/data"), or can be specified as 
absolute paths themselves.  The optional "results_folder" path
(default "micall-results"), will be taken as relative to the path specified by
"--run_folder" if it is not an absolute path.

As this will typically be run in a Dockerized setting, you should make sure
your bind mounts are appropriately set; also note that any input/output paths
will be handled as if they are *inside* the container.

For example: if you built your container as "micall:v0.1.2-3", your MiSeq data is
on your host machine as "/path/on/host/" as "1234A_forward.fastq" and 
"1234A_reverse.fastq", and you want the outputs to be written to 
"/path/on/host/Results", you would invoke it as follows:

    docker run \\ 
        --mount type=bind,source=/path/on/host/,destination=/data \\ 
        micall:v0.1.2-3 sample 1234A_forward.fastq 1234A_reverse.fastq

If you wanted the outputs to be written to a path that is *not* a subfolder of
the run directory, you can achieve this with an additional bind mount, and 
specifying both of the aforementioned optional parameters.  For example,
if you want to write the outputs to "/host/output/path/":

    docker run \\
        --mount type=bind,source=/path/on/host/,destination=/data \\
        --mount type=bind,source=/host/output/path/,destination=/results \\
        micall:v0.1.2-3 folder 1234A_forward.fastq 1234A_reverse.fastq \\
        /data /results
"""

HCV_SAMPLE_DESCRIPTION = """\
Map FASTQ files to references for a single HCV sample.

This will process the HCV sample specified by the command-line arguments.
Input paths specified by these arguments will be taken as relative to the 
"base" path specified by "--run_folder" (default "/data"), or can be specified 
as absolute paths themselves.  The optional "results_folder" path
(default "micall-results"), will be taken as relative to the path specified by
"--run_folder" if it is not an absolute path.

As this will typically be run in a Dockerized setting, you should make sure
your bind mounts are appropriately set; also note that any input/output paths
will be handled as if they are *inside* the container.

For example: if you built your container as "micall:v0.1.2-3", your MiSeq data is
on your host machine as "/path/on/host/" as "1234A_forward.fastq", 
"1234A_reverse.fastq", "1234A_MIDI_forward.fastq", and 
"1234A_MIDI_reverse.fastq", and you want the outputs to be written to 
"/path/on/host/micall-results", you would invoke it as follows:

    docker run \\ 
        --mount type=bind,source=/path/on/host/,destination=/data \\ 
        micall:v0.1.2-3 hcv_sample 1234A_forward.fastq 1234A_reverse.fastq \\ 
        1234A_MIDI_forward.fastq 1234A_MIDI_reverse.fastq

If you wanted the outputs to be written to a path that is *not* a subfolder of
the run directory, you can achieve this with an additional bind mount, and 
specifying both of the aforementioned optional parameters.  For example,
if you want to write the outputs to "/host/output/path/":

    docker run \\ 
        --mount type=bind,source=/path/on/host/,destination=/data \\ 
        --mount type=bind,source=/host/output/path/,destination=/results \\
        micall:v0.1.2-3 folder 1234A_forward.fastq 1234A_reverse.fastq \\ 
        1234A_MIDI_forward.fastq 1234A_MIDI_reverse.fastq \\ 
        /data /results
"""


def get_parser(default_max_active):
    parser = ArgumentParser(
        description=MAIN_DESCRIPTION,
        formatter_class=MiCallFormatter)
    subparsers = parser.add_subparsers(
        title="Sub-commands (i.e. modes of operation)",
    )
    default_folder = "micall-results-" + load_git_version()

    commands = [add_folder_parser(subparsers, default_max_active, default_folder),
                add_sample_parser(subparsers, default_folder),
                add_hcv_sample_parser(subparsers,
                                      default_max_active,
                                      default_folder),
                add_basespace_parser(subparsers, default_max_active)]
    for command_parser in commands:
        command_parser.add_argument(
            "--all_projects",
            "-a",
            action="store_true",
            help="Don't exclude any projects or seeds.")
        command_parser.add_argument(
            "--debug_remap",
            "-d",
            action="store_true",
            help="Write debug files for remapping steps.")
        command_parser.add_argument(
            "--denovo",
            action="store_true",
            help="Use de novo assembly instead of mapping to reference sequences.")
        command_parser.add_argument(
            '--skip',
            nargs='*',
            choices=TrimSteps.all,
            help='Skip over some steps to speed up the analysis.')
        command_parser.add_argument(
            "--keep_scratch",
            "-s",
            action='store_true',
            help="Don't delete the scratch folder when the run is complete.")
        command_parser.add_argument(
            "--project_code",
            "-p",
            help="Select primers to trim: HCV, HIVB, HIVGHA, or SARSCOV2.")

    return parser


def get_available_memory():
    with open('/proc/meminfo') as meminfo:
        for line in meminfo:
            label, value, *units = line.split()
            if label == 'MemAvailable:':
                assert units == ['kB'], units
                return int(value) * 1024
    raise ValueError('MemAvailable not found in /proc/meminfo.')


def add_basespace_parser(subparsers, default_max_active):
    # ####
    # BaseSpace mode.
    # ####
    basespace_parser = subparsers.add_parser(
        "basespace",
        description=BASESPACE_DESCRIPTION,
        help="Used by BaseSpace; if invoking manually you will typically not use this.",
        formatter_class=MiCallFormatter,
    )
    basespace_parser.add_argument('run_folder',
                                  default='/data',
                                  nargs='?',
                                  help='data folder filled in by BaseSpace')
    basespace_parser.add_argument(
        "--max_active",
        "-m",
        type=int,
        default=default_max_active,
        help="Maximum number of samples to process at once.",
    )
    basespace_parser.set_defaults(func=basespace_run)
    return basespace_parser


def add_folder_parser(subparsers, default_max_max_active, default_results_folder):
    # ####
    # The "process a full directory" subcommand.
    # ####
    folder_parser = subparsers.add_parser(
        "folder",
        description=FOLDER_DESCRIPTION,
        help="Process an entire run",
        formatter_class=MiCallFormatter,
    )
    folder_parser.add_argument(
        "run_folder",
        nargs="?",
        help="Path to the folder containing the MiSeq data",
        default=".",
    )
    folder_parser.add_argument(
        "results_folder",
        nargs="?",
        help="Directory to write outputs to (if Dockerized, this is "
             "the path *inside* the container; point a bind mount here).  If "
             "this path is not absolute, it will be taken as relative to "
             "--run_folder.",
        default=default_results_folder,
    )
    folder_parser.add_argument(
        "--max_active",
        "-m",
        type=int,
        default=default_max_max_active,
        help="Maximum number of samples to process at once.",
    )
    folder_parser.add_argument(
        "--fastq1s",
        nargs="*",
        help="Forward-read files to be paired with their reverse counterparts,"
             "specified with the --fastq2s option.  If these are not specified, "
             "all files with a .fastq or .fastq.gz extension will be paired "
             "by sample sheet names, or alphabetically."
    )
    folder_parser.add_argument(
        "--fastq2s",
        nargs="*",
        help="Reverse-read counterparts to the files specified by --fastq1s; "
             "ignored if --fastq1s is not specified."
    )
    folder_parser.set_defaults(func=process_folder)
    return folder_parser


def add_sample_parser(subparsers, default_results_folder):
    # ####
    # The "process a single sample" subcommand.
    # ####
    # First, the inputs.
    single_sample_parser = subparsers.add_parser(
        "sample",
        description=SAMPLE_DESCRIPTION,
        help="Process a single sample",
        formatter_class=MiCallFormatter,
    )
    single_sample_parser.add_argument(
        "fastq1",
        help="FASTQ file containing forward reads (either an absolute path or relative"
             " to the bind-mounted input directory)",
    )
    single_sample_parser.add_argument(
        "fastq2",
        help="FASTQ file containing reverse reads (either an absolute path or relative"
             " to the bind-mounted input directory)",
    )
    single_sample_parser.add_argument(
        "--bad_cycles_csv",
        help="list of tiles and cycles rejected for poor quality (either an absolute "
             "path or relative to the bind-mounted input directory)",
    )

    single_sample_parser.add_argument(
        "--run_folder",
        help="Directory to look for input files in (if Dockerized, this is the "
             "path *inside* the container; point a bind mount here).  If input "
             "file paths are not specified as absolute paths, they will be taken "
             "as being relative to this path.",
        default="/data"
    )
    # Next, the outputs.
    single_sample_parser.add_argument(
        "results_folder",
        nargs="?",
        help="(Optional) directory to write outputs to (if Dockerized, this is "
             "the path *inside* the container; point a bind mount here).  If "
             "this path is not absolute, it will be taken as relative to "
             "--run_folder.",
        default=default_results_folder,
    )
    single_sample_parser.set_defaults(func=single_sample,
                                      max_active=1)
    return single_sample_parser


def add_hcv_sample_parser(subparsers, default_max_active, default_results_folder):
    # ####
    # The "process a single HCV sample" subcommand.
    # ####
    # First, the inputs.
    hcv_sample_parser = subparsers.add_parser(
        "hcv_sample",
        description=HCV_SAMPLE_DESCRIPTION,
        help="Process a single HCV sample",
        formatter_class=MiCallFormatter,
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
        "midi_fastq1",
        help="FASTQ file containing HCV MIDI forward reads (either an absolute path "
             "or relative to the bind-mounted input directory)"
    )
    hcv_sample_parser.add_argument(
        "midi_fastq2",
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

    # Optionally customize your input and output directories.
    hcv_sample_parser.add_argument(
        "--run_folder",
        help="Directory to look for input files in (if Dockerized, this is the "
             "path *inside* the container; point a bind mount here).  If input "
             "file paths are not specified as absolute paths, they will be taken "
             "as being relative to this path.",
        default="/data"
    )
    hcv_sample_parser.add_argument(
        "results_folder",
        nargs="?",
        help="(Optional) directory to write outputs to (if Dockerized, this is "
             "the path *inside* the container; point a bind mount here).  If "
             "this path is not absolute, it will be taken as relative to "
             "--run_folder.",
        default=default_results_folder,
    )
    hcv_sample_parser.set_defaults(func=hcv_sample,
                                   max_active=min(default_max_active, 2))
    return hcv_sample_parser


def basespace_run(args):
    resolved_args = MiCallArgs(args)
    run_info = load_samples(resolved_args.run_folder)
    process_run(run_info, args)
    zip_folder(run_info.output_path, 'resistance_reports')
    zip_folder(run_info.output_path, 'coverage_maps')


def process_folder(args):
    resolved_args = MiCallArgs(args)
    run_info = link_samples(
        resolved_args.run_folder,
        resolved_args.results_folder,
        args.denovo,
        args.fastq1s,
        args.fastq2s,
        project_code=args.project_code,
        skip_censor=TrimSteps.censor in args.skip)
    process_run(run_info, args)


def process_run(run_info, args):
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

    def runner(func, inputs):
        inputs = list(inputs)
        if args.max_active > 1 and len(inputs) > 1:
            with ProcessPoolExecutor(max_workers=args.max_active) as pool:
                list(pool.map(func, inputs))
        else:
            list(map(func, inputs))

    runner(functools.partial(process_sample,
                             args=args,
                             pssm=pssm,
                             use_denovo=run_info.is_denovo),
           run_info.get_all_samples())

    runner(functools.partial(process_resistance, run_info=run_info),
           run_info.sample_groups)

    collate_samples(run_info)
    if run_summary is not None:
        summarize_samples(run_info, run_summary)
    if not args.keep_scratch:
        shutil.rmtree(run_info.scratch_path, ignore_errors=True)
    logger.info('Done.')


def single_sample(args):
    resolved_args = MiCallArgs(args)
    scratch_path = os.path.join(resolved_args.results_folder, "scratch")
    makedirs(scratch_path)

    sample_groups = []
    run_info = RunInfo(sample_groups,
                       reports=['PR_RT', 'IN', 'NS3', 'NS5a', 'NS5b'],
                       output_path=resolved_args.results_folder,
                       scratch_path=scratch_path,
                       is_denovo=args.denovo)
    sample = Sample(fastq1=resolved_args.fastq1,
                    fastq2=resolved_args.fastq2,
                    bad_cycles_csv=resolved_args.bad_cycles_csv,
                    scratch_path=scratch_path)
    sample.project_code = args.project_code
    sample_group = SampleGroup(sample)
    sample_groups.append(sample_group)

    process_run(run_info, args)


def hcv_sample(args):
    resolved_args = MiCallArgs(args)
    midi_args = MiCallArgs(args, map_midi=True)
    scratch_path = os.path.join(args.results_folder, "scratch")
    midi_scratch_path = os.path.join(args.results_folder, "scratch_midi")
    makedirs(scratch_path)
    shutil.rmtree(midi_scratch_path, ignore_errors=True)

    sample_groups = []
    run_info = RunInfo(sample_groups,
                       reports=['PR_RT', 'IN', 'NS3', 'NS5a', 'NS5b'],
                       output_path=args.results_folder,
                       scratch_path=scratch_path,
                       is_denovo=args.denovo)
    main_sample = Sample(fastq1=resolved_args.fastq1,
                         fastq2=resolved_args.fastq2,
                         bad_cycles_csv=resolved_args.bad_cycles_csv,
                         scratch_path=scratch_path)
    midi_sample = Sample(fastq1=midi_args.fastq1,
                         fastq2=midi_args.fastq2,
                         bad_cycles_csv=resolved_args.bad_cycles_csv,
                         scratch_path=midi_scratch_path)
    main_and_midi = SampleGroup(main_sample, midi_sample)
    sample_groups.append(main_and_midi)

    process_run(run_info, args)


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
        primer_node = arg_map.get('Input.project_code')
        if primer_node is None:
            project_code = None
        else:
            project_code = primer_node['Content']

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
                                                   scratch_path,
                                                   project_code),
                                       load_sample(midi_sample_json,
                                                   data_path,
                                                   scratch_path,
                                                   project_code))
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


def load_sample(sample_json,
                data_path,
                scratch_path,
                project_code: typing.Optional[str]):
    if sample_json is None:
        return None
    # Use list_fastq_files to find R1 files in BaseCalls, Alignment_*/*/Fastq, or sample folder
    sample_base = os.path.join(data_path, 'input', 'samples', sample_json['Id'])
    fastq_files = list_fastq_files(sample_base, '*_R1_*', fallback_to_run_path=True)

    if not fastq_files:
        raise RuntimeError(
            'No R1 file found for sample id {}.'.format(sample_json['Id']))

    fastq1 = str(fastq_files[0])
    sample = Sample(fastq1=fastq1,
                    basespace_id=sample_json['Id'],
                    basespace_href=sample_json['Href'])
    sample.scratch_path = os.path.join(scratch_path, sample.name)
    sample.project_code = project_code
    return sample


def link_samples(
        run_path: str,
        output_path: str,
        is_denovo: bool,
        fastq1s: typing.Sequence[str] = None,
        fastq2s: typing.Sequence[str] = None,
        project_code: str = None,
        skip_censor=False):
    """ Load the data from a run folder. """

    shutil.rmtree(output_path, ignore_errors=True)
    makedirs(output_path)

    scratch_path = os.path.join(output_path, 'scratch')
    makedirs(scratch_path)

    sample_groups = []
    run_info_path = os.path.join(run_path, 'RunInfo.xml')
    interop_path = os.path.join(run_path, 'InterOp')
    if skip_censor:
        read_sizes = None
    elif not os.path.exists(run_info_path):
        raise FileNotFoundError(
            f'Cannot censor without {run_info_path}, use "--skip trim.censor".')
    else:
        read_sizes = parse_read_sizes(run_info_path)
    run_info = RunInfo(sample_groups,
                       reports=['PR_RT', 'IN', 'NS3', 'NS5a', 'NS5b'],
                       interop_path=interop_path,
                       scratch_path=scratch_path,
                       output_path=output_path,
                       read_sizes=read_sizes,
                       is_denovo=is_denovo)

    sample_sheet_path = os.path.join(run_path, "SampleSheet.csv")
    if (fastq1s is not None and len(fastq1s) > 0
            or not os.path.exists(sample_sheet_path)):
        if fastq1s is not None and len(fastq1s) > 0:  # forward files are specified
            if fastq2s is None:
                raise ValueError("Reverse read files must also be specified.")
            elif len(fastq2s) != len(fastq1s):
                raise ValueError(
                    "The same number of forward and reverse read files must be "
                    "specified."
                )
            forward_reverse_pairs = zip(fastq1s, fastq2s)

        else:  # there is no sample sheet
            # Sort the FASTQ files alphabetically and run them in pairs.
            logger.info(
                "No sample sheet found; running on all FASTQ files in folder {}".format(
                    run_path
                )
            )
            fastq_files = (list(glob(os.path.join(run_path, "*.fastq")))
                           + list(glob(os.path.join(run_path, "*.fastq.gz"))))
            fastq_files.sort()
            forward_reverse_pairs = []
            for idx in range(0, len(fastq_files), 2):
                forward = fastq_files[idx]
                if idx == len(fastq_files) - 1:
                    # We have an odd number of FASTQ files; ignore this last one.
                    logger.info(
                        "File {} appears extraneous; omitting.".format(forward)
                    )
                    break
                reverse = fastq_files[idx + 1]
                logger.info(
                    "Pairing files {} and {}.".format(forward, reverse)
                )
                forward_reverse_pairs.append((forward, reverse))

        for forward, reverse in forward_reverse_pairs:
            sample = Sample(
                fastq1=os.path.join(run_path, forward),
                fastq2=os.path.join(run_path, reverse),
            )
            sample.project_code = project_code
            sample_groups.append(SampleGroup(sample, midi_sample=None))

    else:  # a sample sheet is specified
        # Find FASTQ files in BaseCalls, Alignment_*/*/Fastq, or run_path
        fastq_files = [str(f) for f in list_fastq_files(run_path, '*_R1_*', fallback_to_run_path=True)]
        if not fastq_files:
            file_names = []
            source_folder = None
        else:
            source_folder = os.path.dirname(fastq_files[0])
            file_names = [os.path.basename(fastq_file) for fastq_file in fastq_files]
        
        groups = find_groups(file_names, sample_sheet_path)
        for group in groups:
            main_file, midi_file = group.names
            if main_file.startswith('Undetermined'):
                continue
            main_sample = Sample(fastq1=os.path.join(str(source_folder), main_file))
            main_sample.project_code = project_code
            if midi_file is None:
                midi_sample = None
            else:
                midi_sample = Sample(fastq1=os.path.join(str(source_folder), midi_file))
                midi_sample.project_code = project_code
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
    sample.skip = args.skip
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
        read_lengths = ReadLengths4(
            forward_read=run_info.read_sizes.read1,
            index1=run_info.read_sizes.index1,
            index2=run_info.read_sizes.index2,
            reverse_read=run_info.read_sizes.read2,
        )
        phix_path = os.path.join(run_info.interop_path, 'ErrorMetricsOut.bin')
        if not os.path.exists(phix_path):
            raise FileNotFoundError(
                f'Cannot censor without {phix_path}, use "--skip trim.censor".')

        # InterOpReader expects the run folder (parent of InterOp), not the InterOp folder itself
        run_path = Path(run_info.interop_path).parent
        reader = InterOpReader(run_path)
        records = reader.read_generic_records(MetricFile.ERROR_METRICS)
        with open(run_info.quality_csv, 'w') as quality:
            write_phix_csv(quality, records, read_lengths, summary)
        with open(run_info.quality_csv) as quality, \
                open(run_info.bad_cycles_csv, 'w') as bad_cycles, \
                open(run_info.bad_tiles_csv, 'w') as bad_tiles:
            report_bad_cycles(quality, bad_cycles, bad_tiles)

        quality_metrics_path = Path(run_info.interop_path) / 'QMetricsOut.bin'
        summarize_quality(quality_metrics_path, summary, read_lengths)

        tile_metrics_path = Path(run_info.interop_path) / 'TileMetricsOut.bin'
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
        dirpath: str
        for dirpath, dirnames, filenames in os.walk(source_path):
            relpath = os.path.relpath(dirpath, source_path)
            filename: str
            for filename in filenames:
                zip_file.write(os.path.join(dirpath, filename),
                               os.path.join(folder_name, relpath, filename))


def collate_samples(run_info: RunInfo):
    """ Combine all the sample files into run files.

    :param run_info: details of the run and samples
    """
    filenames = ['remap_counts.csv',
                 'remap_conseq.csv',
                 'insertions.csv',
                 'failed_read.csv',
                 'nuc.csv',
                 'amino.csv',
                 'conseq.csv',
                 'conseq_all.csv',
                 'failed_align.csv',
                 'coverage_scores.csv',
                 'g2p.csv',
                 'g2p_summary.csv',
                 'resistance.csv',
                 'mutations.csv',
                 'nuc_mutations.csv',
                 'resistance_fail.csv',
                 'resistance_consensus.csv',
                 'cascade.csv',
                 'merge_lengths.csv',
                 'concordance.csv',
                 'concordance_seed.csv']
    if run_info.is_denovo:
        filenames += ['conseq_stitched.csv', 'conseq_region.csv',
                      'unstitched_cascade.csv', 'unstitched_conseq.csv', 'unstitched_contigs.csv']
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
        if os.path.exists(sample_info.genome_concordance_svg):
            safe_file_move(sample_info.genome_concordance_svg,
                           os.path.join(genome_coverage_path,
                                        sample_info.name + '_genome_concordance.svg'))
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


# noinspection PyTypeChecker,PyUnresolvedReferences
def main():
    available_memory = get_available_memory()
    recommended_memory = 1 << 30  # 1GB
    default_max_active = max(1, available_memory // recommended_memory)
    cpu_count = len(os.sched_getaffinity(0))
    default_max_active = min(default_max_active, cpu_count)

    parser = get_parser(default_max_active)
    args = parser.parse_args()
    if not hasattr(args, "func"):
        # No valid subcommand was given; print the help message and exit.
        parser.print_help()
        return
    args.skip = args.skip or ()
    if available_memory < recommended_memory:
        logger.warning("Available memory is less than 1GB, processing may fail.")
    logger.info("Starting on %s with %d CPU's.",
                socket.gethostname(),
                cpu_count)
    args.func(args)


if __name__ == '__main__':
    main()
