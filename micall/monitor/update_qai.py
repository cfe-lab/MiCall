# Script to update QAI with information from the
# conseq.csv files produced by the MiSeq pipeline.
# To execute as a script, run python -m micall.monitor.update_qai
#
# RETRY IMPLEMENTATION:
# This script implements robust retry logic for both network and disk operations:
# - Network failures: QAI login, run lookups, data uploads
# - Disk failures: File reads, sample sheet parsing, flag writing
# - Uses exponential backoff with the same retry parameters as kive_watcher
# - Retries up to MAX_RETRY_ATTEMPTS (10) times by default
# - Logs warnings for each retry attempt, errors after exhausting retries
# - Always creates done_qai_upload flag, even when coverage_scores.csv is missing

import csv
import typing
from argparse import SUPPRESS
from collections import defaultdict
from datetime import datetime
import logging
from functools import partial
from pathlib import Path
from typing import Dict, Sequence, Tuple, List, Set, Optional, TextIO, TypedDict
import argparse
from queue import Queue
from urllib.parse import urlparse

from micall.monitor.sample_watcher import PipelineType
from operator import itemgetter, getitem
import os

from micall.monitor import qai_helper
from micall.utils import sample_sheet_parser
from micall.core.project_config import ProjectConfig, G2P_SEED_NAME
from micall.monitor import disk_operations

logger = logging.getLogger('update_qai')

# Use the same retry configuration as kive_watcher
MAX_RETRY_ATTEMPTS = 100 


def retry_operation(operation, operation_name: str):
    """Execute an operation with retry logic for both network and disk failures.
    
    This function handles:
    - Network failures (requests exceptions, connection errors)
    - Disk failures (IO errors, permission errors)
    - QAI API failures (authentication, server errors)
    - File system issues (missing files, corrupted data)
    
    Args:
        operation: Callable that performs the operation
        operation_name: Human-readable name for logging
        max_attempts: Maximum number of retry attempts
        
    Returns:
        The result of the successful operation
        
    Raises:
        RuntimeError: If all retry attempts fail
    """

    from micall.monitor.kive_watcher import wait_for_retry

    attempt_count = 0
    start_time = None
    last_exception = None

    while attempt_count < MAX_RETRY_ATTEMPTS:
        try:
            return operation()
        except Exception as ex:
            attempt_count += 1
            last_exception = ex
            
            # Record start time on first failure
            if start_time is None:
                start_time = datetime.now()
            
            if attempt_count >= MAX_RETRY_ATTEMPTS:
                logger.error(f'{operation_name} failed after {MAX_RETRY_ATTEMPTS} attempts', exc_info=True)
                break
            
            logger.warning(f'{operation_name} failed (attempt {attempt_count}/{MAX_RETRY_ATTEMPTS})', exc_info=True)
            wait_for_retry(attempt_count, start_time)
    
    # If we get here, all attempts failed
    raise RuntimeError(f'{operation_name} failed after {MAX_RETRY_ATTEMPTS} attempts') from last_exception


Session = qai_helper.Session


class Sequencing(TypedDict):
    id: str
    tag: str
    target_project: str


class Run(TypedDict):
    runname: str
    sequencing_summary: Sequence[Sequencing]
    id: str


class ProjectRegion(TypedDict):
    id: int
    project_name: str
    seed_region_names: Sequence[str]
    coordinate_region_name: str


def retry_qai_login(session: Session, qai_server: str, qai_user: str, qai_password: str, context: str = "") -> None:
    """Login to QAI with retry logic."""
    retry_operation(
        lambda: session.login(qai_server, qai_user, qai_password),
        f"QAI login{' for ' + context if context else ''}"
    )


def retry_find_run(session: Session, experiment_name: str) -> Run:
    """Find QAI run with retry logic."""
    runs = retry_operation(
        lambda: find_run(session, experiment_name),
        f"Finding QAI run for {experiment_name}"
    )

    rowcount: int = len(runs)
    if rowcount == 0:
        raise RuntimeError(
            "No run found with runname {!r}.".format(experiment_name))
    if rowcount != 1:
        raise RuntimeError("Found {} runs with runname {!r}.".format(
            rowcount, experiment_name))
    return runs[0]


def parse_args() -> argparse.Namespace:
    pipeline_parser = partial(getitem, PipelineType)
    parser = argparse.ArgumentParser(
        description="Update the Oracle database with conseq information",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("result_folder",
                        help="Result folder that holds the conseq.csv file")
    parser.add_argument('--pipeline_version',
                        default='0-dev',
                        help='version suffix for batch names and folder names')
    parser.add_argument('--pipeline_group',
                        default=PipelineType.MAIN,
                        type=pipeline_parser,
                        choices=(PipelineType.MAIN,
                                 PipelineType.DENOVO_MAIN,
                                 PipelineType.PROVIRAL),
                        help='group of results to upload')
    parser.add_argument('--qai_server',
                        default=os.environ.get('MICALL_QAI_SERVER',
                                               'http://localhost:4567'),
                        help='server to post reviews on')
    parser.add_argument('--qai_user',
                        default=os.environ.get('MICALL_QAI_USER', 'bob'),
                        help='user name for QAI server')
    parser.add_argument('--qai_password',
                        default=SUPPRESS,
                        help='password for QAI server (default not shown)')

    verbosity_group = parser.add_mutually_exclusive_group()
    verbosity_group.add_argument('--verbose', action='store_true', help='Increase output verbosity.')
    verbosity_group.add_argument('--no-verbose', action='store_true', help='Normal output verbosity.', default=True)
    verbosity_group.add_argument('--debug', action='store_true', help='Maximum output verbosity.')
    verbosity_group.add_argument('--quiet', action='store_true', help='Minimize output verbosity.')

    args = parser.parse_args()
    if not hasattr(args, 'qai_password'):
        args.qai_password = os.environ.get('MICALL_QAI_PASSWORD', 'testing')
    return args


def build_conseqs(conseqs_file: TextIO, run_name: str, sample_sheet: Dict[str, object]) -> List[Dict[str, object]]:
    """
    Parses a Pipeline-produced conseq file and builds JSON objects to send
    to QAI.

    @param conseqs_file: An open file that contains the consensus sequences
        from the counts2csf step for all samples in the run.
    @param run_name: date and machine number for the run
    @param sample_sheet: The data parsed from the sample sheet.
    @return an array of JSON hashes, one for each conseq.
    """

    result: List[Dict[str, object]] = []
    ss = sample_sheet
    conseqs_csv = csv.DictReader(conseqs_file)
    # ss["Data"] is keyed by (what should be) the FASTQ
    # filename, which looks like
    #
    # [sample name with ; and _ replaced by -]_S[sample number].
    #
    # Meanwhile, entries in conseqs_file have a "sample" field holding
    # just the sample name (also with ; and _ replaced).  We make a
    # lookup table to get the FASTQ filename just from the first part.
    # This will make subsequent steps easier (avoids having to do a
    # search through a list/dict of dicts).
    # FASTQ_lookup = {}
    # filename_re = re.compile("(.+)_S.+")
    # for fastq_filename in ss["Data"]:
    #     sample_name = filename_re.match(fastq_filename).group(1)
    #     FASTQ_lookup[sample_name] = fastq_filename

    projects = ProjectConfig.loadDefault()
    sample_seeds: typing.Dict[str, typing.Set[str]] = {}

    for row in conseqs_csv:
        # Each row of this file looks like:
        # sample,region,q-cutoff,s-number,consensus-percent-cutoff,sequence
        # We want to take the "sample" entry and get the corresponding
        # original Sample_Name from the sample sheet. In version 2, this
        # looks like [sample name]~[project name]#[...]
        # In version 1, this looked like [sample name]~[project name]#[...]
        # but both ; and _ got garbled by the MiSeq instrument itself.
        # Thus we have to work around it.
        fastq_filename = row["sample"]

        data_raw = ss.get("Data")
        if data_raw is None:
            raise RuntimeError('Sample sheet is missing "Data" section')
        if hasattr(data_raw, '__getitem__'):
            data = data_raw
        else:
            raise RuntimeError('Sample sheet "Data" section is not indexable')

        data_split_raw = ss.get("DataSplit")
        if data_split_raw is None:
            raise RuntimeError('Sample sheet is missing "DataSplit" section')
        if hasattr(data_split_raw, '__getitem__') and hasattr(data_split_raw, '__iter__'):
            data_split = data_split_raw
        else:
            raise RuntimeError('Sample sheet "DataSplit" section is not indexable and iterable')

        sample_info = data[fastq_filename]
        orig_sample_name = sample_info["orig_sample_name"]
        # FIXME if row["sequence"] is blank we replace it with a dash.
        # Need Conan to make that row blank-able.
        curr_seq = row["sequence"] if len(row["sequence"]) > 0 else "-"
        seeds = sample_seeds.get(fastq_filename)
        if seeds is None:
            seeds = set()
            for sample_row in data_split:
                if sample_row['filename'] != fastq_filename:
                    continue
                target_project = sample_row['project']
                try:
                    seeds |= projects.getProjectSeeds(target_project)
                except KeyError:
                    if target_project != 'Unknown':
                        logger.warning('Failed to load project seeds for %s in %s.',
                                       fastq_filename,
                                       run_name,
                                       exc_info=True)
            sample_seeds[fastq_filename] = seeds
        ok_for_release = row["region"] in seeds
        result.append({
            "samplename": orig_sample_name,
            # July 9, 2014: we can't do this properly right now
            # without a lookup table that is yet to be fully
            # defined.
            "testcode": None,
            "conseq_cutoff": row["consensus-percent-cutoff"],
            "region": row["region"],
            "qcutoff": float(row["q-cutoff"]),
            "snum": fastq_filename.split('_')[-1],
            "seq": curr_seq,
            "ok_for_release": ok_for_release
        })
    return result


def build_review_decisions(coverage_file: TextIO, collated_counts_file: TextIO, cascade_file: TextIO,
                           sample_sheet: Dict[str, object], sequencings: Sequence[Sequencing], project_regions: Sequence[ProjectRegion],
                           regions: List[Dict[str, object]]) -> List[Dict[str, object]]:
    """ Build a list of request objects that will create the review decision
    records.

    @param coverage_file: CSV file with coverage scores
    @param collated_counts_file: CSV file with read counts
    @param cascade_file: CSV file with read counts throughout the pipeline
    @param sample_sheet: the sample sheet for the run
    @param sequencings: the sequencing records from QAI
    @param project_regions: [{"id": project_region_id,
                              "project_name": project_name,
                              "seed_region_names": [seed_region_name],
                              "coordinate_region_name": coordinate_region_name}]
    @param regions: [{"id": region_id, "name": region_name}]
    """

    data_split_raw = sample_sheet.get("DataSplit")
    if data_split_raw is None:
        raise RuntimeError('Sample sheet is missing "DataSplit" section')
    if hasattr(data_split_raw, '__iter__'):
        data_split = data_split_raw
    else:
        raise RuntimeError('Sample sheet "DataSplit" section is not iterable')

    project_region_map = dict([((entry['project_name'],
                                 entry['coordinate_region_name']), entry['id'])
                               for entry in project_regions])
    region_map = dict([(entry['name'], entry['id']) for entry in regions])
    # noinspection PyTypeChecker
    sample_tags = dict(map(itemgetter('filename', 'tags'), data_split))
    # noinspection PyTypeChecker
    sample_names = dict(map(itemgetter('tags', 'filename'), data_split))

    def read_int(table, name):
        ret = float(table[name])
        if float(int(ret)) != ret:
            raise ValueError(f"Bad value for {name!r}: {ret!r}. Expected an integer.")
        return int(ret)

    counts_map = {}  # {tags: raw, (tags, seed): mapped]}
    # sample,type,count
    for counts in csv.DictReader(collated_counts_file):
        count = read_int(counts, 'count')
        tags = sample_tags[counts['sample']]
        count_type = counts['type']
        if count_type not in ('raw', 'unmapped'):
            seed = count_type.split(' ', 1)[1]
            key = tags, seed
            counts_map[key] = count

    unreported_tags = set()
    for counts in csv.DictReader(cascade_file):
        tags = sample_tags[counts['sample']]
        counts_map[tags] = read_int(counts, 'demultiplexed') * 2
        unreported_tags.add(tags)

        key = tags, G2P_SEED_NAME
        counts_map[key] = read_int(counts, 'v3loop') * 2

    sequencing_map: Dict[str, Dict[str, Sequencing]] = defaultdict(dict)  # {tags: {project: sequencing}}
    for sequencing in sequencings:
        sequencing_map[sequencing['tag']][
            sequencing['target_project']] = sequencing

    targeted_projects = set(map(itemgetter('target_project'), sequencings))

    decisions: Dict[Tuple[str, str], Dict[str, object]] = {}  # {(sample_name, region): decision}
    # sample,project,region,q.cut,min.coverage,which.key.pos,off.score,on.score
    for coverage in csv.DictReader(coverage_file):
        tags = sample_tags[coverage['sample']]
        project_map = sequencing_map.get(tags)
        if project_map is None:
            raise KeyError("No sequencing found with tags '%s' for %s. Are "
                           "tagged layouts missing?" % (tags, coverage_file.name))
        sequencing: Optional[Dict[str, object]] = project_map.get(coverage['project'])  # type: ignore
        if sequencing is not None:
            score = read_int(coverage, 'on.score')
        else:
            score = read_int(coverage, 'off.score')
            sequencing = project_map[sorted(project_map.keys())[0]]
        project_region_id = project_region_map[(coverage['project'],
                                                coverage['region'])]
        raw_count = counts_map[tags]
        seed = coverage['seed']
        mapped_count = counts_map.get((tags, seed))
        seed_region_id = region_map[seed]

        decision_key = (coverage['sample'], coverage['region'])
        previous_decision = decisions.get(decision_key)
        is_replacement = (previous_decision is None
                          or score > previous_decision['score']
                          or (score == previous_decision['score']
                              and coverage['project'] in targeted_projects))
        if is_replacement:
            decisions[decision_key] = {
                'sequencing_id': sequencing['id'],
                'project_region_id': project_region_id,
                'seed_region_id': seed_region_id,
                'sample_name': coverage['sample'],
                'score': score,
                'min_coverage': read_int(coverage, 'min.coverage'),
                'min_coverage_pos': read_int(coverage, 'which.key.pos'),
                'raw_reads': raw_count,
                'mapped_reads': mapped_count
            }
            unreported_tags.discard(tags)
    for tags in unreported_tags:
        sample_name = sample_names[tags]
        project_map = sequencing_map.get(tags)
        if project_map is None:
            raise KeyError("No sequencing found with tags '%s'." % tags)
        first_project = sorted(project_map.keys())[0]
        sequencing = project_map[first_project]
        decision_key = sample_name
        decisions[decision_key] = {
            'sequencing_id': sequencing['id'],
            'sample_name': sample_name,
            'raw_reads': counts_map[tags],
            'mapped_reads': 0
        }
    return list(decisions.values())


def upload_review_to_qai(coverage_file: TextIO, collated_counts_file: TextIO, cascade_file: TextIO,
                         run: Run, sample_sheet: Dict[str, object], conseqs: List[Dict[str, object]], session: Session,
                         pipeline_version: str) -> None:
    """ Create a review.

    @param coverage_file: the coverage scores to upload
    @param collated_counts_file: CSV file of read counts to upload
    @param cascade_file: CSV file of read counts throughout the pipeline
    @param run: a hash with the attributes of the run record, including a
        sequencing summary of all the samples and their target projects
    @param sample_sheet: details of the run so we can tell which sample used
        which tags
    @param conseqs: an array of JSON hashes to pass to QAI for the conseq
        child records
    @param session: the QAI session
    @param str pipeline_version: 'X.Y' describing the current version
    """

    runid = run['id']
    sequencings = run['sequencing_summary']
    project_regions = retry_operation(
        lambda: session.get_json("/lab_miseq_project_regions?pipeline=" + pipeline_version),
        f"Getting project regions for pipeline {pipeline_version}"
    )
    if not project_regions:
        raise RuntimeError('Unknown pipeline: ' + pipeline_version)

    regions = retry_operation(
        lambda: session.get_json("/lab_miseq_regions"),
        "Getting lab MiSeq regions"
    )

    decisions = build_review_decisions(coverage_file, collated_counts_file,
                                       cascade_file, sample_sheet, sequencings,
                                       project_regions, regions)

    retry_operation(
        lambda: session.post_json(
            "/lab_miseq_reviews", {
                'runid': runid,
                'pipeline_id': find_pipeline_id(session, pipeline_version),
                'lab_miseq_review_decisions': decisions,
                'lab_miseq_conseqs': conseqs
            }),
        f"Posting lab MiSeq review for run {runid}"
    )


def upload_proviral_tables(session: Session, result_folder: str, run: Run) -> None:
    proviral_file: str = os.path.join(result_folder, "proviral",
                                 "table_precursor.csv")
    aln_proviral_file: str = os.path.join(result_folder, "proviral",
                                     "aligned_table_precursor.csv")
    upload_proviral_csv(session, run['id'], proviral_file, '/proviral/create')
    logger.info('proviral upload success!')
    upload_proviral_csv(session, run['id'],
                        aln_proviral_file,
                        '/aligned_proviral/create')
    logger.info('aligned proviral upload success!')


def upload_proviral_csv(session: Session, run_id: str, csv_path: str, endpoint: str) -> List[object]:
    responses: List[object] = []
    with open(csv_path, newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            row['lab_miseq_run_id'] = run_id
            response = retry_operation(
                lambda: session.post_json(endpoint, row),
                f"Posting proviral data to {endpoint}"
            )
            responses.append(response)
    return responses


def clean_runname(runname: str) -> str:
    try:
        rundate = datetime.strptime(runname, '%d-%b-%y')
        cleaned_runname = datetime.strftime(rundate, '%d-%b-%Y')
    except ValueError:
        cleaned_runname = runname
    return cleaned_runname


def find_run(session: Session, runname: str) -> Sequence[Run]:
    """ Query QAI to find the run id for a given run name.

    @return: a hash with the attributes of the run record, including a
        sequencing summary of all the samples and their target projects.
    """
    cleaned_runname = clean_runname(runname)
    runs: Sequence[Run] = retry_operation(
        lambda: session.get_json("/lab_miseq_runs?summary=sequencing&runname=" + cleaned_runname),
        f"Finding run with name {cleaned_runname}"
    )
    return runs


def find_pipeline_id(session: Session, pipeline_version: str) -> object:
    """ Query QAI to find the pipeline id for the current version.

    :param session: open session on QAI
    :param str pipeline_version: 'X.Y' description of pipeline version to find
    @return: the pipeline id.
    """
    pipelines: List[Dict[str, object]] = retry_operation(
        lambda: session.get_json("/lab_miseq_pipelines?version=" + pipeline_version),
        f"Finding pipeline with version {pipeline_version}"
    )
    rowcount: int = len(pipelines)
    if rowcount == 0:
        raise RuntimeError(
            "No pipeline found with version {!r}.".format(pipeline_version))
    if rowcount != 1:
        raise RuntimeError("Found {} pipelines with version {!r}.".format(
            rowcount, pipeline_version))
    first_pipeline = pipelines[0]
    pipeline_id = first_pipeline.get('id')
    if pipeline_id is None:
        raise RuntimeError("Pipeline record is missing 'id' field.")
    return pipeline_id


def load_ok_sample_regions(result_folder: str) -> Set[Tuple[str, str, str]]:
    ok_sample_regions: Set[Tuple[str, str, str]] = set()
    coverage_file: str = os.path.join(result_folder, 'coverage_scores.csv')
    with open(coverage_file, "r") as f:
        reader = csv.DictReader(f)
        for row in reader:
            if row['on.score'] == '4':
                ok_sample_regions.add(
                    (row['sample'], row['region'], row['q.cut']))

    return ok_sample_regions


def process_folder(item: Tuple[Path, PipelineType], qai_server: str, qai_user: str, qai_password: str, pipeline_version: str) -> None:
    result_folder: Path
    pipeline_group: PipelineType
    result_folder, pipeline_group = item

    def read_sample_sheet_with_retry():
        run_path: Path = Path(result_folder).parent.parent
        return sample_sheet_parser.read_sample_sheet_and_overrides(run_path / 'SampleSheet.csv')

    # Read sample sheet with retry for disk operations
    sample_sheet: Dict[str, object] = retry_operation(
        read_sample_sheet_with_retry,
        f"Reading sample sheet for {result_folder}"
    )

    experiment_name = sample_sheet.get("Experiment Name")
    if experiment_name is None:
        raise RuntimeError("Sample sheet is missing 'Experiment Name' field.")
    if not isinstance(experiment_name, str):
        raise RuntimeError("'Experiment Name' field is not a string.")

    with qai_helper.Session() as session:
        retry_qai_login(session, qai_server, qai_user, qai_password, f"{result_folder}")    
        run = retry_find_run(session, experiment_name)

        # PROVIRAL
        if pipeline_group == PipelineType.PROVIRAL:
            # upload_proviral_tables(session, result_folder, run)
            pass
        # DENOVO
        elif pipeline_group in (PipelineType.DENOVO_MAIN,
                                PipelineType.DENOVO_MIDI,
                                PipelineType.DENOVO_RESISTANCE):
            # Do nothing
            pass
        # REMAPPED
        else:
            retry_operation(
                lambda: process_remapped(result_folder, session, run, pipeline_version),
                f"Processing remapped data for {result_folder}"
            )


def process_remapped(result_folder: Path, session: Session, run: Run, pipeline_version: str) -> None:
    logger.info('Uploading data to Oracle from {}'.format(result_folder))
    result_folder = Path(result_folder)
    collated_conseqs: Path = result_folder / 'conseq.csv'
    collated_counts: Path = result_folder / 'remap_counts.csv'
    cascade: Path = result_folder / 'cascade.csv'
    coverage_scores: Path = result_folder / 'coverage_scores.csv'
    run_path: Path = result_folder.parent.parent

    def read_sample_sheet_with_retry():
        return sample_sheet_parser.read_sample_sheet_and_overrides(run_path / 'SampleSheet.csv')

    def read_conseqs_with_retry():
        with open(collated_conseqs) as f:
            return build_conseqs(f, run['runname'], sample_sheet)

    def upload_with_retry():
        with open(coverage_scores) as f, \
                open(collated_counts) as f2, \
                open(cascade) as f3:
            upload_review_to_qai(f, f2, f3, run, sample_sheet, conseqs, session,
                                 pipeline_version)

    # Read sample sheet with retry
    sample_sheet: Dict[str, object] = retry_operation(
        read_sample_sheet_with_retry,
        f"Reading sample sheet for {result_folder}"
    )

    # Read conseqs with retry
    conseqs: List[Dict[str, object]] = retry_operation(
        read_conseqs_with_retry,
        f"Reading conseqs file for {result_folder}"
    )

    # Upload to QAI with retry
    retry_operation(
        upload_with_retry,
        f"Uploading to QAI for {result_folder}"
    )

    logger.info('Remapped upload success!')


def upload_loop(qai_server: str, qai_user: str, qai_password: str, pipeline_version: str,
                upload_queue: Queue[Optional[Tuple[Path, PipelineType]]]) -> None:

    # Check if `qai_server` is a valid URL.
    parsed_url = urlparse(qai_server)
    if not (parsed_url.scheme and parsed_url.netloc):
        raise ValueError(f"Invalid QAI server URL: {qai_server}")

    while True:
        item = upload_queue.get()
        if item is None:
            break
        results_path, pipeline_group = item
        results_path = Path(results_path)
        done_qai_flag = results_path / 'done_qai_upload'
        coverage_scores = results_path / 'coverage_scores.csv'

        def check_coverage_file():
            return coverage_scores.exists()

        # Check coverage file with retry (in case of temporary file system issues)
        has_coverage_file = retry_operation(
            check_coverage_file,
            f"Checking coverage file for {results_path}",
        )

        # Skip upload if coverage_scores.csv is missing, but still mark as done.
        if has_coverage_file:
            try:
                process_folder((results_path, pipeline_group),
                               qai_server, qai_user, qai_password,
                               pipeline_version)
                status = 'succeeded'
            except Exception:
                logger.error("QAI upload failed for %s", results_path, exc_info=True)
                status = 'failed'
        else:
            logger.debug("Skipping QAI upload for %s (coverage_scores.csv missing)", results_path)
            status = 'skipped'

        # Mark completed to prevent re-upload with retry for disk operations
        def write_done_flag():
            disk_operations.write_text(done_qai_flag, status)

        retry_operation(write_done_flag, f"Writing done_qai_upload flag for {results_path}")


def main() -> None:
    args = parse_args()

    if args.quiet:
        logger.setLevel(logging.ERROR)
    elif args.verbose:
        logger.setLevel(logging.INFO)
    elif args.debug:
        logger.setLevel(logging.DEBUG)
    else:
        logger.setLevel(logging.WARN)

    logging.basicConfig(level=logger.level)

    process_folder((args.result_folder, args.pipeline_group),
                   args.qai_server,
                   args.qai_user,
                   args.qai_password,
                   args.pipeline_version)

    logger.info('Completed upload to Oracle.')


if __name__ == "__main__":
    main()
