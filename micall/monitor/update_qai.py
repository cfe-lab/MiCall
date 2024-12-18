# Script to update QAI with information from the
# conseq.csv files produced by the MiSeq pipeline.
# To execute as a script, run python -m micall.monitor.update_qai

import csv
import typing
from argparse import SUPPRESS
from collections import defaultdict
from datetime import datetime
import logging
from functools import partial
from pathlib import Path

from micall.monitor.sample_watcher import PipelineType
from operator import itemgetter, getitem
import os

from micall.monitor import qai_helper
from micall.utils import sample_sheet_parser
from micall.core.project_config import ProjectConfig, G2P_SEED_NAME
from .kive_watcher import wait_for_retry

logger = logging.getLogger('update_qai')


def parse_args():
    import argparse
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
    args = parser.parse_args()
    if not hasattr(args, 'qai_password'):
        args.qai_password = os.environ.get('MICALL_QAI_PASSWORD', 'testing')
    return args


def build_conseqs(conseqs_file, run_name, sample_sheet):
    """
    Parses a Pipeline-produced conseq file and builds JSON objects to send
    to QAI.

    @param conseqs_file: An open file that contains the consensus sequences
        from the counts2csf step for all samples in the run.
    @param run_name: date and machine number for the run
    @param sample_sheet: The data parsed from the sample sheet.
    @return an array of JSON hashes, one for each conseq.
    """

    result = []
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
        sample_info = ss["Data"][fastq_filename]
        orig_sample_name = sample_info["orig_sample_name"]
        # FIXME if row["sequence"] is blank we replace it with a dash.
        # Need Conan to make that row blank-able.
        curr_seq = row["sequence"] if len(row["sequence"]) > 0 else "-"
        seeds = sample_seeds.get(fastq_filename)
        if seeds is None:
            seeds = set()
            for sample_row in ss["DataSplit"]:
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


def build_review_decisions(coverage_file, collated_counts_file, cascade_file,
                           sample_sheet, sequencings, project_regions,
                           regions):
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

    project_region_map = dict([((entry['project_name'],
                                 entry['coordinate_region_name']), entry['id'])
                               for entry in project_regions])
    region_map = dict([(entry['name'], entry['id']) for entry in regions])
    # noinspection PyTypeChecker
    sample_tags = dict(map(itemgetter('filename', 'tags'), sample_sheet['DataSplit']))
    # noinspection PyTypeChecker
    sample_names = dict(map(itemgetter('tags', 'filename'), sample_sheet['DataSplit']))

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

    sequencing_map = defaultdict(dict)  # {tags: {project: sequencing}}
    for sequencing in sequencings:
        sequencing_map[sequencing['tag']][
            sequencing['target_project']] = sequencing

    targeted_projects = set(map(itemgetter('target_project'), sequencings))

    decisions = {}  # {(sample_name, region): decision}
    # sample,project,region,q.cut,min.coverage,which.key.pos,off.score,on.score
    for coverage in csv.DictReader(coverage_file):
        tags = sample_tags[coverage['sample']]
        project_map = sequencing_map.get(tags)
        if project_map is None:
            raise KeyError("No sequencing found with tags '%s' for %s. Are "
                           "tagged layouts missing?" % (tags, coverage_file.name))
        sequencing = project_map.get(coverage['project'])
        if sequencing is not None:
            score = read_int(coverage, 'on.score')
        else:
            score = read_int(coverage, 'off.score')
            first_project = sorted(project_map.keys())[0]
            sequencing = project_map[first_project]
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


def upload_review_to_qai(coverage_file, collated_counts_file, cascade_file,
                         run, sample_sheet, conseqs, session,
                         pipeline_version):
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

    project_regions = session.get_json("/lab_miseq_project_regions?pipeline=" +
                                       pipeline_version)
    if not project_regions:
        raise RuntimeError('Unknown pipeline: ' + pipeline_version)

    regions = session.get_json("/lab_miseq_regions")

    decisions = build_review_decisions(coverage_file, collated_counts_file,
                                       cascade_file, sample_sheet, sequencings,
                                       project_regions, regions)

    session.post_json(
        "/lab_miseq_reviews", {
            'runid': runid,
            'pipeline_id': find_pipeline_id(session, pipeline_version),
            'lab_miseq_review_decisions': decisions,
            'lab_miseq_conseqs': conseqs
        })


def upload_proviral_tables(session, result_folder, run):
    proviral_file = os.path.join(result_folder, "proviral",
                                 "table_precursor.csv")
    aln_proviral_file = os.path.join(result_folder, "proviral",
                                     "aligned_table_precursor.csv")
    upload_proviral_csv(session, run['id'], proviral_file, '/proviral/create')
    logger.info('proviral upload success!')
    upload_proviral_csv(session, run['id'],
                        aln_proviral_file,
                        '/aligned_proviral/create')
    logger.info('aligned proviral upload success!')


def upload_proviral_csv(session, run_id, csv_path, endpoint):
    responses = []
    with open(csv_path, newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            row['lab_miseq_run_id'] = run_id
            response = session.post_json(endpoint, row)
            responses.append(response)
    return responses


def clean_runname(runname):
    try:
        rundate = datetime.strptime(runname, '%d-%b-%y')
        cleaned_runname = datetime.strftime(rundate, '%d-%b-%Y')
    except ValueError:
        cleaned_runname = runname
    return cleaned_runname


def find_run(session, runname):
    """ Query QAI to find the run id for a given run name.

    @return: a hash with the attributes of the run record, including a
        sequencing summary of all the samples and their target projects.
    """
    cleaned_runname = clean_runname(runname)
    runs = session.get_json("/lab_miseq_runs?summary=sequencing&runname=" +
                            cleaned_runname)
    rowcount = len(runs)
    if rowcount == 0:
        raise RuntimeError(
            "No run found with runname {!r}.".format(cleaned_runname))
    if rowcount != 1:
        raise RuntimeError("Found {} runs with runname {!r}.".format(
            rowcount, cleaned_runname))
    return runs[0]


def find_pipeline_id(session, pipeline_version):
    """ Query QAI to find the pipeline id for the current version.

    :param session: open session on QAI
    :param str pipeline_version: 'X.Y' description of pipeline version to find
    @return: the pipeline id.
    """
    pipelines = session.get_json("/lab_miseq_pipelines?version=" +
                                 pipeline_version)
    rowcount = len(pipelines)
    if rowcount == 0:
        raise RuntimeError(
            "No pipeline found with version {!r}.".format(pipeline_version))
    if rowcount != 1:
        raise RuntimeError("Found {} pipelines with version {!r}.".format(
            rowcount, pipeline_version))
    return pipelines[0]['id']


def load_ok_sample_regions(result_folder):
    ok_sample_regions = set()
    coverage_file = os.path.join(result_folder, 'coverage_scores.csv')
    with open(coverage_file, "r") as f:
        reader = csv.DictReader(f)
        for row in reader:
            if row['on.score'] == '4':
                ok_sample_regions.add(
                    (row['sample'], row['region'], row['q.cut']))

    return ok_sample_regions


def process_folder(item, qai_server, qai_user, qai_password, pipeline_version):
    result_folder, pipeline_group = item

    run_path = Path(result_folder).parent.parent
    sample_sheet = sample_sheet_parser.read_sample_sheet_and_overrides(
        run_path / 'SampleSheet.csv')

    attempt_count = 0
    with qai_helper.Session() as session:
        # noinspection PyBroadException
        try:
            session.login(qai_server, qai_user, qai_password)
            run = find_run(session, sample_sheet["Experiment Name"])

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
                process_remapped(result_folder, session, run, pipeline_version)

        except Exception:
            attempt_count += 1
            wait_for_retry(attempt_count)


def process_remapped(result_folder, session, run, pipeline_version):
    logger.info('Uploading data to Oracle from {}'.format(result_folder))
    collated_conseqs = os.path.join(result_folder, 'conseq.csv')
    collated_counts = os.path.join(result_folder, 'remap_counts.csv')
    cascade = os.path.join(result_folder, 'cascade.csv')
    coverage_scores = os.path.join(result_folder, 'coverage_scores.csv')
    run_path = Path(result_folder).parent.parent
    sample_sheet = sample_sheet_parser.read_sample_sheet_and_overrides(
        run_path / 'SampleSheet.csv')

    with open(collated_conseqs) as f:
        conseqs = build_conseqs(f, run['runname'], sample_sheet)

    with open(coverage_scores) as f, \
            open(collated_counts) as f2, \
            open(cascade) as f3:
        upload_review_to_qai(f, f2, f3, run, sample_sheet, conseqs, session,
                             pipeline_version)
        logger.info('Remapped upload success!')


def upload_loop(qai_server, qai_user, qai_password, pipeline_version,
                upload_queue):
    # noinspection PyBroadException
    try:
        with qai_helper.Session() as session:
            # Try logging in to QAI, just so we learn about problems at launch.
            session.login(qai_server, qai_user, qai_password)
    except Exception:
        logger.error('Unable to log in to QAI.', exc_info=True)

    while True:
        item = upload_queue.get()
        if item is None:
            break
        process_folder(item, qai_server, qai_user, qai_password,
                       pipeline_version)


def main():
    args = parse_args()

    process_folder((args.result_folder, args.pipeline_group),
                   args.qai_server,
                   args.qai_user,
                   args.qai_password,
                   args.pipeline_version)

    logger.info('Completed upload to Oracle.')


if __name__ == "__main__":
    main()
