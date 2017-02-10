# Script to update QAI with information from the
# conseq.csv files produced by the MiSeq pipeline.
# To execute as a script, run python -m micall.monitor.update_qai

import csv
from collections import defaultdict
from datetime import datetime
from glob import glob
import logging
from operator import itemgetter
import os

from micall import settings  # Import first for logging configuration.

from micall.monitor import qai_helper
from micall.utils import sample_sheet_parser
from micall.core.project_config import ProjectConfig

logger = logging.getLogger('update_qai')


def parse_args():
    import argparse
    parser = argparse.ArgumentParser(
        version=settings.pipeline_version,
        description="Update the Oracle database with conseq information")
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--result_folder",
                       "-r",
                       help="Result folder that holds the conseq.csv file")
    group.add_argument("--load_all",
                       "-a",
                       action="store_true",
                       help="load all folders under RAW_DATA that have results.")
    args = parser.parse_args()
    return args, parser


def build_conseqs(conseqs_file,
                  run,
                  sample_sheet,
                  ok_sample_regions):
    """
    Parses a Pipeline-produced conseq file and builds JSON objects to send
    to QAI.

    @param conseqs_file: An open file that contains the consensus sequences
        from the counts2csf step for all samples in the run.
    @param run: a hash with the attributes of the run record, including a
        sequencing summary of all the samples and their target projects
    @param sample_sheet: The data parsed from the sample sheet.
    @param ok_sample_regions: A set of (sample_name, region, qcut) tuples that
        were given a good score by the pipeline.
    @return an array of JSON hashes, one for each conseq.
    """

    result = []
    ss = sample_sheet
    sequencings = run['sequencing_summary']
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
    target_regions = set()  # set([(project_name, tags)])
    for entry in sequencings:
        seeds = projects.getProjectSeeds(entry['target_project'])
        for seed in seeds:
            target_regions.add((entry['tag'], seed))

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
        sample_tags = sample_info["tags"]
        # FIXME if row["sequence"] is blank we replace it with a dash.
        # Need Conan to make that row blank-able.
        curr_seq = row["sequence"] if len(row["sequence"]) > 0 else "-"
        sample_region = (fastq_filename, row["region"], row["q-cutoff"])
        ok_region = sample_region in ok_sample_regions
        is_target_region = (sample_tags, row["region"]) in target_regions
        ok_for_release = ok_region and is_target_region
        result.append({"samplename": orig_sample_name,
                       # July 9, 2014: we can't do this properly right now
                       # without a lookup table that is yet to be fully
                       # defined.
                       "testcode": None,
                       "conseq_cutoff": row["consensus-percent-cutoff"],
                       "region": row["region"],
                       "qcutoff": float(row["q-cutoff"]),
                       "snum": fastq_filename.split('_')[-1],
                       "seq": curr_seq,
                       "ok_for_release": ok_for_release})
    return result


def build_hla_b_seqs(sample_file):
    """
    Build JSON hashes for HLA-B sequence records

    @param sample_file: open file that holds the variant info
    """

    result = []
    expected_exon_prefix = 'HLA-B-exon'
    # sample,seed,qcut,region,index,count,seq
    rows = csv.DictReader(sample_file)
    for row in rows:
        ind = int(row['index'])

        sample_name = row['sample']
        exon = row['region']
        if exon.startswith(expected_exon_prefix):
            exon_number = int(exon[len(expected_exon_prefix):])
        else:
            raise ValueError('Unexpected exon {!r}', exon)
        qcutoff = row['qcut']
        cnt = row['count']
        curr_seq = row['seq']

        result.append({'samplename': sample_name,
                       'testcode': None,
                       'exon': exon_number,
                       'qcutoff': qcutoff,
                       'ind': ind,
                       'cnt': cnt,
                       'string': curr_seq})

    return result


def build_review_decisions(coverage_file,
                           collated_counts_file,
                           sample_sheet,
                           sequencings,
                           project_regions,
                           regions):
    """ Build a list of request objects that will create the review decision
    records.

    @param coverage_file: CSV file with coverage scores
    @param collated_counts_file: CSV file with read counts
    @param sample_sheet: the sample sheet for the run
    @param sequencings: the sequencing records from QAI
    @param project_regions: [{"id": project_region_id,
                              "project_name": project_name,
                              "seed_region_names": [seed_region_name],
                              "coordinate_region_name": coordinate_region_name}]
    @param regions: [{"id": region_id, "name": region_name}]
    """

    project_region_map = dict(
        [((entry['project_name'], entry['coordinate_region_name']), entry['id'])
         for entry in project_regions])
    region_map = dict([(entry['name'], entry['id']) for entry in regions])
    sample_tags = dict(map(itemgetter('filename', 'tags'), sample_sheet['DataSplit']))
    sample_names = dict(map(itemgetter('tags', 'filename'), sample_sheet['DataSplit']))

    counts_map = {}  # {tags: raw, (tags, seed): mapped]}
    unreported_tags = set()
    # sample,type,count
    for counts in csv.DictReader(collated_counts_file):
        count = int(counts['count'])
        tags = sample_tags[counts['sample']]
        count_type = counts['type']
        if count_type == 'raw':
            counts_map[tags] = count
            unreported_tags.add(tags)
        elif count_type != 'unmapped':
            seed = count_type.split(' ', 1)[1]
            key = tags, seed
            current_count = counts_map.get(key, 0)
            counts_map[key] = max(current_count, count)

    sequencing_map = defaultdict(dict)  # {tags: {project: sequencing}}
    for sequencing in sequencings:
        sequencing_map[sequencing['tag']][sequencing['target_project']] = sequencing

    targeted_projects = set(map(itemgetter('target_project'), sequencings))

    decisions = {}  # {(sample_name, region): decision}
    # sample,project,region,q.cut,min.coverage,which.key.pos,off.score,on.score
    for coverage in csv.DictReader(coverage_file):
        tags = sample_tags[coverage['sample']]
        project_map = sequencing_map.get(tags)
        if project_map is None:
            raise KeyError("No sequencing found with tags '%s'." % tags)
        sequencing = project_map.get(coverage['project'])
        if sequencing is not None:
            score = int(coverage['on.score'])
        else:
            score = int(coverage['off.score'])
            first_project = sorted(project_map.iterkeys())[0]
            sequencing = project_map[first_project]
        project_region_id = project_region_map[(
            coverage['project'],
            coverage['region'])]
        raw_count = counts_map[tags]
        seed = coverage['seed']
        mapped_count = counts_map.get((tags, seed))
        seed_region_id = region_map[seed]

        decision_key = (coverage['sample'], coverage['region'])
        previous_decision = decisions.get(decision_key)
        is_replacement = (previous_decision is None or
                          score > previous_decision['score'] or
                          (score == previous_decision['score'] and
                           coverage['project'] in targeted_projects))
        if is_replacement:
            decisions[decision_key] = {
                'sequencing_id': sequencing['id'],
                'project_region_id': project_region_id,
                'seed_region_id': seed_region_id,
                'sample_name': coverage['sample'],
                'score': score,
                'min_coverage': int(coverage['min.coverage']),
                'min_coverage_pos': int(coverage['which.key.pos']),
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
    return decisions.values()


def upload_review_to_qai(coverage_file,
                         collated_counts_file,
                         run,
                         sample_sheet,
                         conseqs,
                         hla_b_seqs,
                         session):
    """ Create a review.

    @param coverage_file: the coverage scores to upload
    @param collated_counts_file: CSV file of read counts to upload
    @param run: a hash with the attributes of the run record, including a
        sequencing summary of all the samples and their target projects
    @param sample_sheet: details of the run so we can tell which sample used
        which tags
    @param conseqs: an array of JSON hashes to pass to QAI for the conseq
        child records
    @param hla_b_seqs: an array of JSON hashes to pass to QAI for the hla_b_seq
        child records
    @param session: the QAI session
    """

    runid = run['id']
    sequencings = run['sequencing_summary']

    project_regions = session.get_json(
        "/lab_miseq_project_regions?pipeline=" + settings.pipeline_version)
    if not project_regions:
        raise RuntimeError('Unknown pipeline: ' + settings.pipeline_version)

    regions = session.get_json("/lab_miseq_regions")

    decisions = build_review_decisions(coverage_file,
                                       collated_counts_file,
                                       sample_sheet,
                                       sequencings,
                                       project_regions,
                                       regions)

    session.post_json("/lab_miseq_reviews",
                      {'runid': runid,
                       'pipeline_id': find_pipeline_id(session),
                       'lab_miseq_review_decisions': decisions,
                       'lab_miseq_conseqs': conseqs,
                       'lab_miseq_hla_b_seqs': hla_b_seqs})


def clean_runname(runname):
    try:
        rundate = datetime.strptime(runname, '%d-%b-%y')
        clean_runname = datetime.strftime(rundate, '%d-%b-%Y')
    except ValueError:
        clean_runname = runname
    return clean_runname


def find_run(session, runname):
    """ Query QAI to find the run id for a given run name.

    @return: a hash with the attributes of the run record, including a
        sequencing summary of all the samples and their target projects.
    """
    cleaned_runname = clean_runname(runname)
    runs = session.get_json(
        "/lab_miseq_runs?summary=sequencing&runname=" + cleaned_runname)
    rowcount = len(runs)
    if rowcount == 0:
        raise RuntimeError("No run found with runname {!r}.".format(cleaned_runname))
    if rowcount != 1:
        raise RuntimeError("Found {} runs with runname {!r}.".format(
            rowcount,
            cleaned_runname))
    return runs[0]


def find_pipeline_id(session):
    """ Query QAI to find the pipeline id for the current version.

    @return: the pipeline id.
    """
    pipelines = session.get_json(
        "/lab_miseq_pipelines?version=" + settings.pipeline_version)
    rowcount = len(pipelines)
    if rowcount == 0:
        raise RuntimeError("No pipeline found with version {!r}.".format(
            settings.pipeline_version))
    if rowcount != 1:
        raise RuntimeError("Found {} pipelines with version {!r}.".format(
            rowcount,
            settings.pipeline_version))
    return pipelines[0]['id']


def load_ok_sample_regions(result_folder):
    ok_sample_regions = set()
    coverage_file = os.path.join(result_folder, 'coverage_scores.csv')
    with open(coverage_file, "rU") as f:
        reader = csv.DictReader(f)
        for row in reader:
            if row['on.score'] == '4':
                ok_sample_regions.add((row['sample'], row['region'], row['q.cut']))

    return ok_sample_regions


def process_folder(result_folder):
    logger.info('Uploading data to Oracle from {}'.format(result_folder))
    collated_conseqs = os.path.join(result_folder, 'conseq.csv')
    collated_counts = os.path.join(result_folder, 'remap_counts.csv')
    nuc_variants = os.path.join(result_folder, 'nuc_variants.csv')
    coverage_scores = os.path.join(result_folder, 'coverage_scores.csv')
    all_results_path, _ = os.path.split(os.path.normpath(result_folder))
    run_path, _ = os.path.split(all_results_path)
    sample_sheet_file = os.path.join(run_path, "SampleSheet.csv")
    with open(sample_sheet_file, "rU") as f:
        sample_sheet = sample_sheet_parser.sample_sheet_parser(f)

    ok_sample_regions = load_ok_sample_regions(result_folder)

    with qai_helper.Session() as session:
        session.login(settings.qai_path,
                      settings.qai_user,
                      settings.qai_password)

        run = find_run(session, sample_sheet["Experiment Name"])

        with open(collated_conseqs, "rU") as f:
            conseqs = build_conseqs(f,
                                    run,
                                    sample_sheet,
                                    ok_sample_regions)
        with open(nuc_variants, "rU") as f:
            hla_b_seqs = build_hla_b_seqs(f)
        with open(coverage_scores, "rU") as f, open(collated_counts, "rU") as f2:
            upload_review_to_qai(f,
                                 f2,
                                 run,
                                 sample_sheet,
                                 conseqs,
                                 hla_b_seqs,
                                 session)


def main():
    args, parser = parse_args()

    if args.result_folder:
        process_folder(args.result_folder)
    elif not args.load_all:
        parser.print_usage()
        exit(0)
    else:
        runs = glob(settings.rawdata_mount + 'MiSeq/runs/*/{}'.format(
            settings.NEEDS_PROCESSING))
        runs.sort()

        for run in runs:
            run_folder, _ = os.path.split(run)
            disabled_marker = os.path.join(run_folder, settings.ERROR_PROCESSING)
            if os.path.exists(disabled_marker):
                continue

            result_path = os.path.join(run_folder,
                                       'Results/version_' + settings.pipeline_version)

            if os.path.exists(result_path):
                # noinspection PyBroadException
                try:
                    process_folder(result_path)
                except Exception:
                    logger.error('Failed to process %s',
                                 result_path,
                                 exc_info=True)

    logger.info('Completed all uploads to Oracle.')


if __name__ == "__main__":
    main()
