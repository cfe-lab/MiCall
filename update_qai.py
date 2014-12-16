#!/usr/bin/env python

# Script to update QAI with information from the
# collated_conseqs.csv files produced by the MiSeq pipeline.

import csv
from datetime import datetime
from glob import glob
import json
import logging
import os
import requests

import miseq_logging
import sample_sheet_parser
import settings
from project_config import ProjectConfig


def parse_args():
    import argparse
    parser = argparse.ArgumentParser(
        version=settings.pipeline_version, 
        description="Update the Oracle database with conseq information")
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--result_folder", "-r", help="Result folder that holds collated_conseqs.csv file")
    group.add_argument("--load_all", "-a", action="store_true", help="load all folders under RAW_DATA that have results.")
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
    target_regions = set() # set([(project_name, tags)])
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
        fastq_filename = "{}_{}".format(row["sample"], row["s-number"])
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
                       "snum": row["s-number"],
                       "seq": curr_seq,
                       "ok_for_release": ok_for_release})
    return result
        
def post_json(session, path, data):
    response = session.post(
        settings.qai_path + path,
        data=json.dumps(data),
        headers={'Content-Type': 'application/json',
                 'Accept': 'application/json'})
    response.raise_for_status()
    return response

def get_json(session, path):
    response = session.get(
        settings.qai_path + path,
        headers={'Accept': 'application/json'})
    response.raise_for_status()
    return response

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
                           project_regions):
    """ Build a list of request objects that will create the review decision
    records.
    
    @param coverage_file: CSV file with coverage scores
    @param collated_counts_file: CSV file with read counts
    @param sample_sheet: the sample sheet for the run
    @param sequencings: the sequencing records from QAI
    @param project_regions: [{"id": project_region_id,
                              "project_name": project_name,
                              "seed_region_name": seed_region_name,
                              "coordinate_region_name": coordinate_region_name}]
    """
    
    project_region_map = dict(
        [((entry['project_name'], entry['coordinate_region_name']),
          (entry['id'], entry['seed_region_name']))
         for entry in project_regions])
    sample_tags = {}
    for sample in sample_sheet['DataSplit']:
        sample_tags[sample['filename']] = sample['tags']
    
    counts_map = {} # {tags: raw, (tags, seed): mapped]}
    # sample_name,type,count
    for counts in csv.DictReader(collated_counts_file):
        count = int(counts['count'])
        tags = sample_tags[counts['sample_name']]
        count_type = counts['type']
        if count_type == 'raw':
            counts_map[tags] = count
        elif count_type != 'unmapped':
            seed = count_type.split(' ', 1)[1]
            key = tags, seed
            current_count = counts_map.get(key, 0)
            counts_map[key] = max(current_count, count)
    
    decisions = {} # {(sample_name, region): decision}
    #sample,project,region,q.cut,min.coverage,which.key.pos,off.score,on.score
    for coverage in csv.DictReader(coverage_file):
        tags = sample_tags[coverage['sample']]
        score = int(coverage['off.score'])
        sequencing_with_target = None
        sequencing_with_tags = None
        for sequencing in sequencings:
            if sequencing['tag'] != tags:
                continue
            sequencing_with_tags = sequencing
            if sequencing['target_project'] != coverage['project']:
                continue
            sequencing_with_target = sequencing
            score = int(coverage['on.score'])
            break
        sequencing = sequencing_with_target or sequencing_with_tags
        if sequencing is None:
            raise StandardError("No sequencing found with tags '%s'." % tags)
        project_region_id, seed = project_region_map[(
            coverage['project'],
            coverage['region'])]
        raw_count = counts_map[tags]
        mapped_count = counts_map[(tags, seed)]
        decision_key = (coverage['sample'], coverage['region'])
        previous_decision = decisions.get(decision_key)
        if previous_decision is None or score > previous_decision['score']:
            decisions[decision_key] = {
                'sequencing_id': sequencing['id'],
                'project_region_id': project_region_id,
                'sample_name': coverage['sample'],
                'score': score,
                'min_coverage': int(coverage['min.coverage']),
                'min_coverage_pos': int(coverage['which.key.pos']),
                'raw_reads': raw_count,
                'mapped_reads': mapped_count
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
    
    response = get_json(session, "/lab_miseq_project_regions")
    project_regions = response.json()
    
    decisions = build_review_decisions(coverage_file,
                                       collated_counts_file,
                                       sample_sheet,
                                       sequencings,
                                       project_regions)
    
    post_json(session,
              "/lab_miseq_reviews",
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
    response = get_json(
        session,
        "/lab_miseq_runs?summary=sequencing&runname=" + cleaned_runname)
    runs = response.json()
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
    response = session.get(
        settings.qai_path + "/lab_miseq_pipelines.json?version=" + settings.pipeline_version)
    pipelines = response.json()
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

def process_folder(result_folder, logger):
    logger.info('Uploading data to Oracle from {}'.format(result_folder))
    collated_conseqs = os.path.join(result_folder, 'collated_conseqs.csv')
    collated_counts = os.path.join(result_folder, 'collated_counts.csv')
    nuc_variants = os.path.join(result_folder, 'nuc_variants.csv')
    coverage_scores = os.path.join(result_folder, 'coverage_scores.csv')
    all_results_path, _ = os.path.split(result_folder)
    run_path, _ = os.path.split(all_results_path)
    sample_sheet_file = os.path.join(run_path, "SampleSheet.csv")
    with open(sample_sheet_file, "rU") as f:
        sample_sheet = sample_sheet_parser.sample_sheet_parser(f)
        
    ok_sample_regions = load_ok_sample_regions(result_folder)
    
    with requests.Session() as session:
        response = session.post(settings.qai_path + "/account/login",
                                data={'user_login': settings.qai_user,
                                      'user_password': settings.qai_password})
        if response.status_code == requests.codes.forbidden:  # @UndefinedVariable
            exit('Login failed, check qai_user in settings.py')
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
    
    logger = miseq_logging.init_logging_console_only(logging.INFO)

    if args.result_folder:
        process_folder(args.result_folder, logger)
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
                try:
                    process_folder(result_path, logger)
                except:
                    logger.error('Failed to process %s',
                                 result_path,
                                 exc_info=True)    
    
    logger.info('Completed all uploads to Oracle.')
        

if __name__ == "__main__":
    main()
