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

def upload_conseqs_to_qai(conseqs_file,
                          runid,
                          sample_sheet,
                          ok_sample_regions,
                          version,
                          logger,
                          session):
    """
    Parses a Pipeline-produced conseq file and uploads to Oracle.

    By default the Oracle user is specified in settings_default.py.
    @param conseqs_file: An open file that contains the consensus sequences
        from the counts2csf step for all samples in the run.
    @param sample_sheet: The data parsed from the sample sheet.
    @param ok_sample_regions: A set of (sample_name, region, qcut) tuples that
        were given a good score by the pipeline.
    @param version: the current version of the pipeline.
    @param logger: the logger to record events in.
    @param session: QAI session.
    """
    
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
    target_regions = set() # set([(project_name, tags)])
    for entry in ss["DataSplit"]:
        seeds = projects.getProjectSeeds(entry['project'])
        for seed in seeds:
            target_regions.add((entry['tags'], seed))
        
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
        send_json(session,
                  "/lab_miseq_runs/%d/lab_miseq_conseqs" % runid,
                  {"samplename": orig_sample_name,
                   # July 9, 2014: we can't do this properly right now 
                   # without a lookup table that is yet to be fully 
                   # defined.
                   "testcode": None,
                   "pipelineversion": version,
                   "conseq_cutoff": row["consensus-percent-cutoff"],
                   "region": row["region"],
                   "qcutoff": float(row["q-cutoff"]),
                   "snum": row["s-number"],
                   "seq": curr_seq,
                   "ok_for_release": ok_for_release})
        
def send_json(session, path, data):
    response = session.post(
        settings.qai_path + path,
        data=json.dumps(data),
        headers={'Content-Type': 'application/json',
                 'Accept': 'application/json'})
    response.raise_for_status()
    return response

def upload_hla_to_qai(sample_file, runid, version, session):
    """
    Insert HLA-B sequence records
    
    @param sample_file: open file that holds the variant info
    @param runid: the run id to insert all records under
    @param version: the current version of the pipeline.
    @param session: QAI session.
    """
    
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
    
        send_json(session,
                  "/lab_miseq_runs/%d/lab_miseq_hla_b_seqs" % runid,
                  {'samplename': sample_name,
                   'testcode': None,
                   'pipelineversion': version,
                   'exon': exon_number,
                   'qcutoff': qcutoff,
                   'ind': ind,
                   'cnt': cnt,
                   'string': curr_seq})
        
def build_review_decisions(coverage_file,
                           collated_counts_file,
                           sample_sheet,
                           sequencings,
                           projects):
    """ Build a list of request objects that will create the review decision
    records.
    
    @param coverage_file: CSV file with coverage scores
    @param collated_counts_file: CSV file with read counts
    @param sample_sheet: the sample sheet for the run
    @param sequencings: the sequencing records from QAI
    @param projects: ProjectConfig object
    """
    
    decisions = []
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
        project_region_id, seed = projects.findProjectRegion(
            coverage['project'],
            coverage['region'])
        raw_count = counts_map[tags]
        mapped_count = counts_map[(tags, seed)]
        decisions.append({'sequencing_id': sequencing['id'],
                          'project_region_id': project_region_id,
                          'score': score,
                          'min_coverage': int(coverage['min.coverage']),
                          'min_coverage_pos': int(coverage['which.key.pos']),
                          'raw_reads': raw_count,
                          'mapped_reads': mapped_count})
    return decisions

def upload_review_to_qai(coverage_file,
                         collated_counts_file,
                         run,
                         sample_sheet,
                         session):
    """ Create a review.
    
    @param coverage_file: the coverage scores to upload
    @param collated_counts_file: CSV file of read counts to upload
    @param run: a hash with the attributes of the run record, including a
        sequencing summary of all the samples and their target projects
    @param sample_sheet: details of the run so we can tell which sample used
        which tags
    @param session: the QAI session
    """
    
    runid = run['id']
    sequencings = run['sequencing_summary']
    
    response = send_json(session,
                         "/lab_miseq_reviews",
                         {'runid': runid,
                          'pipeline_id': find_pipeline_id(session)})
    review_id = response.json()['id']
    
    decisions = build_review_decisions(coverage_file,
                                       collated_counts_file,
                                       sample_sheet,
                                       sequencings,
                                       ProjectConfig.loadDefault())
    for decision in decisions:
        send_json(session,
                  "/lab_miseq_reviews/%d/lab_miseq_review_decisions" % review_id,
                  decision)

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
    response = session.get(
        settings.qai_path + "/lab_miseq_runs.json?summary=sequencing&runname=" + cleaned_runname)
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
    
def check_existing_data(session, sample_sheet, logger):
    """ Find the run id, and check for existing detail records.
    
    Delete the detail records if there are any.
    @param session: open QAI session
    @param sample_sheet: parsed data from the sample sheet
    @param logger: for reporting any existing records
    @return: a hash with the attributes of the run record, including a
        sequencing summary of all the samples and their target projects.
    """
    run = find_run(session, sample_sheet["Experiment Name"])
    runid = run['id']
    
    response = session.get(
        settings.qai_path + "/lab_miseq_runs/%d/lab_miseq_conseqs.json?pipelineversion=%s" %
        (runid, settings.pipeline_version))
    conseqs = response.json()
    for conseq in conseqs:
        response = session.delete(
            settings.qai_path + "/lab_miseq_runs/%d/lab_miseq_conseqs/%d" % (runid, conseq['id']))
        response.raise_for_status()
    
    response = session.get(
        settings.qai_path + "/lab_miseq_runs/%d/lab_miseq_hla_b_seqs.json?pipelineversion=%s" %
        (runid, settings.pipeline_version))
    hla_b_seqs = response.json()
    for hla_b_seq in hla_b_seqs:
        response = session.delete(
            settings.qai_path + "/lab_miseq_runs/%d/lab_miseq_hla_b_seqs/%d" % (runid, hla_b_seq['id']))
        response.raise_for_status()
    
    total_rowcount = len(conseqs) + len(hla_b_seqs)
    if total_rowcount > 0:
        logger.warn('Deleted {} existing records for run id {}'.format(
            total_rowcount,
            runid))
    
    return run

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
        run = check_existing_data(session, sample_sheet, logger)
        runid = run['id']

        with open(collated_conseqs, "rU") as f:
            upload_conseqs_to_qai(f,
                                  runid,
                                  sample_sheet,
                                  ok_sample_regions,
                                  settings.pipeline_version,
                                  logger,
                                  session)
        with open(nuc_variants, "rU") as f:
            upload_hla_to_qai(f, runid, settings.pipeline_version, session)
        with open(coverage_scores, "rU") as f, open(collated_counts, "rU") as f2:
            upload_review_to_qai(f,
                                 f2,
                                 run,
                                 sample_sheet,
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
