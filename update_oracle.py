#!/usr/bin/env python

# Script to update the database with information from the
# collated_conseqs.csv files produced by the MiSeq pipeline.

# Usage:
# ./update_oracle.py [collated conseq CSV file] [sample sheet CSV] [version name]

import csv
import cx_Oracle as DB
from datetime import datetime
from glob import glob
import os
import settings

import miseq_logging
import sample_sheet_parser


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

def insert_conseq_records(conseqs_file,
                   version,
                   target_regions,
                   ok_sample_regions,
                   ss,
                   curs,
                   runid):
    """ Insert consensus sequence records.
    
    @param conseqs_file: open file that holds the consensus sequence for each
    sample and region.
    @param version: the current version of the pipeline.
    @param target_regions: a set of (tag, region) tuples.
    @param ok_sample_regions: A set of (sample_name, region, qcut) tuples that
        were given a good score by the pipeline.
    @param ss: the sample sheet data.
    @param curs: the open database cursor.
    @param runid: the run id to insert the records under.
    """
    insert_statement = """
insert
into    LAB_MISEQ_CONSEQ
        (
        id,
        runid,
        runname,
        samplename,
        testcode,
        pipelineversion,
        conseq_cutoff,
        region,
        qcutoff,
        snum,
        seq,
        ok_for_release,
        enterdate
        )
values  (
        lab_miseq_conseq_seq.nextval,
        :runid,
        (
        SELECT  run.runname
        FROM    lab_miseq_run run
        WHERE   run.id = :runid
        ),
        :samplename,
        :testcode,
        :pipelineversion,
        :conseq_cutoff,
        :region,
        :qcutoff,
        :snum,
        :seq,
        :ok_for_release,
        SYSDATE
        )
"""
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
        curs.execute(insert_statement,
                     {"runid": runid,
                      "samplename": orig_sample_name,
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

def insert_hla_records(sample_file, runid, version, conn):
    """
    Insert HLA-B sequence records
    
    @param sample_file: open file that holds the variant info
    @param runid: the run id to insert all records under
    @param version: the current version of the pipeline.
    @param conn: the database connection. All changes will be committed before
        this method returns.
    """
    
    insert_statement = """
insert
into LAB_MISEQ_HLA_B_SEQ (
    id,
    run_id,
    samplename,
    testcode,
    pipelineversion,
    exon,
    qcutoff,
    ind,
    cnt,
    string,
    enterdate
)
values (
    lab_miseq_hla_b_seq_seq.nextval,
    :runid,
    :samplename,
    :testcode,
    :pipelineversion,
    :exon,
    :qcutoff,
    :ind,
    :cnt,
    :string,
    SYSDATE
)
"""

    expected_exon_prefix = 'exon'
    curs = conn.cursor()
    try:
        rows = csv.DictReader(sample_file)
        for row in rows:
            if row['refname'] != 'HLA-B':
                continue
                
            sample_name = row['sample']
            exon = row['subregion']
            if exon.startswith(expected_exon_prefix):
                exon_number = int(exon[len(expected_exon_prefix):])
            else:
                raise ValueError('Unexpected exon {!r}', exon)
            qcutoff = row['qcut']
            ind = row['index']
            cnt = row['count']
            curr_seq = row['seq']
    
            curs.execute(insert_statement, {'runid': runid,
                                            'samplename': sample_name,
                                            'testcode': None,
                                            'pipelineversion': version,
                                            'exon': exon_number,
                                            'qcutoff': qcutoff,
                                            'ind': ind,
                                            'cnt': cnt,
                                            'string': curr_seq})
        conn.commit()
    finally:
        curs.close()


def clean_runname(runname):
    try:
        rundate = datetime.strptime(runname, '%d-%b-%y')
        clean_runname = datetime.strftime(rundate, '%d-%b-%Y')
    except ValueError:
        clean_runname = runname
    return clean_runname

def find_runid(curs, runname):
    """ Query the database to find the run id for a given run name.
    
    Returns the run id.
    """
    sql = """
    SELECT  id
    FROM    lab_miseq_run
    WHERE   UPPER(runname) = UPPER(:runname)
    """
    cleaned_runname = clean_runname(runname)
    curs.execute(sql, {'runname': cleaned_runname})
    rows = curs.fetchall()
    rowcount = len(rows)
    if rowcount == 0:
        raise RuntimeError("No run found with runname {!r}.".format(cleaned_runname))
    if rowcount != 1:
        raise RuntimeError("Found {} runs with runname {!r}.".format(
            rowcount,
            cleaned_runname))
    return rows[0][0]

def find_target_regions(curs, runname):
    """ Query the database to find which sample regions were targeted.
    
    Finds rows in lab_miseq_sequencing that match runname and returns a set
    of tuples (tag, region) that were targeted. Tags are found in the sample
    sheet. For example, N501-N701. 
    """
    sql = """
    SELECT  tag,
            project_name
    FROM    lab_miseq_sequencing
    WHERE   UPPER(runname) = UPPER(:runname)
    AND     project_name != 'EMPTY'
    """
    project_regions = {'PR-RT': ['PR', 'RT'],
                       'HCV': ['HCV1A-H77-core',
                               'HCV1A-H77-E1',
                               'HCV1A-H77-E2',
                               'HCV1A-H77-p7',
                               'HCV1A-H77-NS2',
                               'HCV1A-H77-NS3',
                               'HCV1A-H77-NS4a',
                               'HCV1A-H77-NS4b',
                               'HCV1A-H77-NS5a',
                               'HCV1A-H77-NS5b']
                       }
    target_regions = set()
    cleaned_runname = clean_runname(runname)
    curs.execute(sql, {'runname': cleaned_runname})
    while True:
        row = curs.fetchone()
        if row is None:
            break
        
        tag, project = row
        regions = project_regions.get(project, [project])
        for region in regions:
            target_regions.add((tag, region))
    return target_regions

def check_existing_data(conn, sample_sheet, logger):
    """ Find the run id, and check for existing detail records.
    
    Delete the detail records if there are any.
    @param conn: open database connection
    @param sample_sheet: parsed data from the sample sheet
    @param logger: for reporting any existing records
    """
    
    conseq_sql = """
    DELETE
    FROM   lab_miseq_conseq
    WHERE  runid = :runid
    """
    hla_sql = """
    DELETE
    FROM   lab_miseq_hla_b_seq
    WHERE  run_id = :runid
    """
    curs = conn.cursor()
    try:
        runid = find_runid(curs, sample_sheet["Experiment Name"])
        curs.execute(conseq_sql, {'runid': runid})
        total_rowcount = curs.rowcount
        curs.execute(hla_sql, {'runid': runid})
        total_rowcount += curs.rowcount
        if total_rowcount > 0:
            logger.warn('Deleted {} existing records for run id {}'.format(
                total_rowcount,
                runid))
        conn.commit()
    finally:
        curs.close()
        
    return runid


def upload_conseqs_to_Oracle(conseqs_file,
                             runid,
                             sample_sheet,
                             ok_sample_regions,
                             version,
                             logger,
                             conn):
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
    @param conn: database connection. All changes will be committed before this
        method returns.
    """
    
    curs = conn.cursor()
    try:
        target_regions = find_target_regions(curs, 
                                             sample_sheet["Experiment Name"])
        
        insert_conseq_records(conseqs_file,
                       version,
                       target_regions,
                       ok_sample_regions,
                       sample_sheet,
                       curs,
                       runid)
        conn.commit()
    finally:
        curs.close()


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
    nuc_variants = os.path.join(result_folder, 'nuc_variants.csv')
    all_results_path, _ = os.path.split(result_folder)
    run_path, _ = os.path.split(all_results_path)
    sample_sheet_file = os.path.join(run_path, "SampleSheet.csv")
    with open(sample_sheet_file, "rU") as f:
        sample_sheet = sample_sheet_parser.sample_sheet_parser(f)
        
    ok_sample_regions = load_ok_sample_regions(result_folder)
    
    conn = DB.connect(settings.oracle_uploader,
                      settings.oracle_uploader_pass,
                      settings.oracle_db)
    try:
        runid = check_existing_data(conn, sample_sheet, logger)

        with open(collated_conseqs, "rU") as f:
            upload_conseqs_to_Oracle(f,
                                     runid,
                                     sample_sheet,
                                     ok_sample_regions,
                                     settings.pipeline_version,
                                     logger,
                                     conn)
        with open(nuc_variants, "rU") as f:
            insert_hla_records(f, runid, settings.pipeline_version, conn)
    finally:
        conn.close()


def main():
    args, parser = parse_args()
    
    logger = miseq_logging.init_logging_console_only()

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
