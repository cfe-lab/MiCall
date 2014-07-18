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


def insert_records(conseqs_file, version, ss, curs, runid):
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
    # for FASTQ_filename in ss["Data"]:
    #     sample_name = filename_re.match(FASTQ_filename).group(1)
    #     FASTQ_lookup[sample_name] = FASTQ_filename
    for row in conseqs_csv:
        # Each row of this file looks like:
        # sample,region,q-cutoff,s-number,consensus-percent-cutoff,sequence
        # We want to take the "sample" entry and get the corresponding
        # original Sample_Name from the sample sheet. In version 2, this
        # looks like [sample name]~[project name]#[...]
        # In version 1, this looked like [sample name]~[project name]#[...]
        # but both ; and _ got garbled by the MiSeq instrument itself.
        # Thus we have to work around it.
        FASTQ_filename = "{}_{}".format(row["sample"], row["s-number"])
        orig_sample_name = ss["Data"][FASTQ_filename]["orig_sample_name"]
        # FIXME if row["sequence"] is blank we replace it with a dash.
        # Need Conan to make that row blank-able.
        curr_seq = row["sequence"] if len(row["sequence"]) > 0 else "-"
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
                      "seq": curr_seq})


def find_runid(curs, runname):
    """ Query the database to find the run id for a given run name.
    
    Returns the run id.
    """
    sql = """
    SELECT  id
    FROM    lab_miseq_run
    WHERE   UPPER(runname) = UPPER(:runname)
    """
    try:
        rundate = datetime.strptime(runname, '%d-%b-%y')
        clean_runname = datetime.strftime(rundate, '%d-%b-%Y')
    except ValueError:
        clean_runname = runname
    curs.execute(sql, {'runname': clean_runname})
    rows = curs.fetchall()
    rowcount = len(rows)
    if rowcount == 0:
        raise RuntimeError("No run found with runname {!r}.".format(clean_runname))
    if rowcount != 1:
        raise RuntimeError("Found {} runs with runname {!r}.".format(
            rowcount,
            clean_runname))
    return rows[0][0]


def check_existing_data(curs, runid, logger):
    """ Check if there is already data for runid in the table.
    
    Delete it if there is any.
    """
    
    sql = """
    DELETE
    FROM   lab_miseq_conseq
    WHERE  runid = :runid
    """
    curs.execute(sql, {'runid': runid})
    if curs.rowcount > 0:
        logger.warn('Deleted {} existing records for run id {}'.format(
            curs.rowcount,
            runid))


def upload_conseqs_to_Oracle(conseqs_file,
                             sample_sheet_file,
                             version,
                             logger,
                             user=settings.oracle_uploader,
                             password=settings.oracle_uploader_pass,
                             dbname=settings.oracle_db):
    """
    Parses a Pipeline-produced conseq file and uploads to Oracle.

    By default the Oracle user is specified in settings_default.py.
    """
    # Parse the sample sheet into a dictionary.
    ss = sample_sheet_parser.sample_sheet_parser(sample_sheet_file)
    
    conn = DB.connect(user, password, dbname)
    curs = conn.cursor()
    try:
        
        runid = find_runid(curs, ss["Experiment Name"])
        check_existing_data(curs, runid, logger)
    
        insert_records(conseqs_file, version, ss, curs, runid)
    finally:
        curs.close()
        conn.commit()
        conn.close()

def process_folder(result_folder, logger):
    logger.info('Uploading data to Oracle from {}'.format(result_folder))
    collated_conseqs = os.path.join(result_folder, 'collated_conseqs.csv')
    all_results_path, _ = os.path.split(result_folder)
    run_path, _ = os.path.split(all_results_path)
    sample_sheet = os.path.join(run_path, "SampleSheet.csv")
    with open(collated_conseqs, "rU") as f, open(sample_sheet, "rU") as g:
        upload_conseqs_to_Oracle(f, g, settings.pipeline_version, logger)

def main():
    import argparse
    
    parser = argparse.ArgumentParser(
        version=settings.pipeline_version,
        description="Update the Oracle database with conseq information")
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--result_folder",
                        "-r",
                        help="Result folder that holds collated_conseqs.csv file")
    group.add_argument("--load_all",
                        "-a",
                        action="store_true",
                        help="load all folders under RAW_DATA that have results.")

    args = parser.parse_args()
    
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
