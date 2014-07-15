#!/usr/bin/env python

# Script to update the database with information from the
# collated_conseqs.csv files produced by the MiSeq pipeline.

# Usage:
# ./update_oracle.py [collated conseq CSV file] [sample sheet CSV] [version name]

import re
import sys
import csv
import cx_Oracle as DB
import settings

import sample_sheet_parser

def upload_conseqs_to_Oracle(conseqs_file, sample_sheet_file, version,
                             user=settings.oracle_uploader, password=settings.oracle_uploader_pass,
                             dbname=settings.oracle_db):
    """
    Parses a Pipeline-produced conseq file and uploads to Oracle.

    By default the Oracle user is specified in settings_default.py.
    """
    # Parse the sample sheet into a dictionary.
    ss = sample_sheet_parser.sample_sheet_parser(sample_sheet_file)
    
    conn = DB.connect(user, password, dbname)
    curs = conn.cursor()

    insert_statement = """
insert into specimen.LAB_MISEQ_CONSEQ
    (runname, samplename, testcode, pipelineversion, conseq_cutoff, region,
     qcutoff, snum, seq, enterdate)
values (:runname, :samplename, :testcode, :pipelineversion, :conseq_cutoff, :region,
     :qcutoff, :snum, :seq, SYSDATE)
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

    print("Uploading records...")
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

        curs.execute(
            insert_statement,
            {
                "runname": ss["Experiment Name"],
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
                "seq": curr_seq
            })
        
    print("... done.")
    
    curs.close()
    conn.commit()
    conn.close()


def main():
    import argparse
    
    parser = argparse.ArgumentParser(
        description="Update the Oracle database with conseq information")
    parser.add_argument("collated_conseqs",
                        help="CSV file containing all conseqs")
    parser.add_argument("sample_sheet",
                        help="sample sheet")
    parser.add_argument("version", type=str, help="Pipeline version used")

    args = parser.parse_args()

    ss_info = None
    with open(args.collated_conseqs, "rb") as f, open(args.sample_sheet, "rb") as g:
        upload_conseqs_to_Oracle(f, g, args.version)

if __name__ == "__main__":
    main()
