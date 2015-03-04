#! /usr/bin/env python

import argparse
import logging

import miseq_logging
import csv

def parseArgs():
    parser = argparse.ArgumentParser(
        description='Post-processing of short-read alignments.')
    
    parser.add_argument('quality_csv',
                        type=argparse.FileType('rU'),
                        help='QC error rate data')
    parser.add_argument('bad_cycles_csv',
                        type=argparse.FileType('w'),
                        help='List of tiles and cycles rejected for poor quality')
    
    return parser.parse_args()

logger = miseq_logging.init_logging_console_only(logging.DEBUG)

def main():
    args = parseArgs()
    with args.quality_csv, args.bad_cycles_csv:
        reader = csv.DictReader(args.quality_csv)
        writer = csv.DictWriter(args.bad_cycles_csv, ['tile', 'cycle'])
        writer.writeheader()
        for row in reader:
            errorrate = float(row['errorrate'])
            if errorrate > 7.5:
                writer.writerow(dict(tile=row['tile'],
                                     cycle=row['cycle']))

if __name__ == '__main__':
    main()