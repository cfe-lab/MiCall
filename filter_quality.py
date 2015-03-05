#! /usr/bin/env python

import argparse
import logging

import miseq_logging
import csv
import itertools
import math

BAD_ERROR_RATE = 7.5

def parseArgs():
    parser = argparse.ArgumentParser(
        description='Post-processing of short-read alignments.')
    
    parser.add_argument('quality_csv',
                        type=argparse.FileType('rU'),
                        help='QC error rate data, grouped by tile')
    parser.add_argument('bad_cycles_csv',
                        type=argparse.FileType('w'),
                        help='List of tiles and cycles rejected for poor quality')
    
    return parser.parse_args()

logger = miseq_logging.init_logging_console_only(logging.DEBUG)

def grouper(cycle):
    return (cycle['tile'], math.copysign(1, int(cycle['cycle'])))

def main():
    args = parseArgs()
    with args.quality_csv, args.bad_cycles_csv:
        reader = csv.DictReader(args.quality_csv)
        writer = csv.DictWriter(args.bad_cycles_csv,
                                ['tile', 'cycle', 'errorrate'])
        writer.writeheader()
        for _tile_direction, cycles in itertools.groupby(reader, grouper):
            is_bad = False
            for cycle in cycles:
                is_bad = is_bad or float(cycle['errorrate']) >= BAD_ERROR_RATE
                if is_bad:
                    writer.writerow(cycle)

if __name__ == '__main__':
    main()