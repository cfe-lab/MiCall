#! /usr/bin/env python3.4

import argparse
import csv
import itertools
import logging
import math
from operator import itemgetter
import os

from micall.core import miseq_logging

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


def direction_grouper(cycle):
    return math.copysign(1, int(cycle['cycle']))


def report_bad_cycles(quality_csv, bad_cycles_csv, bad_tiles_csv=None):
    reader = csv.DictReader(quality_csv)
    writer = csv.DictWriter(bad_cycles_csv,
                            ['tile', 'cycle', 'errorrate'],
                            lineterminator=os.linesep)
    writer.writeheader()
    if bad_tiles_csv is not None:
        tile_writer = csv.DictWriter(bad_tiles_csv,
                                     ['tile', 'bad_cycles'],
                                     lineterminator=os.linesep)
        tile_writer.writeheader()
    for tile, tile_cycles in itertools.groupby(reader, itemgetter('tile')):
        bad_cycle_count = 0
        for _direction, cycles in itertools.groupby(tile_cycles,
                                                    direction_grouper):
            is_bad = False
            for cycle in cycles:
                errorrate = cycle['errorrate']
                is_bad = (is_bad or
                          errorrate is None or
                          errorrate == '' or
                          float(errorrate) >= BAD_ERROR_RATE)
                if is_bad:
                    writer.writerow(cycle)
                    bad_cycle_count += 1
        if bad_tiles_csv is not None:
            tile_writer.writerow(dict(tile=tile, bad_cycles=bad_cycle_count))


def main():
    args = parseArgs()
    with args.quality_csv, args.bad_cycles_csv:
        report_bad_cycles(args.quality_csv, args.bad_cycles_csv)

if __name__ == '__main__':
    main()
elif __name__ == '__live_coding__':
    import unittest
    from micall.tests.filter_quality_test import FilterQualityTest

    suite = unittest.TestSuite()
    suite.addTest(FilterQualityTest("test_tile_count"))
    test_results = unittest.TextTestRunner().run(suite)

    print(test_results.errors)
    print(test_results.failures)
