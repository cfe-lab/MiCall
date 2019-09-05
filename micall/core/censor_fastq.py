#!/usr/bin/env python

import argparse
import csv
from gzip import GzipFile
import itertools
import math
import os


def parseArgs():
    parser = argparse.ArgumentParser(
        description='Censor tiles and cycles from a FASTQ file.')

    parser.add_argument('original_fastq',
                        type=argparse.FileType('rb'),
                        help='<input> FASTQ.gz containing original reads')
    parser.add_argument('bad_cycles_csv',
                        type=argparse.FileType('rU'),
                        help='<input> List of tiles and cycles rejected for poor quality')
    parser.add_argument('censored_fastq',
                        type=argparse.FileType('w'),
                        help='<output> uncompressed FASTQ containing censored reads')
    parser.add_argument('--unzipped',
                        '-u',
                        action='store_true',
                        help='Set if the FASTQ file is not compressed')

    return parser.parse_args()


def censor(src,
           bad_cycles_reader,
           dest,
           use_gzip=True,
           summary_file=None):
    """ Censor bases from a FASTQ file that were read in bad cycles.

    @param src: an open FASTQ file to read from
    @param bad_cycles_reader: an iterable collection of bad cycle entries:
        {'tile': tile, 'cycle': cycle}
    @param dest: an open FASTQ file to write to: censored bases will
        be written as 'N' with a quality '#'.
    @param summary_file: an open CSV file to write to: write a single row
        with the average read quality for the whole sample
    """
    bad_cycles = set()
    for cycle in bad_cycles_reader:
        bad_cycles.add((cycle['tile'], int(cycle['cycle'])))

    base_count = 0
    score_sum = 0.0
    if use_gzip:
        src = GzipFile(fileobj=src)
        dest = GzipFile(fileobj=dest)

    for ident, seq, opt, qual in itertools.zip_longest(src, src, src, src):
        # returns an aggregate of 4 lines per call
        ident_fields, read_fields = map(str.split, ident.decode('utf-8').split(' '), '::')
        tile = ident_fields[4]
        read_direction = read_fields[0]
        cycle_sign = 1 if read_direction == '1' else -1
        dest.write(ident)

        bad_count = 0
        for cycle, base in enumerate(seq.decode('utf-8').rstrip(), start=1):
            cycle = math.copysign(cycle, cycle_sign)
            if (tile, cycle) in bad_cycles:
                bad_count += 1
            else:
                if bad_count:
                    dest.write(b'N' * bad_count)
                    bad_count = 0
                dest.write(bytes(base, 'utf-8'))
        dest.write(b'\n')
        dest.write(opt)

        bad_count = 0
        for cycle, score in enumerate(qual.decode('utf-8').rstrip(), start=1):
            float_score = ord(score) - 33
            score_sum += float_score
            base_count += 1
            cycle = math.copysign(cycle, cycle_sign)
            if (tile, cycle) in bad_cycles:
                bad_count += 1
            else:
                if bad_count:
                    # write cumulative output of bad cycles
                    dest.write(b'#' * bad_count)
                    bad_count = 0
                dest.write(bytes(score, 'utf-8'))
        dest.write(b'\n')

    if summary_file is not None:
        avg_quality = score_sum/base_count if base_count > 0 else None
        summary = dict(base_count=base_count,
                       avg_quality=avg_quality)
        summary_writer = csv.DictWriter(summary_file,
                                        ['avg_quality', 'base_count'],
                                        lineterminator=os.linesep)
        summary_writer.writeheader()
        summary_writer.writerow(summary)


if __name__ == '__main__':
    args = parseArgs()

    censor(src=args.original_fastq,
           bad_cycles_reader=csv.DictReader(args.bad_cycles_csv),
           dest=args.censored_fastq,
           use_gzip=not args.unzipped)
elif __name__ == '__live_coding__':
    import unittest
    from micall.tests.censor_fastq_test import CensorTest

    suite = unittest.TestSuite()
    suite.addTest(CensorTest("testSummaryEmpty"))
    test_results = unittest.TextTestRunner().run(suite)

    print(test_results.errors)
    print(test_results.failures)
