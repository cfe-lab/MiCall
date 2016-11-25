#!/usr/bin/env python

""" Censor reads based on phiX quality data, and also trim adapter sequences. """

import argparse
import csv
from gzip import GzipFile
import itertools
import math
import os

from micall.utils.externals import CutAdapt

# version of bowtie2, used for version control
CUT_ADAPT_VERSION = '1.11'
# path to executable, so you can install more than one version
CUT_ADAPT_PATH = 'cutadapt-' + CUT_ADAPT_VERSION


def parse_args():
    parser = argparse.ArgumentParser(
        description='Censor tiles and cycles from a FASTQ file.')

    parser.add_argument('original1_fastq',
                        type=argparse.FileType('rb'),
                        help='<input> FASTQ.gz containing original reads (read 1)')
    parser.add_argument('original2_fastq',
                        type=argparse.FileType('rb'),
                        help='<input> FASTQ.gz containing original reads (read 2)')
    parser.add_argument('bad_cycles_csv',
                        type=argparse.FileType('rU'),
                        help='<input> List of tiles and cycles rejected for poor quality')
    parser.add_argument('trimmed1_fastq',
                        type=argparse.FileType('w'),
                        help='<output> uncompressed FASTQ containing trimmed reads (read 1)')
    parser.add_argument('trimmed2_fastq',
                        type=argparse.FileType('w'),
                        help='<output> uncompressed FASTQ containing trimmed reads (read 2)')
    parser.add_argument('--unzipped',
                        '-u',
                        action='store_true',
                        help='Set if the original FASTQ files are not compressed')

    return parser.parse_args()


def trim(original_fastq_filenames,
         bad_cycles_filename,
         trimmed_fastq_filenames,
         use_gzip=True,
         summary_file=None):
    """

    :param original_fastq_filenames: sequence of two filenames, containing
        read 1 and read 2 in FASTQ format
    :param bad_cycles_filename: a CSV file with 'tile' and 'cycle' columns
    :param trimmed_fastq_filenames: sequence of two filenames, will write
        trimmed original reads in FASTQ format: adapters will be trimmed and
        censored bases will be written as 'N' with a quality '#'.
    :param use_gzip: True if the original file should be unzipped
    :param summary_file: an open CSV file to write to: write one row
        with the average read quality for each file
    """
    if summary_file is None:
        summary_writer = None
    else:
        summary_writer = csv.DictWriter(summary_file,
                                        ['avg_quality', 'base_count'],
                                        lineterminator=os.linesep)
        summary_writer.writeheader()
    cut_adapt = CutAdapt(CUT_ADAPT_VERSION, CUT_ADAPT_PATH)

    censored_filenames = [filename + '.censored.fastq'
                          for filename in trimmed_fastq_filenames]
    if not os.path.exists(bad_cycles_filename):
        bad_cycles = []
    else:
        with open(bad_cycles_filename, 'rU') as bad_cycles:
            bad_cycles = list(csv.DictReader(bad_cycles))
    for src_name, dest_name in zip(original_fastq_filenames, censored_filenames):
        with open(src_name, 'rb') as src, open(dest_name, 'w') as dest:
            censor(src, bad_cycles, dest, use_gzip, summary_writer)

    script_path = os.path.dirname(__file__)
    adapter_files = [os.path.join(script_path, 'adapters_read{}.fasta'.format(i))
                     for i in (1, 2)]
    cutadapt_args = ['-a', 'file:' + adapter_files[0],
                     '-A', 'file:' + adapter_files[1],
                     '-o', trimmed_fastq_filenames[0],
                     '-p', trimmed_fastq_filenames[1],
                     '--quiet',
                     censored_filenames[0],
                     censored_filenames[1]]
    cut_adapt.check_output(cutadapt_args)
    for filename in censored_filenames:
        os.remove(filename)


def censor(original_file,
           bad_cycles_reader,
           censored_file,
           use_gzip=True,
           summary_writer=None):
    """ Censor bases from a FASTQ file that were read in bad cycles.

    @param original_file: an open FASTQ file to read from
    @param bad_cycles_reader: an iterable collection of bad cycle entries:
        {'tile': tile, 'cycle': cycle}
    @param censored_file: an open FASTQ file to write to: censored bases will
        be written as 'N' with a quality '#'.
    @param use_gzip: True if the original file should be unzipped
    @param summary_writer: an open CSV DictWriter to write to: write a single row
        with the average read quality for the whole sample
    """
    bad_cycles = set()
    for cycle in bad_cycles_reader:
        bad_cycles.add((cycle['tile'], int(cycle['cycle'])))

    src = original_file
    dest = censored_file
    base_count = 0
    score_sum = 0.0
    if use_gzip:
        src = GzipFile(fileobj=original_file)

    for ident, seq, opt, qual in itertools.izip_longest(src, src, src, src):
        # returns an aggregate of 4 lines per call
        ident_fields, read_fields = map(str.split, ident.split(' '), '::')
        tile = ident_fields[4]
        read_direction = read_fields[0]
        cycle_sign = 1 if read_direction == '1' else -1
        dest.write(ident)
        bad_count = 0
        for cycle, base in enumerate(seq.rstrip(), start=1):
            cycle = math.copysign(cycle, cycle_sign)
            if (tile, cycle) in bad_cycles:
                bad_count += 1
            else:
                if bad_count:
                    dest.write('N' * bad_count)
                    bad_count = 0
                dest.write(base)
        dest.write('\n')
        dest.write(opt)
        bad_count = 0
        for cycle, score in enumerate(qual.rstrip(), start=1):
            float_score = ord(score) - 33
            score_sum += float_score
            base_count += 1
            cycle = math.copysign(cycle, cycle_sign)
            if (tile, cycle) in bad_cycles:
                bad_count += 1
            else:
                if bad_count:
                    dest.write('#' * bad_count)
                    bad_count = 0
                dest.write(score)
        dest.write('\n')
    if summary_writer is not None:
        avg_quality = score_sum/base_count if base_count > 0 else None
        summary = dict(base_count=base_count,
                       avg_quality=avg_quality)
        summary_writer.writerow(summary)


if __name__ == '__main__':
    args = parse_args()

    censor(args.original_fastq,
           csv.DictReader(args.bad_cycles_csv),
           args.censored_fastq,
           not args.unzipped)
elif __name__ == '__live_coding__':
    import unittest
    from micall.tests.trim_fastqs_test import CensorTest

    suite = unittest.TestSuite()
    suite.addTest(CensorTest("testSummaryEmpty"))
    test_results = unittest.TextTestRunner().run(suite)

    print(test_results.errors)
    print(test_results.failures)
