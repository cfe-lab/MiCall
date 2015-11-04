#!/usr/bin/env python

import itertools
import argparse
import csv
import math
from gzip import GzipFile


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


def censor(original_file, bad_cycles_reader, censored_file, use_gzip=True):
    """ Censor bases from a FASTQ file that were read in bad cycles.

    @param original_file: an open FASTQ file to read from
    @param bad_cycles_reader: an iterable collection of bad cycle entries:
        {'tile': tile, 'cycle': cycle}
    @param censored_file: an open FASTQ file to write to: censored bases will
        be written as 'N' with a quality '#'.
    """
    bad_cycles = set()
    for cycle in bad_cycles_reader:
        bad_cycles.add((cycle['tile'], int(cycle['cycle'])))

    src = original_file
    dest = censored_file
    if use_gzip:
        src = GzipFile(fileobj=original_file)

    for ident, seq, opt, qual in itertools.izip_longest(src, src, src, src):
        # returns an aggregate of 4 lines per call
        ident_fields, read_fields = map(str.split, ident.split(' '), '::')
        tile = ident_fields[4]
        read_direction = read_fields[0]
        cycle_sign = 1 if read_direction == '1' else -1
        dest.write(ident)
        for cycle, base in enumerate(seq, start=1):
            cycle = math.copysign(cycle, cycle_sign)
            dest.write('N' if (tile, cycle) in bad_cycles else base)
        dest.write(opt)
        for cycle, score in enumerate(qual, start=1):
            cycle = math.copysign(cycle, cycle_sign)
            dest.write('#' if (tile, cycle) in bad_cycles else score)


if __name__ == '__main__':
    args = parseArgs()

    censor(args.original_fastq,
           csv.DictReader(args.bad_cycles_csv),
           args.censored_fastq,
           not args.unzipped)
