#!/usr/bin/env python3.6

""" Censor reads based on phiX quality data, and also trim adapter sequences. """

import argparse
import csv
import logging
from gzip import GzipFile
import itertools
import math
import os

from io import TextIOWrapper
from pathlib import Path

from cutadapt.__main__ import main as cutadapt_main


def parse_args():
    parser = argparse.ArgumentParser(
        description='Censor tiles and cycles from a FASTQ file.')

    parser.add_argument('original1_fastq',
                        help='<input> FASTQ.gz containing original reads (read 1)')
    parser.add_argument('original2_fastq',
                        help='<input> FASTQ.gz containing original reads (read 2)')
    parser.add_argument('bad_cycles_csv',
                        help='<input> List of tiles and cycles rejected for poor quality')
    parser.add_argument('trimmed1_fastq',
                        help='<output> uncompressed FASTQ containing trimmed reads (read 1)')
    parser.add_argument('trimmed2_fastq',
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

    dedapted_filenames = [filename + '.dedapted.fastq'
                          for filename in trimmed_fastq_filenames]
    censored_filenames = [filename + '.censored.fastq'
                          for filename in trimmed_fastq_filenames]
    if not os.path.exists(bad_cycles_filename):
        bad_cycles = []
    else:
        with open(bad_cycles_filename, 'r') as bad_cycles:
            bad_cycles = list(csv.DictReader(bad_cycles))
    for i, (src_name, dest_name) in enumerate(
            zip(original_fastq_filenames, censored_filenames)):
        cycle_sign = 1 - 2*i
        with open(src_name, 'rb') as src, open(dest_name, 'w') as dest:
            censor(src, bad_cycles, dest, use_gzip, summary_writer, cycle_sign)

    cut_adapters(censored_filenames[0],
                 censored_filenames[1],
                 dedapted_filenames[0],
                 dedapted_filenames[1])
    cut_primers(dedapted_filenames[0],
                dedapted_filenames[1],
                trimmed_fastq_filenames[0],
                trimmed_fastq_filenames[1])
    for filename in censored_filenames + dedapted_filenames:
        try:
            os.remove(filename)
        except OSError:
            # We tried to tidy up a temporary file, but it's not critical.
            pass


def cut_adapters(original_fastq1: Path,
                 original_fastq2: Path,
                 trimmed_fastq1: Path,
                 trimmed_fastq2: Path):
    script_path = os.path.dirname(__file__)
    adapter_files = [os.path.join(script_path, 'adapters_read{}.fasta'.format(i))
                     for i in (1, 2)]
    run_cut_adapt(['--quiet',
                   '-a', 'file:' + adapter_files[0],
                   '-A', 'file:' + adapter_files[1],
                   '-o', str(trimmed_fastq1),
                   '-p', str(trimmed_fastq2),
                   str(original_fastq1),
                   str(original_fastq2)])


def cut_primers(original_fastq1: Path,
                original_fastq2: Path,
                trimmed_fastq1: Path,
                trimmed_fastq2: Path):
    script_path = Path(__file__).parent
    run_cut_adapt(['--quiet',
                   '-g', f'file:{script_path/"primers_left.fasta"}',
                   '-a', f'file:{script_path/"primers_right.fasta"}',
                   '-G', f'file:{script_path/"primers_right_rev.fasta"}',
                   '-A', f'file:{script_path/"primers_left_rev.fasta"}',
                   '-o', str(trimmed_fastq1),
                   '-p', str(trimmed_fastq2),
                   '--overlap=5',
                   str(original_fastq1),
                   str(original_fastq2)])


def run_cut_adapt(cut_adapt_args):
    # Instead of launching a child process, run it in the current process.
    # Requires messing with the default logging.
    root_logger = logging.getLogger()
    original_level = root_logger.getEffectiveLevel()
    root_logger.setLevel(logging.CRITICAL)
    try:
        cutadapt_main(cut_adapt_args)
    finally:
        root_logger.setLevel(original_level)


def censor(original_file,
           bad_cycles_reader,
           censored_file,
           use_gzip=True,
           summary_writer=None,
           cycle_sign=1):
    """ Censor bases from a FASTQ file that were read in bad cycles.

    @param original_file: an open FASTQ file to read from
    @param bad_cycles_reader: an iterable collection of bad cycle entries:
        {'tile': tile, 'cycle': cycle}
    @param censored_file: an open FASTQ file to write to: censored bases will
        be written as 'N' with a quality '#'.
    @param use_gzip: True if the original file should be unzipped
    @param summary_writer: an open CSV DictWriter to write to: write a single row
        with the average read quality for the whole sample
    @param cycle_sign: +1 or -1, shows which direction to run cycles when
        checking bad_cycles.
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
    src = TextIOWrapper(src)

    for ident, seq, opt, qual in itertools.zip_longest(src, src, src, src):
        # returns an aggregate of 4 lines per call
        dest.write(ident)
        if not bad_cycles:
            dest.write(seq)
            tile = None
        else:
            ident_fields, read_fields = map(str.split, ident.split(' '), '::')
            tile = ident_fields[4]
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

    trim((args.original1_fastq, args.original2_fastq),
         args.bad_cycles_csv,
         (args.trimmed1_fastq, args.trimmed2_fastq),
         use_gzip=not args.unzipped)
