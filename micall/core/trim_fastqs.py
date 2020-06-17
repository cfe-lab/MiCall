#!/usr/bin/env python3.6

""" Censor reads based on phiX quality data, and also trim adapter sequences. """

import argparse
import csv
import logging
import shutil
from datetime import datetime
from gzip import GzipFile
import itertools
import math
import os

from io import TextIOWrapper
from itertools import chain
from pathlib import Path

import typing

from Bio import SeqIO

from micall.utils.externals import CutAdapt

logger = logging.getLogger(__name__)


class TrimSteps:
    censor = 'trim.censor'
    adapters = 'trim.adapters'
    primers = 'trim.primers'
    all = (censor, adapters, primers)


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


def trim(original_fastq_filenames: typing.Sequence[str],
         bad_cycles_filename: str,
         trimmed_fastq_filenames: typing.Sequence[str],
         use_gzip: bool = True,
         summary_file: typing.TextIO = None,
         skip: typing.Tuple[str] = ()):
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
    :param skip: names of steps to skip
    """
    if summary_file is None:
        summary_writer = None
    else:
        summary_writer = csv.DictWriter(summary_file,
                                        ['avg_quality', 'base_count'],
                                        lineterminator=os.linesep)
        summary_writer.writeheader()

    censored_filenames = [filename + '.censored.fastq'
                          for filename in trimmed_fastq_filenames]
    if (TrimSteps.censor in skip) or not os.path.exists(bad_cycles_filename):
        bad_cycles = []
    else:
        with open(bad_cycles_filename, 'r') as bad_cycles:
            bad_cycles = list(csv.DictReader(bad_cycles))
    start_time = datetime.now()
    for i, (src_name, dest_name) in enumerate(
            zip(original_fastq_filenames, censored_filenames)):

        cycle_sign = 1 - 2*i
        with open(src_name, 'rb') as src, open(dest_name, 'w') as dest:
            censor(src, bad_cycles, dest, use_gzip, summary_writer, cycle_sign)
    logger.debug('Finished censoring in %s for %s.',
                 datetime.now() - start_time,
                 trimmed_fastq_filenames[0])

    cut_all(Path(censored_filenames[0]),
            Path(censored_filenames[1]),
            Path(trimmed_fastq_filenames[0]),
            Path(trimmed_fastq_filenames[1]),
            skip)
    purge_temp_files(censored_filenames)


def cut_all(censored_fastq1: Path,
            censored_fastq2: Path,
            trimmed_fastq1: Path,
            trimmed_fastq2: Path,
            skip: typing.Tuple[str] = ()):
    dedapted_filenames = [Path(str(filename) + '.dedapted.fastq')
                          for filename in (censored_fastq1, censored_fastq2)]
    ltrimmed_filenames = [Path(str(filename) + '.ltrimmed.fastq')
                          for filename in (censored_fastq1, censored_fastq2)]
    rtrimmed_filenames = [Path(str(filename) + '.rtrimmed.fastq')
                          for filename in (censored_fastq1, censored_fastq2)]
    if TrimSteps.adapters in skip:
        dedapted_filenames[0].symlink_to(censored_fastq1)
        dedapted_filenames[1].symlink_to(censored_fastq2)
    else:
        start_time = datetime.now()
        cut_adapters(censored_fastq1,
                     censored_fastq2,
                     dedapted_filenames[0],
                     dedapted_filenames[1])
        logger.debug('Trimmed adapters in %s for %s.',
                     datetime.now() - start_time,
                     trimmed_fastq1)
    if TrimSteps.primers in skip:
        shutil.copy(str(dedapted_filenames[0]), str(trimmed_fastq1))
        shutil.copy(str(dedapted_filenames[1]), str(trimmed_fastq2))
    else:
        start_time = datetime.now()
        cut_left_primers(dedapted_filenames[0],
                         dedapted_filenames[1],
                         ltrimmed_filenames[0],
                         ltrimmed_filenames[1])
        logger.debug('Trimmed left primers in %s for %s.',
                     datetime.now() - start_time,
                     trimmed_fastq1)
        right_start_time = datetime.now()
        cut_right_primers(dedapted_filenames[0],
                          dedapted_filenames[1],
                          rtrimmed_filenames[0],
                          rtrimmed_filenames[1])
        logger.debug('Trimmed right primers in %s for %s.',
                     datetime.now() - right_start_time,
                     trimmed_fastq1)
        dimer_start_time = datetime.now()
        combine_primer_trimming(dedapted_filenames[0],
                                dedapted_filenames[1],
                                ltrimmed_filenames[0],
                                ltrimmed_filenames[1],
                                rtrimmed_filenames[0],
                                rtrimmed_filenames[1],
                                trimmed_fastq1,
                                trimmed_fastq2)
        logger.debug('Combined primer trimming in %s for %s.',
                     datetime.now() - dimer_start_time,
                     trimmed_fastq1)
        logger.debug('Trimmed all primers in %s for %s.',
                     datetime.now() - start_time,
                     trimmed_fastq1)
    purge_temp_files(dedapted_filenames)
    purge_temp_files(ltrimmed_filenames)
    purge_temp_files(rtrimmed_filenames)


def purge_temp_files(filenames):
    for filename in filenames:
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


def cut_left_primers(original_fastq1: Path,
                     original_fastq2: Path,
                     trimmed_fastq1: Path,
                     trimmed_fastq2: Path):
    """ Trim left primers off both sets of reads.

    The filtering options are a little tricky. --trimmed-only means that only
    the reads that got filtered are written into the ltrimmed.fastq file.
    --pair-filter=both means that both reads have to be untrimmed before they
    will get excluded from the ltrimmed.fastq file. In other words, if either
    read had a primer trimmed off, then both reads will get written to the
    ltrimmed.fastq file.
    """
    script_path = Path(__file__).parent
    run_cut_adapt(['--quiet',
                   '-g', f'file:{script_path/"primers_left.fasta"}',
                   '-A', f'file:{script_path/"primers_left_end.fasta"}',
                   '-o', str(trimmed_fastq1),
                   '-p', str(trimmed_fastq2),
                   '--overlap=5',
                   '--trimmed-only',
                   '--pair-filter=both',
                   str(original_fastq1),
                   str(original_fastq2)])


def cut_right_primers(original_fastq1: Path,
                      original_fastq2: Path,
                      trimmed_fastq1: Path,
                      trimmed_fastq2: Path):
    script_path = Path(__file__).parent
    run_cut_adapt(['--quiet',
                   '-a', f'file:{script_path/"primers_right_end.fasta"}',
                   '-G', f'file:{script_path/"primers_right.fasta"}',
                   '-o', str(trimmed_fastq1),
                   '-p', str(trimmed_fastq2),
                   '--overlap=5',
                   '--trimmed-only',
                   '--pair-filter=both',
                   str(original_fastq1),
                   str(original_fastq2)])


def combine_primer_trimming(original_fastq1: Path,
                            original_fastq2: Path,
                            ltrimmed_fastq1: Path,
                            ltrimmed_fastq2: Path,
                            rtrimmed_fastq1: Path,
                            rtrimmed_fastq2: Path,
                            trimmed_fastq1: Path,
                            trimmed_fastq2: Path):
    trimmed_sequences1 = cut_primer_dimer_sequences(
        original_fastq1,
        start_trimmed_fastq=ltrimmed_fastq1,
        end_trimmed_fastq=rtrimmed_fastq1)
    trimmed_sequences2 = cut_primer_dimer_sequences(
        original_fastq2,
        start_trimmed_fastq=rtrimmed_fastq2,
        end_trimmed_fastq=ltrimmed_fastq2)
    with trimmed_fastq1.open('w') as f1, trimmed_fastq2.open('w') as f2:
        for seq1, seq2 in zip(trimmed_sequences1, trimmed_sequences2):
            if len(seq1) > 0 and len(seq2) > 0:
                SeqIO.write([seq1], f1, 'fastq')
                SeqIO.write([seq2], f2, 'fastq')


def cut_primer_dimer_sequences(original_fastq: Path,
                               start_trimmed_fastq: Path,
                               end_trimmed_fastq: Path):
    start_trimmed_seqs = chain(SeqIO.parse(start_trimmed_fastq, 'fastq'), [None])
    end_trimmed_seqs = chain(SeqIO.parse(end_trimmed_fastq, 'fastq'), [None])
    next_start_trimmed_seq = next(start_trimmed_seqs)
    next_end_trimmed_seq = next(end_trimmed_seqs)
    for original_seq in SeqIO.parse(original_fastq, 'fastq'):
        if (next_start_trimmed_seq is None or
                next_start_trimmed_seq.id != original_seq.id):
            start = 0
        else:
            start = len(original_seq) - len(next_start_trimmed_seq)
            next_start_trimmed_seq = next(start_trimmed_seqs)

        if (next_end_trimmed_seq is None or
                next_end_trimmed_seq.id != original_seq.id):
            end = None
        else:
            end = len(next_end_trimmed_seq)
            next_end_trimmed_seq = next(end_trimmed_seqs)
        yield original_seq[start:end]


def run_cut_adapt(cut_adapt_args):
    cut_adapt = CutAdapt()
    cut_adapt.check_output(cut_adapt_args)


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
