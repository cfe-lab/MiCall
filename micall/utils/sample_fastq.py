#!/usr/bin/env python

import argparse
import random


def parse_args():
    parser = argparse.ArgumentParser(
        description="Randomly sample reads from FASTQ files for quick processing.")
    parser.add_argument('fastq1',
                        type=argparse.FileType('r'),
                        help='original FASTQ file of forward reads')
    parser.add_argument('fastq2',
                        type=argparse.FileType('r'),
                        help='original FASTQ file of reverse reads')
    parser.add_argument('short_fastq1',
                        type=argparse.FileType('w'),
                        help='FASTQ file to write forward reads in')
    parser.add_argument('short_fastq2',
                        type=argparse.FileType('w'),
                        help='FASTQ file to write reverse reads in')
    parser.add_argument('--count',
                        '-n',
                        type=float,
                        default=10000.0,
                        help='approximate number of read pairs to write')
    return parser.parse_args()


def get_reads(fastq_file):
    """ Yield reads as tuples of four lines: header, sequence, '+', quality. """
    for read in zip(fastq_file, fastq_file, fastq_file, fastq_file):
        yield read


def get_named_reads(fastq_file):
    """ Yield (name, read) pairs. """
    for read in get_reads(fastq_file):
        header = read[0]
        name = header.split(' ')[0]
        yield (name, read)


def process_read(name, read, out_file, odds, skipped_names, chosen_names):
    """ Write a read to the out_file if it is chosen.

    @param name: the name of the read that is used to match forward and reverse
    reads
    @param read: a tuple of four lines that makes up a read
    @param out_file: an open file to write the chosen reads to
    @param odds: a float between zero and one that sets the odds of choosing
    a read
    @param skipped_names: a set of names that have already been skipped, but
    their partners have not been seen yet
    @param chosen_names: a set of names that have already been chosen and
    written, but their partners have not been seen yet
    """
    try:
        skipped_names.remove(name)
        # This name was skipped, and we've seen both reads. We're done.
        return
    except KeyError:
        pass
    try:
        chosen_names.remove(name)
        is_chosen = True
    except KeyError:
        # Haven't seen this name yet, decide whether to choose it.
        is_chosen = random.uniform(0, 1) < odds
        if is_chosen:
            chosen_names.add(name)
        else:
            skipped_names.add(name)

    if is_chosen:
        for line in read:
            out_file.write(line)


def main():
    args = parse_args()
    for line_count, _ in enumerate(args.fastq1, 1):
        pass
    args.fastq1.seek(0)
    read_count = line_count/4
    odds = args.count/read_count
    rev_reads = get_named_reads(args.fastq2)
    skipped_names = set()
    chosen_names = set()
    for fwd_name, fwd_read in get_named_reads(args.fastq1):
        rev_name, rev_read = rev_reads.next()
        process_read(fwd_name,
                     fwd_read,
                     args.short_fastq1,
                     odds,
                     skipped_names,
                     chosen_names)
        process_read(rev_name,
                     rev_read,
                     args.short_fastq2,
                     odds,
                     skipped_names,
                     chosen_names)

main()
