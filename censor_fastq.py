#!/usr/bin/env python

import itertools
import argparse
import csv
import math

def parseArgs():
    parser = argparse.ArgumentParser(
        description='Censor tiles and cycles from a FASTQ file.')
    
    parser.add_argument('original_fastq',
                        type=argparse.FileType('rU'),
                        help='<input> FASTQ containing original reads')
    parser.add_argument('bad_cycles_csv',
                        type=argparse.FileType('rU'),
                        help='<input> List of tiles and cycles rejected for poor quality')
    parser.add_argument('censored_fastq',
                        type=argparse.FileType('w'),
                        help='<output> FASTQ containing censored reads')
    
    return parser.parse_args()

def censor(original_file, bad_cycles_reader, censored_file):
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

    for ident, seq, opt, qual in itertools.izip_longest(original_file,
                                                        original_file,
                                                        original_file,
                                                        original_file):
        ident_fields, read_fields = map(str.split, ident.split(' '), '::')
        tile = ident_fields[4]
        read_direction = read_fields[0]
        cycle_sign = 1 if read_direction == '1' else -1
        censored_file.write(ident)
        for cycle, base in enumerate(seq, start=1):
            cycle = math.copysign(cycle, cycle_sign)
            censored_file.write('N' if (tile, cycle) in bad_cycles else base)
        censored_file.write(opt)
        for cycle, score in enumerate(qual, start=1):
            cycle = math.copysign(cycle, cycle_sign)
            censored_file.write('#' if (tile, cycle) in bad_cycles else score)

if __name__ == '__main__':
    args = parseArgs()

    censor(args.original_fastq,
           csv.DictReader(args.bad_cycles_csv),
           args.censored_fastq)
#     args.censored_fastq.write('test\n')
#     args.censored_fastq.close()
#     print 'Done with {}.'.format(args.censored_fastq.name)
