#! /usr/bin/env python3.4

import argparse
from collections import Counter
import csv
from operator import itemgetter
import os
import re

from gotoh import align_it

from micall.core.sam2aln import merge_pairs
from micall.utils.translation import translate, reverse_and_complement

# screen for in-frame deletions
pat = re.compile('([A-Z])(---)+([A-Z])')

QMIN = 20   # minimum base quality within insertions
QCUT = 10   # minimum base quality to not be censored
QDELTA = 5
DEFAULT_MIN_COUNT = 3
GAP_OPEN_COST = 10
GAP_EXTEND_COST = 3
USE_TERMINAL_COST = 1
MIN_PAIR_ALIGNMENT_SCORE = 20
# V3LOOP region of HIV HXB2, GenBank accession number K03455
V3LOOP_REF = ('TGTACAAGACCCAACAACAATACAAGAAAAAGAATCCGTATCCAGAGAGGACCAGGGA'
              'GAGCATTTGTTACAATAGGAAAAATAGGAAATATGAGACAAGCACATTGT')
# V3LOOP_AA = translate(V3LOOP_REF)
# PSSM_AREF = 'CTRPNXNNTXXRKSIRIXXXGPGQXXXAFYATXXXXGDIIGDIXXRQAHC'.replace('X', '')
MIN_V3_ALIGNMENT_SCORE = len(V3LOOP_REF) // 2


def parse_args():
    parser = argparse.ArgumentParser(description='Calculate g2p scores from amino acid sequences.')

    parser.add_argument('fastq1', type=argparse.FileType('rU'),
                        help='<input> FASTQ file containing read 1 reads')
    parser.add_argument('fastq2', type=argparse.FileType('rU'),
                        help='<input> FASTQ file containing read 2 reads')
    parser.add_argument('g2p_csv', type=argparse.FileType('w'),
                        help='<output> CSV containing g2p predictions.')
    parser.add_argument('g2p_summary_csv', type=argparse.FileType('w'),
                        help='<output> CSV containing overall call for the sample.')
    return parser.parse_args()


def fastq_g2p(pssm, fastq1, fastq2, g2p_csv, g2p_summary_csv=None, min_count=1):
    reader = FastqReader(fastq1, fastq2)
    merged_reads = merge_reads(reader)
    trimmed_reads = trim_reads(merged_reads)
    read_counts = count_reads(trimmed_reads)

    write_rows(pssm, read_counts, g2p_csv, g2p_summary_csv, min_count)


def write_rows(pssm, read_counts, g2p_csv, g2p_summary_csv=None, min_count=1):
    # apply g2p algorithm to merged reads
    g2p_writer = csv.DictWriter(
        g2p_csv,
        ['rank',
         'count',
         'g2p',
         'fpr',
         'call',
         'seq',
         'aligned',
         'error',
         'comment'],
        lineterminator=os.linesep)
    g2p_writer.writeheader()
    counts = Counter()
    skip_count = 0
    for s, count in read_counts:
        if count < min_count:
            skip_count += count
            continue
        # remove in-frame deletions
        seq = re.sub(pat, r'\g<1>\g<3>', s)

        row = _build_row(seq, count, counts, pssm)
        g2p_writer.writerow(row)
    if skip_count:
        counts['mapped'] += skip_count
        g2p_writer.writerow(dict(rank=counts['rank'] + 1,
                                 count=skip_count,
                                 error='count < {}'.format(min_count)))

    if g2p_summary_csv is not None:
        if counts['valid'] == 0:
            x4_pct_display = ''
            final_call = ''
        else:
            x4_pct = 100.0 * counts['x4'] / counts['valid']
            final_call = 'X4' if x4_pct >= 2.0 else 'R5'
            x4_pct_display = '{:0.2f}'.format(x4_pct)
        summary_writer = csv.writer(g2p_summary_csv, lineterminator=os.linesep)
        summary_writer.writerow(['mapped', 'valid', 'X4calls', 'X4pct', 'final'])
        summary_writer.writerow([counts['mapped'],
                                 counts['valid'],
                                 counts['x4'],
                                 x4_pct_display,
                                 final_call])


def _build_row(seq, count, counts, pssm):
    counts['rank'] += 1
    counts['mapped'] += count
    row = {'rank': counts['rank'],
           'count': count}
    seqlen = len(seq)
    if seq.upper().count('N') > (0.5*seqlen):
        # if more than 50% of the sequence is garbage
        row['error'] = 'low quality'
        return row

    if seqlen == 0:
        row['error'] = 'zerolength'
        return row

    stats = {}
    prot = translate(seq,
                     ambig_char='X',
                     stats=stats)
    row['seq'] = prot

    if seqlen % 3 != 0:
        row['error'] = 'notdiv3'
        if prot == 'CIRPNNNTRKVCI*DQDKHSMQQVK**GI*DEHI':
            exit(seq)
        return row

    # sanity check 1 - bounded by cysteines
    if not prot.startswith('C') or not prot.endswith('C'):
        row['error'] = 'cysteines'
        return row

    # sanity check 2 - too many ambiguous codons
    if stats['ambiguous'] > 1 or stats['max_aminos'] > 2:
        row['error'] = '> 2 ambiguous'
        return row

    # sanity check 3 - no stop codons
    if prot.count('*') > 0:
        row['error'] = 'stop codons'
        return row

    # sanity check 4 - V3 length in range 32-40 inclusive
    if stats['length'] < 32 or stats['length'] > 40:
        row['error'] = 'length'
        return row

    score, aligned = pssm.run_g2p(seq)
    if score is None:
        row['error'] = 'failed to align'
        return row

    aligned2 = ''.join([aa_list[0]
                        if len(aa_list) == 1 else '[%s]' % ''.join(aa_list)
                        for aa_list in aligned])

    fpr = pssm.g2p_to_fpr(score)
    if fpr > 3.5:
        call = 'R5'
    else:
        call = 'X4'
        counts['x4'] += count
    counts['valid'] += count
    row['g2p'] = '{:06f}'.format(score)
    row['fpr'] = fpr
    row['call'] = call
    row['seq'] = aligned2.replace('-', '')
    row['aligned'] = aligned2
    row['comment'] = 'ambiguous' if '[' in aligned2 else ''
    return row


class FastqError(Exception):
    pass


class FastqReader:
    def __init__(self, fastq1, fastq2):
        self.fastq1 = fastq1
        self.fastq2 = fastq2

    def __iter__(self):
        cache = {}  # {pair_name: (read_name, bases, qual)}
        read2_iter = iter(self.get_reads(self.fastq2))
        for pair_name, read1 in self.get_reads(self.fastq1):
            read2 = cache.pop(pair_name, None)
            while read2 is None:
                try:
                    pair_name2, read2 = next(read2_iter)
                except StopIteration:
                    raise FastqError('No match for read {}.'.format(pair_name))
                if pair_name2 != pair_name:
                    cache[pair_name2] = read2
                    read2 = None
            yield pair_name, read1, read2

    def get_reads(self, fastq):
        for header, bases, sep, quality in zip(fastq, fastq, fastq, fastq):
            assert header.startswith('@'), header
            assert sep.startswith('+'), sep
            pair_name, read_name = header[1:].split()
            yield pair_name, (read_name, bases.rstrip('\n'), quality.rstrip('\n'))


def merge_reads(reads):
    for pair_name, (_, seq1, qual1), (_, seq2, qual2) in reads:
        seq2 = reverse_and_complement(seq2)
        aligned1, aligned2, score = align_it(seq1,
                                             seq2,
                                             GAP_OPEN_COST,
                                             GAP_EXTEND_COST,
                                             USE_TERMINAL_COST)
        if score >= MIN_PAIR_ALIGNMENT_SCORE and aligned1[0] != '-':
            aligned_qual1 = align_quality(aligned1, qual1)
            aligned_qual2 = align_quality(aligned2, reversed(qual2))
            merged = merge_pairs(aligned1, aligned2, aligned_qual1, aligned_qual2)
            yield (pair_name, merged)


def align_quality(nucs, qual):
    qual_iter = iter(qual)
    qual = ''.join('!' if nuc == '-' else next(qual_iter)
                   for nuc in nucs)
    return qual


def trim_reads(reads):
    for pair_name, seq in reads:
        aligned_ref, aligned_seq, score = align_it(V3LOOP_REF,
                                                   seq,
                                                   GAP_OPEN_COST,
                                                   GAP_EXTEND_COST,
                                                   USE_TERMINAL_COST)
        if score < MIN_V3_ALIGNMENT_SCORE:
            continue
        left_padding = right_padding = 0
        for left_padding, nuc in enumerate(aligned_ref):
            if nuc != '-':
                break
        for right_padding, nuc in enumerate(reversed(aligned_ref)):
            if nuc != '-':
                break
        trimmed_read = aligned_seq[left_padding:(-right_padding or None)]
        stripped_read = trimmed_read.replace('-', '')
        yield pair_name, stripped_read


def count_reads(reads):
    counts = Counter(map(itemgetter(1), reads))
    return counts.most_common()


def main():
    args = parse_args()
    from micall.g2p.pssm_lib import Pssm
    pssm = Pssm()
    fastq_g2p(pssm=pssm,
              fastq1=args.fastq1,
              fastq2=args.fastq2,
              g2p_csv=args.g2p_csv,
              g2p_summary_csv=args.g2p_summary_csv,
              min_count=DEFAULT_MIN_COUNT)

if __name__ == '__main__':
    # note, must be called from project root if executing directly
    # i.e., python micall/g2p/fastq_g2p.py -h
    main()
