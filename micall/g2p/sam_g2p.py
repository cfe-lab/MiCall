#! /usr/bin/env python3.4

import argparse
from collections import Counter
import csv
import os
import re

from micall.core.sam2aln import merge_pairs, apply_cigar
from micall.utils.translation import translate

# screen for in-frame deletions
pat = re.compile('([A-Z])(---)+([A-Z])')

QMIN = 20   # minimum base quality within insertions
QCUT = 10   # minimum base quality to not be censored
QDELTA = 5
DEFAULT_MIN_COUNT = 3


def parse_args():
    parser = argparse.ArgumentParser(description='Calculate g2p scores from amino acid sequences.')

    parser.add_argument('remap_csv', type=argparse.FileType('rU'),
                        help='<input> CSV containing remap output (modified SAM)')
    parser.add_argument('nuc_csv', type=argparse.FileType('rU'),
                        help='<input> CSV containing nucleotide frequency output from aln2counts.py')
    parser.add_argument('g2p_csv', type=argparse.FileType('w'),
                        help='<output> CSV containing g2p predictions.')
    parser.add_argument('g2p_summary_csv', type=argparse.FileType('w'),
                        help='<output> CSV containing overall call for the sample.')
    return parser.parse_args()


class RegionTracker:
    def __init__(self, tracked_region):
        self.tracked_region = tracked_region
        self.ranges = {}

    def add_nuc(self, seed, region, query_pos):
        """
         # Add a nucleotide position to the tracker.
        :param seed: name of the seed region
        :param region: name of the coordinate region
        :param query_pos: query position in the consensus coordinates
        :return: unused
        """
        if region != self.tracked_region:
            return

        if seed in self.ranges:
            range = self.ranges[seed]
            if range[1] < query_pos:
                range[1] = query_pos
            elif query_pos < range[0]:
                range[0] = query_pos
        else:
            self.ranges.update({seed: [query_pos, query_pos]})

    def get_range(self, seed):
        """
        Get the minimum and maximum query positions that were seen for a seed.
        :param seed: name of the seed region
        :return: array of two integers
        """
        return self.ranges.get(seed, [None, None])


def sam_g2p(pssm, remap_csv, nuc_csv, g2p_csv, g2p_summary_csv=None, min_count=1):
    pairs = {}  # cache read for pairing
    merged = Counter()  # { merged_nuc_seq: count }
    tracker = RegionTracker('V3LOOP')

    # look up clipping region for each read
    reader = csv.DictReader(nuc_csv)
    for row in reader:
        if row['query.nuc.pos'] == '':
            # skip deletions in query relative to reference
            continue
        tracker.add_nuc(row['seed'], row['region'], int(row['query.nuc.pos'])-1)

    # parse contents of remap CSV output
    reader = csv.DictReader(remap_csv)
    for row in reader:
        clip_from, clip_to = tracker.get_range(row['rname'])
        if clip_from is None or row['cigar'] == '*':
            # uninteresting region
            continue

        seq2, qual2, ins2 = apply_cigar(row['cigar'],
                                        row['seq'],
                                        row['qual'],
                                        int(row['pos'])-1,
                                        clip_from,
                                        clip_to)

        mate = pairs.pop(row['qname'], None)
        if mate:
            seq1 = mate['seq']
            qual1 = mate['qual']
            ins1 = mate['ins']

            mseq = merge_pairs(seq1, seq2, qual1, qual2, ins1, ins2)

            merged[mseq] += 1

        else:
            pairs.update({row['qname']: {'seq': seq2,
                                         'qual': qual2,
                                         'ins': ins2}})

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
    for s, count in merged.most_common():
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
    row = {}
    row['rank'] = counts['rank']
    row['count'] = count
    seqlen = len(seq)
    if seq.upper().count('N') > (0.5*seqlen):
        # if more than 50% of the sequence is garbage
        row['error'] = 'low quality'
        return row

    if seqlen % 3 != 0:
        row['error'] = 'notdiv3'
        return row

    if seqlen == 0:
        row['error'] = 'zerolength'
        return row

    stats = {}
    prot = translate(seq,
                     ambig_char='X',
                     stats=stats)

    row['seq'] = prot
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


def main():
    args = parse_args()
    from micall.g2p.pssm_lib import Pssm
    pssm = Pssm()
    sam_g2p(pssm=pssm,
            remap_csv=args.remap_csv,
            nuc_csv=args.nuc_csv,
            g2p_csv=args.g2p_csv,
            g2p_summary_csv=args.g2p_summary_csv,
            min_count=DEFAULT_MIN_COUNT)

if __name__ == '__main__':
    # note, must be called from project root if executing directly
    # i.e., python micall/g2p/sam_g2p.py -h
    main()
elif __name__ == '__live_coding__':
    import unittest
    from micall.tests.sam_g2p_test import SamG2PTest

    suite = unittest.TestSuite()
    suite.addTest(SamG2PTest("testAmbiguousMixture"))
    test_results = unittest.TextTestRunner().run(suite)

    print(test_results.errors)
    print(test_results.failures)
