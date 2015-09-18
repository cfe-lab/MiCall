#! /usr/bin/env python

import argparse
from collections import Counter
from csv import DictReader
import re

from micall.core.sam2aln import merge_pairs, apply_cigar
from micall.utils.translation import translate

# screen for in-frame deletions
pat = re.compile('([A-Z])(---)+([A-Z])')

QMIN = 20   # minimum base quality within insertions
QCUT = 10   # minimum base quality to not be censored
QDELTA = 5

def parse_args():
    parser = argparse.ArgumentParser(description='Calculate g2p scores from amino acid sequences.')

    parser.add_argument('remap_csv', type=argparse.FileType('rU'),
                        help='<input> CSV containing remap output (modified SAM)')
    parser.add_argument('nuc_csv', type=argparse.FileType('rU'),
                        help='<input> CSV containing nucleotide frequency output from aln2counts.py')
    parser.add_argument('g2p_csv', type=argparse.FileType('w'),
                        help='<output> CSV containing g2p predictions.')
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


def sam_g2p(pssm, remap_csv, nuc_csv, g2p_csv):
    is_ruby_compatible = False # Used for comparing results with old Ruby code
    pairs = {}  # cache read for pairing
    merged = Counter()  # { merged_nuc_seq: count }
    tracker = RegionTracker('V3LOOP')

    # look up clipping region for each read
    reader = DictReader(nuc_csv)
    for row in reader:
        if row['query.nuc.pos'] == '':
            # skip deletions in query relative to reference
            continue
        tracker.add_nuc(row['seed'], row['region'], int(row['query.nuc.pos'])-1)

    # parse contents of remap CSV output
    reader = DictReader(remap_csv)
    for row in reader:
        clip_from, clip_to = tracker.get_range(row['rname'])
        if clip_from is None:
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
    g2p_csv.write('rank,count,g2p,fpr,aligned,error\n')  # CSV header
    rank = 0
    for s, count in merged.most_common():
        # remove in-frame deletions
        if is_ruby_compatible:
            seq = re.sub('---', '', s)
        else:
            seq = re.sub(pat, r'\g<1>\g<3>', s)

        rank += 1
        prefix = '%d,%d' % (rank, count)
        seqlen = len(seq)
        if seq.upper().count('N') > (0.5*seqlen):
            # if more than 50% of the sequence is garbage
            g2p_csv.write('%s,,,,low quality\n' % prefix)
            continue
     
        if seqlen % 3 != 0:
            g2p_csv.write('%s,,,,notdiv3\n' % prefix)
            continue
     
        if seqlen == 0:
            g2p_csv.write('%s,,,,zerolength\n' % prefix)
            continue

        prot = translate(seq, 0, ambig_char='X', translate_mixtures=False)

        # sanity check 1 - bounded by cysteines
        if not prot.startswith('C') or not prot.endswith('C'):
            g2p_csv.write('%s,,,%s,cysteines\n' % (prefix, prot))
            continue

        # sanity check 2 - too many ambiguous codons
        if is_ruby_compatible:
            if prot.count('X') > 0:
                g2p_csv.write('%s,,,%s,ambiguous\n' % (prefix, prot))
                continue
        else:
            if prot.count('X') > 2:
                g2p_csv.write('%s,,,%s,> 2 ambiguous\n' % (prefix, prot))
                continue

        # sanity check 3 - no stop codons
        if prot.count('*') > 0:
            g2p_csv.write('%s,,,%s,stop codons\n' % (prefix, prot))
            continue

        # sanity check 4 - V3 length in range 32-40 inclusive
        if len(prot) < 32 or len(prot) > 40:
            g2p_csv.write('%s,,,%s,length\n' % (prefix, prot))
            continue

        score, aligned = pssm.run_g2p(seq)
        if is_ruby_compatible:
            score = round(score, 10) if score is not None else score

        try:
            aligned2 = ''.join([aa_list[0] if len(aa_list) == 1 else '[%s]'%''.join(aa_list)
                                for aa_list in aligned])
        except:
            # sequence failed to align
            g2p_csv.write('%s,%s,,,failed to align\n' % (prefix, score))
            continue

        fpr = pssm.g2p_to_fpr(score)
        g2p_csv.write(','.join(map(str, [prefix,
                                         score,
                                         fpr,
                                         aligned2,
                                         'ambiguous' if '[' in aligned2 else '']))+'\n')


def main():
    args = parse_args()
    from micall.g2p.pssm_lib import Pssm
    pssm = Pssm()
    sam_g2p(pssm=pssm, remap_csv=args.remap_csv, nuc_csv=args.nuc_csv, g2p_csv=args.g2p_csv)

if __name__ == '__main__':
    # note, must be called from project root if executing directly
    # i.e., python micall/g2p/sam_g2p.py -h
    main()
