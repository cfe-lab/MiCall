#! /usr/bin/env python

import re
import argparse
from csv import DictReader
from micall.core.sam2aln import merge_pairs
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


def apply_cigar_and_clip (cigar, seq, qual, pos=0, clip_from=0, clip_to=None):
    """ Applies a cigar string to recreate a read, then clips the read.

    Use CIGAR string (Compact Idiosyncratic Gapped Alignment Report) in SAM data
    to apply soft clips, insertions, and deletions to the read sequence.
    
    @param cigar: a string in the CIGAR format, describing the relationship
        between the read sequence and the consensus sequence
    @param seq: the sequence that was read
    @param qual: quality codes for each base in the read
    @param pos: first position of the read, given in zero-based consensus
        coordinates
    @param clip_from: first position to include after clipping, given in
        zero-based consensus coordinates
    @param clip_to: last position to include after clipping, given in
        zero-based consensus coordinates, None means no clipping at the end
    @return: the new sequence and the new quality string. If none of the read
        was within the clipped range, then both strings will be blank.
    """
    newseq = '-' * int(pos)  # pad on left
    newqual = '!' * int(pos)
    is_valid = re.match(r'^((\d+)([MIDNSHPX=]))*$', cigar)
    tokens =   re.findall(r'(\d+)([MIDNSHPX=])', cigar)
    if not is_valid:
        raise RuntimeError('Invalid CIGAR string: {!r}.'.format(cigar))
    if len(tokens) == 0:
        return None, None, None
    left = 0
    for token in tokens:
        length, operation = token
        length = int(length)
        # Matching sequence: carry it over
        if operation == 'M':
            newseq += seq[left:(left+length)]
            newqual += qual[left:(left+length)]
            left += length
        # Deletion relative to reference: pad with gaps
        elif operation == 'D':
            newseq += '-'*length
            newqual += ' '*length  # Assign fake placeholder score (Q=-1)
        # Insertion relative to reference
        elif operation == 'I':
            its_quals = qual[left:(left+length)]
            if all([(ord(q)-33) >= QMIN for q in its_quals]):
                # accept only high quality insertions relative to sample consensus
                newseq += seq[left:(left+length)]
                newqual += its_quals
                if left + pos <= clip_from:
                    clip_from += length
                if clip_to:
                    clip_to += length  # accommodate insertion
            left += length
        # Soft clipping leaves the sequence in the SAM - so we should skip it
        elif operation == 'S':
            left += length
        else:
            raise RuntimeError('Unsupported CIGAR token: {!r}.'.format(
                ''.join(token)))
        if left > len(seq):
            raise RuntimeError(
                'CIGAR string {!r} is too long for sequence {!r}.'.format(cigar,
                                                                          seq))
    
    if left < len(seq):
        raise RuntimeError(
            'CIGAR string {!r} is too short for sequence {!r}.'.format(cigar,
                                                                       seq))
    
    end = None if clip_to is None else clip_to + 1
    return newseq[clip_from:end], newqual[clip_from:end]


def sam_g2p(pssm, remap_csv, nuc_csv, g2p_csv):
    pairs = {}  # cache read for pairing
    merged = {}  # tabular merged sequence variants
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

        seq2, qual2 = apply_cigar_and_clip(row['cigar'], row['seq'], row['qual'], int(row['pos'])-1, clip_from, clip_to)

        mate = pairs.get(row['qname'], None)
        if mate:
            seq1 = mate['seq']
            qual1 = mate['qual']

            # merge pairs only if they are the same length
            if len(seq1) == len(seq2):
                mseq = merge_pairs(seq1, seq2, qual1, qual2)
            else:
                # implies an insertion not covered in one read mate
                mseq = seq1 if len(seq1) > len(seq2) else seq2

            if mseq not in merged:
                merged.update({mseq: 0})
            merged[mseq] += 1

        else:
            pairs.update({row['qname']: {'seq': seq2, 'qual': qual2}})


    sorted = [(v,k) for k, v in merged.iteritems()]
    sorted.sort(reverse=True)

    # apply g2p algorithm to merged reads
    g2p_csv.write('rank,count,g2p,fpr,aligned,error\n')  # CSV header
    rank = 0
    for count, s in sorted:
        # remove in-frame deletions
        seq = re.sub(pat, r'\g<1>\g<3>', s)

        rank += 1
        prefix = '%d,%d' % (rank, count)
        seqlen = len(seq)
        if seq.upper().count('N') > (0.5*seqlen):
            # if more than 50% of the sequence is garbage
            g2p_csv.write('%s,,,,low quality\n' % prefix)
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
        if prot.count('X') > 2:
            g2p_csv.write('%s,,,%s,>2ambiguous\n' % (prefix, prot))
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
