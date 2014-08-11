#! /usr/bin/env python

"""
Shipyard-style MiSeq pipeline, step 3
Takes SAM in CSV format as input.  Assumes that these data have been
through remap.py (but not strictly necessary!).  Merges paired-end reads
and outputs aligned sequences (with insertions and deletions, minus soft
clips).

Dependencies:
    bowtie2-build
    bowtie2-align
    samtools
    settings.py
"""

import argparse
import os
import sys
import itertools
import re

from settings import max_prop_N, read_mapping_cutoff, sam2csf_q_cutoffs

parser = argparse.ArgumentParser(
    description='Conversion of SAM data into aligned format.')
parser.add_argument('sam_csv', help='<input> SAM output of bowtie2 in CSV format')
parser.add_argument('output_csv', help='<output> CSV containing cleaned and merged reads')
parser.add_argument('failed_csv', help='<output> CSV containing reads that failed to merge')
args = parser.parse_args()


cigar_re = re.compile('[0-9]+[MIDNSHPX=]')  # CIGAR token
gpfx = re.compile('^[-]+')  # length of gap prefix


def apply_cigar (cigar, seq, qual):
    newseq = ''
    newqual = ''
    tokens = cigar_re.findall(cigar)
    if len(tokens) == 0:
        return None, None, None
    # Account for removing soft clipped bases on left
    shift = 0
    if tokens[0].endswith('S'):
        shift = int(tokens[0][:-1])
    left = 0
    for token in tokens:
        length = int(token[:-1])
        # Matching sequence: carry it over
        if token[-1] == 'M':
            newseq += seq[left:(left+length)]
            newqual += qual[left:(left+length)]
            left += length
        # Deletion relative to reference: pad with gaps
        elif token[-1] == 'D':
            newseq += '-'*length
            newqual += ' '*length  # Assign fake placeholder score (Q=-1)
        # Insertion relative to reference: skip it (excise it)
        elif token[-1] == 'I':
            left += length
            continue
        # Soft clipping leaves the sequence in the SAM - so we should skip it
        elif token[-1] == 'S':
            left += length
            continue
        else:
            print "Unable to handle CIGAR token: {} - quitting".format(token)
            sys.exit()
    return shift, newseq, newqual


def merge_pairs (seq1, seq2, qual1, qual2, q_cutoff=10, minimum_q_delta=5):
    mseq = ''
    # force second read to be longest of the two
    if len(seq1) > len(seq2):
        seq1, seq2 = seq2, seq1
        qual1, qual2 = qual2, qual1

    for i, c2 in enumerate(seq2):
        q2 = ord(qual2[i])-33
        if i < len(seq1):
            c1 = seq1[i]
            q1 = ord(qual1[i])-33
            if c1 == '-' and c2 == '-':
                mseq += '-'
                continue
            if c1 == c2:  # Reads agree and at least one has sufficient confidence
                if q1 > q_cutoff or q2 > q_cutoff:
                    mseq += c1
                else:
                    mseq += 'N'  # neither base is confident
            else:
                if abs(q2 - q1) >= minimum_q_delta:
                    if q1 > max(q2, q_cutoff):
                        mseq += c1
                    elif q2 > max(q1, q_cutoff):
                        mseq += c2
                    else:
                        mseq += 'N'
                else:
                    mseq += 'N'  # cannot resolve between discordant bases
        else:
            # past end of read 1
            if c2 == '-' and q2 == 0:
                mseq += 'n'  # interval between reads
                continue
            mseq += c2 if (q2 > q_cutoff) else 'N'
    return mseq


def len_gap_prefix(s):
    hits = gpfx.findall(s)
    if hits:
        return len(hits[0])
    return 0


def main():
    # check that the inputs exist
    if not os.path.exists(args.sam_csv):
        print 'No input CSV found at', args.sam_csv
        sys.exit(1)

    # check that the output paths are valid
    output_path = os.path.split(args.output_csv)[0]
    if not os.path.exists(output_path) and output_path != '':
        print 'Output path does not exist:', output_path
        sys.exit(1)


    handle = open(args.sam_csv, 'rU')
    outfile = open(args.output_csv, 'w')
    failfile = open(args.failed_csv, 'w')

    for refname, group in itertools.groupby(handle, lambda x: x.split(',')[2]):
        aligned = dict([(qcut, {}) for qcut in sam2csf_q_cutoffs])
        cached_reads = {}  # for mate pairing
        for line in group:
            if refname == '0':
                print line
                sys.exit()

            qname, _, _, pos, mapq, cigar, _, _, _, seq, qual = line.strip('\n').split(',')
            if cigar == '*' or int(mapq) < read_mapping_cutoff:
                continue

            _, seq1, qual1 = apply_cigar(cigar, seq, qual)
            pos1 = int(pos)-1
            seq2 = '-'*pos1 + seq1  # pad sequence on left
            qual2 = '!'*pos1 + qual1  # assign lowest quality to gap prefix so it does not override mate

            if qname not in cached_reads:
                cached_reads.update({qname: (seq2, qual2)})
                continue

            # otherwise we have a matched pair
            seq1, qual1 = cached_reads.pop(qname)
            for qcut in sam2csf_q_cutoffs:
                mseq = merge_pairs(seq1, seq2, qual1, qual2, qcut)
                prop_N = mseq.count('N') / float(len(mseq.strip('-')))
                if prop_N > max_prop_N:
                    # merged read is too messy
                    failfile.write(','.join(map(str, [qname, qcut, seq1, qual1, seq2, qual2, prop_N, mseq])))
                    failfile.write('\n')
                    continue

                if mseq in aligned[qcut]:
                    aligned[qcut][mseq] += 1  # compress identical sequences
                else:
                    aligned[qcut].update({mseq: 1})

        # output compressed data
        for qcut, data in aligned.iteritems():
            # sort variants by count
            intermed = [(count, len_gap_prefix(s), s) for s, count in data.iteritems()]
            intermed.sort(reverse=True)
            for rank, (count, offset, seq) in enumerate(intermed):
                outfile.write(','.join(map(str, [refname, qcut, rank, count, offset, seq.strip('-')])))
                outfile.write('\n')

    handle.close()
    outfile.close()
    failfile.close()

if __name__ == '__main__':
    main()
