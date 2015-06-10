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
import sys
import itertools
import re

from micall.settings import max_prop_N, read_mapping_cutoff, sam2aln_q_cutoffs
from csv import DictReader, DictWriter

def parseArgs():
    parser = argparse.ArgumentParser(
        description='Conversion of SAM data into aligned format.')
    parser.add_argument('remap_csv',
                        type=argparse.FileType('rU'),
                        help='<input> SAM output of bowtie2 in CSV format')
    parser.add_argument('aligned_csv',
                        type=argparse.FileType('w'),
                        help='<output> CSV containing cleaned and merged reads')
    parser.add_argument('insert_csv',
                        type=argparse.FileType('w'),
                        help='<output> CSV containing insertions relative to sample consensus')
    parser.add_argument('failed_csv',
                        type=argparse.FileType('w'),
                        help='<output> CSV containing reads that failed to merge')
    parser.add_argument('-p', type=int, default=None, help='(optional) number of threads')
    
    return parser.parse_args()


cigar_re = re.compile('[0-9]+[MIDNSHPX=]')  # CIGAR token
gpfx = re.compile('^[-]+')  # length of gap prefix


def apply_cigar (cigar, seq, qual):
    """
    Use CIGAR string (Compact Idiosyncratic Gapped Alignment Report) in SAM data
    to apply soft clips, insertions, and deletions to the read sequence.
    Any insertions relative to the sample consensus sequence are discarded to
    enforce a strict pairwise alignment, and returned separately in a
    dict object.
    """
    newseq = ''
    newqual = ''
    insertions = {}
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
            insertions.update({left: (seq[left:(left+length)], qual[left:(left+length)])})
            left += length
            continue
        # Soft clipping leaves the sequence in the SAM - so we should skip it
        elif token[-1] == 'S':
            left += length
            continue
        else:
            print "Unable to handle CIGAR token: {} - quitting".format(token)
            sys.exit()

    return shift, newseq, newqual, insertions


def merge_pairs (seq1, seq2, qual1, qual2, q_cutoff=10, minimum_q_delta=5):
    """
    Combine paired-end reads into a single sequence by managing discordant
    base calls on the basis of quality scores.
    """
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

        
def is_first_read(flag):
    """
    Interpret bitwise flag from SAM field.
    Returns True or False indicating whether the read is the first read in a pair.
    """
    IS_FIRST_SEGMENT = 0x40
    return (int(flag) & IS_FIRST_SEGMENT) != 0


def matchmaker(remap_csv):
    """
    An iterator that returns pairs of reads sharing a common qname from a remap CSV.
    Note that unpaired reads will be left in the cached_rows dictionary and
    discarded.
    :param remap_csv: open file handle to CSV generated by remap.py
    :return: yields pairs of rows from DictReader corresponding to paired reads
    """
    reader = DictReader(remap_csv)
    cached_rows = {}
    for row in reader:
        qname = row['qname']
        if qname not in cached_rows:
            cached_rows.update({qname: row})
            continue
        # current row should be the second read of the pair
        yield cached_rows[qname], row


def parse_sam(rows):
    """
    :param rows:  chunk of tuples of paired rows from matchmaker generator
    :return:
    """
    row1, row2 = rows
    failed_list = []
    insert_list = []
    rname = row1['rname']
    qname = row1['qname']

    cigar1 = row1['cigar']
    if cigar1 == '*' or int(row1['mapq']) < read_mapping_cutoff:
        return rname, None, insert_list

    pos1 = int(row1['pos'])-1  # convert 1-index to 0-index
    _, seq1, qual1, inserts = apply_cigar(cigar1, row1['seq'], row1['qual'])

    # report insertions relative to sample consensus
    for left, (iseq, iqual) in inserts.iteritems():
        insert_list.append([qname,
                            'F' if is_first_read(row1['flag']) else 'R',
                            rname,
                            pos1+left,
                            iseq,
                            iqual])

    seq1 = '-'*pos1 + seq1  # pad sequence on left
    qual1 = '!'*pos1 + qual1  # assign lowest quality to gap prefix so it does not override mate


    # now process the mate
    cigar2 = row2['cigar']
    if cigar2 == '*' or int(row2['mapq']) < read_mapping_cutoff:
        return rname, None, insert_list

    pos2 = int(row2['pos'])-1  # convert 1-index to 0-index
    _, seq2, qual2, inserts = apply_cigar(cigar2, row2['seq'], row2['qual'])
    for left, (iseq, iqual) in inserts.iteritems():
        insert_list.append([qname,
                            'F' if is_first_read(row2['flag']) else 'R',
                            rname,
                            pos2+left,
                            iseq,
                            iqual])
    seq2 = '-'*pos2 + seq2
    qual2 = '!'*pos2 + qual2

    # merge reads
    mseqs = {}
    for qcut in sam2aln_q_cutoffs:
        mseq = merge_pairs(seq1, seq2, qual1, qual2, qcut)
        prop_N = mseq.count('N') / float(len(mseq.strip('-')))
        if prop_N > max_prop_N:
            # fail read pair
            failed_list.append([qname, qcut, seq1, qual1, seq2, qual2, prop_N, mseq])
            continue
        mseqs.update({qcut: mseq})

    return rname, mseqs, insert_list, failed_list


def sam2aln(remap_csv, aligned_csv, insert_csv, failed_csv, nthreads=None):
    # prepare outputs
    insert_fields =  ['qname', 'fwd_rev', 'refname', 'pos', 'insert', 'qual']
    insert_writer = DictWriter(insert_csv, insert_fields)
    insert_writer.writeheader()

    failed_fields =  ['qname', 'qcut', 'seq1', 'qual1', 'seq2', 'qual2', 'prop_N', 'mseq']
    failed_writer = DictWriter(failed_csv, failed_fields)
    failed_writer.writeheader()

    aligned = {}
    if nthreads:
        from multiprocessing import Pool
        pool = Pool(processes=nthreads)
        iter = pool.imap(parse_sam, iterable=matchmaker(remap_csv), chunksize=100)
    else:
        iter = itertools.imap(parse_sam, matchmaker(remap_csv))

    for rname, mseqs, insert_list, failed_list in iter:
        if mseqs is None:
            continue

        if rname not in aligned:
            aligned.update({rname: dict([(qcut, {}) for qcut in sam2aln_q_cutoffs])})

        for qcut, mseq in mseqs.iteritems():
            # collect identical merged sequences
            if mseq not in aligned[rname][qcut]:
                aligned[rname][qcut].update({mseq: 0})
            aligned[rname][qcut][mseq] += 1

            # write out inserts to CSV
            for items in insert_list:
                insert_writer.writerow(dict(zip(insert_fields, items)))

            # write out failed read mergers to CSV
            for items in failed_list:
                failed_writer.writerow(dict(zip(failed_fields, items)))

    failed_csv.close()
    insert_csv.close()

    # write out merged sequences to file
    aligned_fields = ['refname', 'qcut', 'rank', 'count', 'offset', 'seq']
    aligned_writer = DictWriter(aligned_csv, aligned_fields)
    aligned_writer.writeheader()
    for rname, data in aligned.iteritems():
        for qcut, data2 in data.iteritems():
            # sort variants by count
            intermed = [(count, len_gap_prefix(s), s) for s, count in data2.iteritems()]
            intermed.sort(reverse=True)
            for rank, (count, offset, seq) in enumerate(intermed):
                aligned_writer.writerow(dict(refname=rname,
                                             qcut=qcut,
                                             rank=rank,
                                             count=count,
                                             offset=offset,
                                             seq=seq.strip('-')))
    aligned_csv.close()



def main():
    args = parseArgs()
    sam2aln(remap_csv=args.remap_csv,
            aligned_csv=args.aligned_csv,
            insert_csv=args.insert_csv,
            failed_csv=args.failed_csv,
            nthreads=args.p)

if __name__ == '__main__':
    main()
