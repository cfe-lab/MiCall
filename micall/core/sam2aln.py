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
    settings.py
"""

import argparse
import itertools
import re

from micall.settings import max_prop_N, read_mapping_cutoff, sam2aln_q_cutoffs
from csv import DictReader, DictWriter
import os
import collections

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
            insertions[left] = (seq[left:(left+length)], qual[left:(left+length)])
            left += length
            continue
        # Soft clipping leaves the sequence in the SAM - so we should skip it
        elif token[-1] == 'S':
            left += length
            continue
        else:
            raise RuntimeError(
                "Unable to handle CIGAR token: {} - quitting".format(token))

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
        old_row = cached_rows.pop(qname, None)
        if old_row is None:
            cached_rows[qname] = row
        else:
            # current row should be the second read of the pair
            yield old_row, row


def parse_sam(rows):
    """ Merge two matched reads into a single aligned read.
    
    Also report insertions and failed merges.
    @param rows: tuple holding a pair of matched rows - forward and reverse reads
    @return: (refname, merged_seqs, insert_list, failed_list) where
        merged_seqs is {qcut: seq} the merged sequence for each cutoff level
        insert_list is [{'qname': query_name,
                         'fwd_rev': 'F' or 'R',
                         'refname': refname,
                         'pos': pos,
                         'insert': insertion_sequence,
                         'qual': insertion_quality_sequence}] insertions
        relative to the reference sequence.
        failed_list is [{'qname': query_name,
                         'qcut': qcut,
                         'seq1': seq1,
                         'qual1': qual1,
                         'seq2': seq2,
                         'qual2': qual2,
                         'prop_N': proportion_of_Ns,
                         'mseq': merged_sequence}] sequences that failed to
        merge.
    """
    row1, row2 = rows
    mseqs = {}
    failed_list = []
    insert_list = []
    rname = row1['rname']
    qname = row1['qname']

    cigar1 = row1['cigar']
    cigar2 = row2['cigar']
    if (cigar1 == '*' or
        int(row1['mapq']) < read_mapping_cutoff or
        cigar2 == '*' or
        int(row2['mapq']) < read_mapping_cutoff or
        row1['rname'] != row2['rname']):
        
        failed_list.append({'qname': qname,
                            'seq1': row1['seq'],
                            'qual1': row1['qual'],
                            'seq2': row2['seq'],
                            'qual2': row2['qual'],
                            'mapq1': row1['mapq'],
                            'mapq2': row2['mapq']})
    else:
        pos1 = int(row1['pos'])-1  # convert 1-index to 0-index
        _, seq1, qual1, inserts = apply_cigar(cigar1, row1['seq'], row1['qual'])
    
        # report insertions relative to sample consensus
        for left, (iseq, iqual) in inserts.iteritems():
            insert_list.append({'qname': qname,
                                'fwd_rev': 'F' if is_first_read(row1['flag']) else 'R',
                                'refname': rname,
                                'pos': pos1+left,
                                'insert': iseq,
                                'qual': iqual})
    
        seq1 = '-'*pos1 + seq1  # pad sequence on left
        qual1 = '!'*pos1 + qual1  # assign lowest quality to gap prefix so it does not override mate
    
    
        # now process the mate
        pos2 = int(row2['pos'])-1  # convert 1-index to 0-index
        _, seq2, qual2, inserts = apply_cigar(cigar2, row2['seq'], row2['qual'])
        for left, (iseq, iqual) in inserts.iteritems():
            insert_list.append({'qname': qname,
                                'fwd_rev': 'F' if is_first_read(row2['flag']) else 'R',
                                'refname': rname,
                                'pos': pos2+left,
                                'insert': iseq,
                                'qual': iqual})
        seq2 = '-'*pos2 + seq2
        qual2 = '!'*pos2 + qual2
    
        # merge reads
        for qcut in sam2aln_q_cutoffs:
            mseq = merge_pairs(seq1, seq2, qual1, qual2, qcut)
            prop_N = mseq.count('N') / float(len(mseq.strip('-')))
            if prop_N > max_prop_N:
                # fail read pair
                failed_list.append({'qname': qname,
                                    'qcut': qcut,
                                    'seq1': seq1,
                                    'qual1': qual1,
                                    'seq2': seq2,
                                    'qual2': qual2,
                                    'prop_N': prop_N,
                                    'mseq': mseq})
                continue
            mseqs[qcut] = mseq

    return rname, mseqs, insert_list, failed_list

def parse_sam_in_threads(remap_csv, nthreads):
    """ Call parse_sam() in multiple processes.
    
    Launch a multiprocessing pool, walk through the iterator, and then be sure
    to close the pool at the end.
    """
    from multiprocessing import Pool
    pool = Pool(processes=nthreads)
    try:
        reads = pool.imap(parse_sam, iterable=matchmaker(remap_csv), chunksize=100)
        for read in reads:
            yield read
    finally:
        pool.close()
        pool.join()

def sam2aln(remap_csv, aligned_csv, insert_csv, failed_csv, nthreads=None):
    # prepare outputs
    insert_fields =  ['qname', 'fwd_rev', 'refname', 'pos', 'insert', 'qual']
    insert_writer = DictWriter(insert_csv, insert_fields, lineterminator=os.linesep)
    insert_writer.writeheader()

    failed_fields =  ['qname', 'qcut', 'seq1', 'qual1', 'seq2', 'qual2', 'prop_N', 'mseq', 'mapq1', 'mapq2']
    failed_writer = DictWriter(failed_csv, failed_fields, lineterminator=os.linesep)
    failed_writer.writeheader()

    empty_region = collections.defaultdict(collections.Counter)
    aligned = collections.defaultdict(empty_region.copy)
    if nthreads:
        iter = parse_sam_in_threads(remap_csv, nthreads)
    else:
        iter = itertools.imap(parse_sam, matchmaker(remap_csv))

    for rname, mseqs, insert_list, failed_list in iter:
        region = aligned[rname]

        for qcut, mseq in mseqs.iteritems():
            # collect identical merged sequences
            mseq_counter = region[qcut]
            mseq_counter[mseq] += 1

        # write out inserts to CSV
        insert_writer.writerows(insert_list)

        # write out failed read mergers to CSV
        failed_writer.writerows(failed_list)

    # write out merged sequences to file
    aligned_fields = ['refname', 'qcut', 'rank', 'count', 'offset', 'seq']
    aligned_writer = DictWriter(aligned_csv, aligned_fields, lineterminator=os.linesep)
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



def main():
    args = parseArgs()
    sam2aln(remap_csv=args.remap_csv,
            aligned_csv=args.aligned_csv,
            insert_csv=args.insert_csv,
            failed_csv=args.failed_csv,
            nthreads=args.p)

if __name__ == '__main__':
    main()
