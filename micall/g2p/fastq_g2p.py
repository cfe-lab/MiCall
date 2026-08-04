#! /usr/bin/env python3.6

import argparse
import contextlib
from collections import Counter
import csv
from csv import DictWriter
from itertools import takewhile
import os
import re

# noinspection PyUnresolvedReferences
from gotoh import align_it, align_it_aa

from micall.core.consensus_builder import ConsensusBuilder
from micall.core.sam2aln import merge_pairs, SAM2ALN_Q_CUTOFFS
from micall.utils.big_counter import BigCounter
from micall.utils.translation import translate, reverse_and_complement
from micall.core.project_config import ProjectConfig, G2P_SEED_NAME

# screen for in-frame deletions
pat = re.compile('([A-Z])(---)+([A-Z])')

Q_CUTOFF = SAM2ALN_Q_CUTOFFS[0]
DEFAULT_MIN_COUNT = 3
GAP_OPEN_COST = 10
GAP_EXTEND_COST = 3
USE_TERMINAL_COST = 1
MIN_PAIR_ALIGNMENT_SCORE = 20
MIN_VALID = 7500
MIN_VALID_PERCENT = 75.0
COORDINATE_REF_NAME = "V3LOOP"
# PSSM_AREF = 'CTRPNXNNTXXRKSIRIXXXGPGQXXXAFYATXXXXGDIIGDIXXRQAHC'.replace('X', '')


def parse_args():
    parser = argparse.ArgumentParser(description='Calculate g2p scores from amino acid sequences.')

    parser.add_argument('fastq1', type=argparse.FileType('r'),
                        help='<input> FASTQ file containing read 1 reads')
    parser.add_argument('fastq2', type=argparse.FileType('r'),
                        help='<input> FASTQ file containing read 2 reads')
    parser.add_argument('g2p_csv', type=argparse.FileType('w'),
                        help='<output> CSV containing g2p predictions.')
    parser.add_argument('g2p_summary_csv', type=argparse.FileType('w'),
                        help='<output> CSV containing overall call for the sample.')
    parser.add_argument('unmapped1',
                        type=argparse.FileType('w'),
                        help='<output> FASTQ file containing unmapped read 1 reads')
    parser.add_argument('unmapped2',
                        type=argparse.FileType('w'),
                        help='<output> FASTQ file containing unmapped read 2 reads')
    parser.add_argument('aligned_csv',
                        type=argparse.FileType('w'),
                        help='<output> CSV containing mapped reads aligned to HIV seed')
    parser.add_argument('merged_contigs_csv',
                        type=argparse.FileType('w'),
                        help='<output> CSV containing amplicon contigs based on length')
    return parser.parse_args()


def fastq_g2p(pssm,
              fastq1,
              fastq2,
              g2p_csv,
              g2p_summary_csv=None,
              unmapped1=None,
              unmapped2=None,
              aligned_csv=None,
              min_count=1,
              min_valid=1,
              min_valid_percent=0.0,
              merged_contigs_csv=None):
    g2p_filename = getattr(g2p_csv, 'name', None)
    if g2p_filename is None:
        count_prefix = None
    else:
        working_path = os.path.dirname(g2p_csv.name)
        count_prefix = os.path.join(working_path, 'read_counts')
    project_config = ProjectConfig.loadDefault()
    hiv_seed = project_config.getReference(G2P_SEED_NAME)
    coordinate_ref = project_config.getReference(COORDINATE_REF_NAME)
    v3loop_ref = extract_target(hiv_seed, coordinate_ref)
    reader = FastqReader(fastq1, fastq2)
    merged_reads = merge_reads(reader)
    consensus_builder = ConsensusBuilder()
    counted_reads = consensus_builder.build(merged_reads)
    trimmed_reads = trim_reads(counted_reads, v3loop_ref)
    mapped_reads = write_unmapped_reads(trimmed_reads, unmapped1, unmapped2)
    read_counts = count_reads(mapped_reads, count_prefix)
    if aligned_csv is not None:
        read_counts = write_aligned_reads(read_counts,
                                          aligned_csv,
                                          hiv_seed,
                                          v3loop_ref)

    write_rows(pssm,
               read_counts,
               g2p_csv,
               g2p_summary_csv,
               min_count,
               min_valid=min_valid,
               min_valid_percent=min_valid_percent)
    if merged_contigs_csv is not None:
        contig_writer = DictWriter(merged_contigs_csv, ['contig'])
        contig_writer.writeheader()
        for consensus in consensus_builder.get_consensus_by_lengths():
            unambiguous_consensus = consensus.replace('N', '').replace('-', '')
            if unambiguous_consensus:
                contig_writer.writerow(dict(contig=consensus))


def write_rows(pssm,
               read_counts,
               g2p_csv,
               g2p_summary_csv=None,
               min_count=1,
               min_valid=1,
               min_valid_percent=0.0):
    """ Apply G2P algorithm to merged reads.

    :param pssm: PSSM library to calculate scores.
    :param read_counts: aligned reads and counts - [(ref, seq), count]
    :param g2p_csv: open CSV file to write results to
    :param g2p_summary_csv: open CSV file to write summary results to
    :param min_count: minimum count for a sequence to include it in the results
    :param min_valid: minimum valid reads to make the R5/X4 call in the summary
    :param min_valid_percent: minimum percentage of valid reads out of all
        reads to make the R5/X4 call in the summary
    """
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
    top_reads, skip_count = get_top_counts(read_counts, min_count)
    counts = Counter()
    for (ref, s), count in top_reads:
        seq = s.replace('-', '')

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
            valid_pct_display = '0.00'
        else:
            x4_pct = 100.0 * counts['x4'] / counts['valid']
            valid_pct = 100.0 * counts['valid'] / counts['mapped']
            if counts['valid'] < min_valid or valid_pct < min_valid_percent:
                final_call = ''
            elif x4_pct >= 2.0:
                final_call = 'X4'
            else:
                final_call = 'R5'
            x4_pct_display = '{:0.2f}'.format(x4_pct)
            valid_pct_display = '{:0.2f}'.format(valid_pct)
        summary_writer = csv.writer(g2p_summary_csv, lineterminator=os.linesep)
        summary_writer.writerow(
            ['mapped', 'valid', 'X4calls', 'X4pct', 'final', 'validpct'])
        summary_writer.writerow([counts['mapped'],
                                 counts['valid'],
                                 counts['x4'],
                                 x4_pct_display,
                                 final_call,
                                 valid_pct_display])


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
    if fpr >= 3.5:
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
        """ Creates an iterator over the reads in a FASTQ file.

        Iterator items:
        (pair_name, (read1_name, bases, quality), (read2_name, bases, quality))
        """
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

    @staticmethod
    def get_reads(fastq):
        for header, bases, sep, quality in zip(fastq, fastq, fastq, fastq):
            assert header.startswith('@'), header
            assert sep.startswith('+'), sep
            pair_name, read_name = header[1:].rstrip().split(maxsplit=1)
            yield pair_name, (read_name, bases.rstrip('\n'), quality.rstrip('\n'))


def extract_target(seed_ref, coordinate_ref):
    """ Extract a portion of the seed that aligns with the coordinate reference.

    :param seed_ref: seed reference (nucleotide sequence)
    :param coordinate_ref: coordinate reference (amino acid sequence)
    :return: subsequence of seed_ref that maps to coordinate_ref
    """
    best_alignment = (-1000000, '', '', 0)
    for frame_index in range(3):
        seed_aminos = translate('-'*frame_index + seed_ref)
        aseed, acoord, score = align_it_aa(seed_aminos,
                                           coordinate_ref,
                                           GAP_OPEN_COST,
                                           GAP_EXTEND_COST,
                                           USE_TERMINAL_COST)
        best_alignment = max(best_alignment, (score, aseed, acoord, frame_index))
    score, aseed, acoord, frame_index = best_alignment
    assert score >= len(coordinate_ref) // 2, score

    target = []
    seed_index = -frame_index
    for s, c in zip(aseed, acoord):
        if s == '-':
            continue
        seed_index += 3
        if c == '-':
            continue
        target.append(seed_ref[seed_index-3:seed_index])
    return ''.join(target)


def merge_reads(reads):
    """ Generator over merged reads.

    :param reads: iterable of reads from FastqReader
    :return: a generator with items (merged_bases may be None if merge fails):
    (pair_name,
     (read1_name, bases, quality),
     (read2_name, bases, quality),
     merged_bases)
    """
    for pair_name, (r1_name, seq1, qual1), (r2_name, seq2, qual2) in reads:
        if not (seq1 and seq2):
            score = -1
            aligned1 = aligned2 = None
        else:
            seq2_rev = reverse_and_complement(seq2)
            aligned1, aligned2, score = align_it(seq1,
                                                 seq2_rev,
                                                 GAP_OPEN_COST,
                                                 GAP_EXTEND_COST,
                                                 USE_TERMINAL_COST)
        if score >= MIN_PAIR_ALIGNMENT_SCORE and aligned1[0] != '-':
            aligned_qual1 = align_quality(aligned1, qual1)
            aligned_qual2 = align_quality(aligned2, reversed(qual2))
            merged = merge_pairs(aligned1,
                                 aligned2,
                                 aligned_qual1,
                                 aligned_qual2,
                                 q_cutoff=Q_CUTOFF)
        else:
            merged = None
        yield (pair_name,
               (r1_name, seq1, qual1),
               (r2_name, seq2, qual2),
               merged)


def align_quality(nucs, qual):
    qual_iter = iter(qual)
    qual = ''.join('!' if nuc == '-' else next(qual_iter)
                   for nuc in nucs)
    return qual


def trim_reads(reads, v3loop_ref, score_counts=None):
    """ Generator over reads that are aligned to the reference and trimmed.

    :param reads: generator from merge_reads()
    :param v3loop_ref: nucleotide sequence for V3LOOP
    :param score_counts: {score: count} to report on the alignment score
        distribution
    :return: Generator items (aligned_ref and aligned_seq may be None if merge
    or trim fails):
    (pair_name,
     (read1_name, bases, quality),
     (read2_name, bases, quality),
     (aligned_ref, aligned_seq))
    """
    # Measured as roughly halfway between HCV reads and V3LOOP reads
    min_v3_alignment_score = 2*len(v3loop_ref)

    for pair_name, read1, read2, seq in reads:
        trimmed_aligned_ref = trimmed_aligned_seq = None
        if seq is not None:
            aligned_ref, aligned_seq, score = align_it(v3loop_ref,
                                                       seq,
                                                       GAP_OPEN_COST,
                                                       GAP_EXTEND_COST,
                                                       USE_TERMINAL_COST)
            if score_counts is not None:
                score_counts[score] += 1
            if score >= min_v3_alignment_score:
                left_padding = right_padding = 0
                for left_padding, nuc in enumerate(aligned_ref):
                    if nuc != '-':
                        break
                for right_padding, nuc in enumerate(reversed(aligned_ref)):
                    if nuc != '-':
                        break
                start, end = left_padding, -right_padding or None
                trimmed_aligned_ref = aligned_ref[start:end]
                trimmed_aligned_seq = aligned_seq[start:end]
        yield pair_name, read1, read2, (trimmed_aligned_ref,
                                        trimmed_aligned_seq)


def write_unmapped_reads(reads, unmapped1, unmapped2):
    """ Write reads that failed to merge or align with V3LOOP reference.

    :param reads: a generator with these items (aligned_ref and aligned_seq
     are None if it failed to merge or align)
    (pair_name,
     (read1_name, bases, quality),
     (read2_name, bases, quality),
     (aligned_ref, aligned_seq))
    :param unmapped1: open FASTQ file that the failed reads will be written to
    :param unmapped2: open FASTQ file that the failed reads will be written to
    :return: a generator with (aligned_ref, aligned_seq) for the reads that
    didn't fail.
    """
    for pair_name, read1, read2, (aligned_ref, aligned_seq) in reads:
        if aligned_ref is not None and aligned_seq is not None:
            yield aligned_ref, aligned_seq
        elif unmapped1 is not None and unmapped2 is not None:
            write_fastq_read(unmapped1, pair_name, read1)
            write_fastq_read(unmapped2, pair_name, read2)


def write_fastq_read(fastq, pair_name, read):
    read_name, bases, quality = read
    fastq.write("""\
@{} {}
{}
+
{}
""".format(pair_name, read_name, bases, quality))


def count_reads(reads, file_prefix):
    """ Count unique sequences in trimmed reads.

    :param reads: a generator with items:
    (aligned_ref, aligned_seq)
    :param file_prefix: used to store temp files for counting sequences
    :return: generator of ((aligned_ref, aligned_seq), count)
    """
    if file_prefix is None:
        all_counts = Counter()
        counts_context = contextlib.suppress()  # Dummy value.
    else:
        all_counts = counts_context = BigCounter(file_prefix)
    with counts_context:
        for aligned_ref, aligned_seq in reads:
            key = aligned_ref + '\t' + aligned_seq
            all_counts[key] += 1
        for key, count in all_counts.items():
            yield tuple(key.split('\t')), count


def get_top_counts(read_counts, min_count=1):
    """ Report most common reads, in descending order.

    :param read_counts: a generator with items ((aligned_ref, aligned_seq), count)
    :param min_count: the minimum count to be included in the result
    :return: [((aligned_ref, aligned_seq), count)], ignored_count where the
        items are by descending count, and the ignored_count is the total of
        all counts less than min_count
    """
    counts = Counter()
    ignored_count = 0
    for read, count in read_counts:
        if count < min_count:
            ignored_count += count
        else:
            counts[read] += count
    return counts.most_common(), ignored_count


def write_aligned_reads(counts, aligned_csv, hiv_seed, v3loop_ref):
    """ Write reads, aligned to the HIV seed sequence.

    :param counts: [((aligned_ref, aligned_seq), count)]
    :param aligned_csv: open CSV file to write the aligned reads to
    :param hiv_seed: seed reference to align the V3LOOP ref to
    :param v3loop_ref: reference the reads were all aligned to
    """
    writer = csv.DictWriter(
        aligned_csv,
        ['refname', 'qcut', 'rank', 'count', 'offset', 'seq'],
        lineterminator=os.linesep)
    writer.writeheader()

    seed_vs_v3, v3_vs_seed, score = align_it(hiv_seed,
                                             v3loop_ref,
                                             GAP_OPEN_COST,
                                             GAP_EXTEND_COST,
                                             USE_TERMINAL_COST)

    # Count dashes at start of aligned_v3loop
    v3_offset = sum(1 for _ in takewhile(lambda c: c == '-', v3_vs_seed))
    seed_positions = list(zip(seed_vs_v3[v3_offset:], v3_vs_seed[v3_offset:]))

    for rank, read_count in enumerate(counts):
        yield read_count
        (v3_vs_read, read_vs_v3), count = read_count
        is_started = False
        seq_offset = 0
        seq = ''
        read_positions = iter(zip(v3_vs_read, read_vs_v3))
        for seed_char, v3_vs_seed_char in seed_positions:
            if v3_vs_seed_char == '-':
                seq += '-'
                continue
            try:
                while True:
                    v3_vs_read_char, read_char = next(read_positions)
                    if v3_vs_read_char != '-':
                        break
            except StopIteration:
                break
            if seed_char != '-':
                if read_char == '-' and not is_started:
                    seq_offset += 1
                else:
                    is_started = True
                    seq += read_char
        seq = seq.rstrip('-')
        writer.writerow(dict(refname=G2P_SEED_NAME,
                             qcut=Q_CUTOFF,
                             rank=rank,
                             count=count,
                             offset=v3_offset + seq_offset,
                             seq=seq))


def main():
    args = parse_args()
    from micall.g2p.pssm_lib import Pssm
    pssm = Pssm()
    fastq_g2p(pssm=pssm,
              fastq1=args.fastq1,
              fastq2=args.fastq2,
              g2p_csv=args.g2p_csv,
              g2p_summary_csv=args.g2p_summary_csv,
              unmapped1=args.unmapped1,
              unmapped2=args.unmapped2,
              aligned_csv=args.aligned_csv,
              min_count=DEFAULT_MIN_COUNT,
              min_valid=MIN_VALID,
              min_valid_percent=MIN_VALID_PERCENT,
              merged_contigs_csv=args.merged_contigs_csv)


if __name__ == '__main__':
    # note, must be called from project root if executing directly
    # i.e., python micall/g2p/fastq_g2p.py -h
    main()
