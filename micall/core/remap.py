"""
Takes preliminary SAM as CSV input.  Iterative re-mapping of reads from
original FASTQ files.
Also report the number of reads mapped before and after processing.
See docs/design/remap.md for a detailed description.
"""

import argparse
import typing
from typing import Callable, Optional
from collections import Counter, defaultdict
import csv
from csv import DictReader, DictWriter
from functools import partial
from logging import getLogger
import os
import re
import shutil
import sys

# noinspection PyUnresolvedReferences
from pathlib import Path

from gotoh import align_it

import Levenshtein

from micall.core import project_config
from micall.core.project_config import ProjectConfig
from micall.core.sam2aln import apply_cigar, merge_pairs, merge_inserts
from micall.core.prelim_map import READ_GAP_OPEN, READ_GAP_EXTEND, REF_GAP_OPEN, \
    REF_GAP_EXTEND, check_fastq
from micall.utils.externals import Bowtie2, Bowtie2Build, LineCounter
from micall.utils.translation import reverse_and_complement
from micall.utils.work_dir import WorkDir
from micall.utils.stderr import Stderr
from micall.utils.remap_callback import RemapCallback
from micall.utils.cache import cached

CONSENSUS_Q_CUTOFF = 20         # Min Q for base to contribute to conseq (pileup2conseq)
MIN_MAPPING_EFFICIENCY = 0.95   # Fraction of fastq reads mapped needed
MAX_REMAPS = 3                  # Number of remapping attempts if mapping efficiency unsatisfied
SAM_FLAG_IS_UNMAPPED = 0x4
SAM_FLAG_IS_MATE_UNMAPPED = 0x8
SAM_FLAG_IS_FIRST_SEGMENT = 0x40
PARTIAL_CONTIG_SUFFIX = 'partial'
REVERSED_CONTIG_SUFFIX = 'reversed'
EXCLUDED_CONTIG_SUFFIX = 'excluded'
ARE_CONTIGS_MERGED = False

# SAM file format
SAM_FIELDS = [
    'qname',
    'flag',
    'rname',
    'pos',
    'mapq',
    'cigar',
    'rnext',
    'pnext',
    'tlen',
    'seq',
    'qual'
]

cigar_re = re.compile('[0-9]+[MIDNSHPX=]')  # CIGAR token

logger = getLogger(__name__)
indel_re = re.compile('[+-][0-9]+')
line_counter = LineCounter()


def is_first_read(flag):
    """
    Interpret bitwise flag from SAM field.
    Returns True or False indicating whether the read is the first read in a pair.
    """
    return (int(flag) & SAM_FLAG_IS_FIRST_SEGMENT) != 0


def is_unmapped_read(flag):
    """
    Interpret bitwise flag from SAM field.
    Returns True if the read is unmapped.
    """
    return (int(flag) & SAM_FLAG_IS_UNMAPPED) != 0


def is_short_read(read_row, max_primer_length):
    """
    Return length of matched intervals in read
    :param read_row:
    :param max_primer_length:
    :return:
    """
    cigar = read_row['cigar']
    sizes = re.findall(r'(\d+)M', cigar)
    match_length = max(map(int, sizes))
    return match_length <= max_primer_length


def merge_reads(quality_cutoff, read_pair):
    """ Merge a pair of reads.

    Also skip reads that don't meet certain criteria.
    @param quality_cutoff: minimum quality score for a base to be counted
    @param read_pair: a sequence of two sequences, each with fields from a
    SAM file record
    @return: (rname, mseq, merged_inserts, qual1, qual2) or None to skip the pair
    """
    read1, read2 = read_pair
    if read2 and read1[2] != read2[2]:
        # region mismatch, ignore the read pair.
        return None
    filtered_reads = []
    rname = None
    for read in read_pair:
        if not read:
            continue
        (_qname,
         flag,
         rname,
         refpos_str,
         _mapq,
         cigar,
         _rnext,
         _pnext,
         _tlen,
         seq,
         qual) = read[:11]  # ignore optional fields
        if is_unmapped_read(flag):
            continue
        filtered_reads.append(dict(rname=rname,
                                   cigar=cigar,
                                   seq=seq,
                                   qual=qual,
                                   pos=int(refpos_str)))
    if not filtered_reads:
        return None
    seq1, qual1, ins1 = apply_cigar(filtered_reads[0]['cigar'],
                                    filtered_reads[0]['seq'],
                                    filtered_reads[0]['qual'],
                                    filtered_reads[0]['pos']-1)
    if len(filtered_reads) == 1:
        seq2 = qual2 = ''
        ins2 = None
    else:
        seq2, qual2, ins2 = apply_cigar(filtered_reads[1]['cigar'],
                                        filtered_reads[1]['seq'],
                                        filtered_reads[1]['qual'],
                                        filtered_reads[1]['pos']-1)
    mseq = merge_pairs(seq1, seq2, qual1, qual2, q_cutoff=quality_cutoff)
    merged_inserts = merge_inserts(ins1, ins2, quality_cutoff)
    return rname, mseq, merged_inserts, qual1, qual2


def extract_relevant_seed(aligned_conseq, aligned_seed):
    """ Extract the portion of a seed that is relevant to the consensus.

    :param str aligned_conseq: the consensus sequence, aligned to the seed
    :param str aligned_seed: the seed reference, aligned to the consensus
    :return: the portion of the seed that mapped to or was surrounded by the
    consensus.
    """
    match = re.match('-*([^-](.*[^-])?)', aligned_conseq)
    return aligned_seed[match.start(1):match.end(1)].replace('-', '')


def sam_to_conseqs(samfile,
                   quality_cutoff=0,
                   debug_reports=None,
                   seeds=None,
                   is_filtered=False,
                   filter_coverage=1,
                   distance_report=None,
                   original_seeds=None):
    """ Build consensus sequences for each reference from a SAM file.

    @param samfile: an open file in the SAM format containing reads with their
        mapped position and quality scores
    @param quality_cutoff: minimum quality score for a base to be counted
    @param debug_reports: {(rname, pos): None} a dictionary with keys for all
        of the regions and positions that you want a report for. The value
        will be set to a string describing the counts and qualities at that
        position.
    @param seeds: {name: sequence} If this is set,
        any positions without coverage will be set to the base from the seed
        reference. If there are no reads mapped to a reference, it will not
        be included as a new consensus.
    @param is_filtered: if True, then any consensus that has migrated so far
        from its seed that it is closer to a different seed, will not be
        included as a new consensus.
    @param filter_coverage: when filtering on consensus distance, only include
        portions with at least this depth of coverage
    @param distance_report: empty dictionary or None. Dictionary will return:
        {rname: {'seed_dist': seed_dist, 'other_dist': other_dist,
        'other_seed': other_seed}}
    @param original_seeds: {name: sequence} Original seed references used in
        the distance report.
    @return: {reference_name: consensus_sequence}
    """

    if debug_reports:
        for key in debug_reports.keys():
            debug_reports[key] = Counter()

    # refmap structure: {refname: {pos: Counter({nuc: count})}}
    refmap = {}

    pairs = matchmaker(samfile, include_singles=True)
    merged_reads = map(partial(merge_reads, quality_cutoff), pairs)
    read_counts = Counter()
    for merged_read in merged_reads:
        if merged_read is None:
            continue
        rname, mseq, merged_inserts, qual1, qual2 = merged_read
        read_counts[rname] += 1
        pos_nucs = refmap.get(rname)
        if pos_nucs is None:
            pos_nucs = refmap[rname] = defaultdict(Counter)
        update_counts(rname,
                      qual1,
                      qual2,
                      mseq,
                      merged_inserts,
                      pos_nucs,
                      debug_reports)

    if debug_reports:
        for key, counts in debug_reports.items():
            mixtures = []
            nucs = set()
            qualities = set()
            for nuc, quality in counts.keys():
                nucs.add(nuc)
                qualities.add(quality)
            qualities = sorted(qualities)
            for min_quality in qualities:
                filtered_counts = Counter()
                for (nuc, nuc_qual), count in counts.items():
                    if nuc_qual >= min_quality:
                        filtered_counts[nuc] += count
                mixture = []
                for nuc, count in filtered_counts.items():
                    mixture.append('{}: {}'.format(nuc, count))
                mixtures.append('{}{{{}}}'.format(min_quality,
                                                  ', '.join(sorted(mixture))))
            debug_reports[key] = ', '.join(sorted(mixtures))

    new_conseqs = counts_to_conseqs(refmap, seeds)
    is_filtering = original_seeds and is_filtered
    if not is_filtering or len(new_conseqs) <= 1:
        return new_conseqs

    relevant_conseqs = {}
    for name in sorted(new_conseqs.keys()):
        conseq = new_conseqs[name]
        counts = refmap[name]
        relevant_conseq = u''
        for pos, c in enumerate(conseq, 1):
            pos_counts = sum(counts[pos].values())
            if pos_counts >= filter_coverage:
                relevant_conseq += c
        relevant_conseqs[name] = relevant_conseq

    drop_drifters(relevant_conseqs, original_seeds, distance_report, read_counts)
    return {seed_name: new_conseqs[seed_name] for seed_name in relevant_conseqs}


def drop_drifters(relevant_conseqs: typing.Dict[str, str],
                  original_seeds: typing.Dict[str, str],
                  distance_report: typing.Dict[str, dict],
                  read_counts: typing.Dict[str, int]):
    gap_open_penalty = 15
    gap_extend_penalty = 3
    use_terminal_gap_penalty = 1
    is_filtering = True
    while is_filtering and len(relevant_conseqs) > 1:
        drifted_seeds = []  # [(count, name)]
        for name in sorted(relevant_conseqs.keys()):
            relevant_conseq = relevant_conseqs[name]
            if not relevant_conseq:
                # None of the coverage was acceptable.
                drifted_seeds.append((read_counts[name], name))
                continue
            if len(relevant_conseq) > 10_000:
                # Gotoh is too expensive on long sequences, and we only really
                # need to check drift between HCV genotypes.
                continue

            other_seed = other_dist = seed_dist = None
            for seed_name in sorted(relevant_conseqs.keys()):
                seed_ref = original_seeds[seed_name]
                aligned_seed, aligned_conseq, _score = align_it(seed_ref,
                                                                relevant_conseq,
                                                                gap_open_penalty,
                                                                gap_extend_penalty,
                                                                use_terminal_gap_penalty)
                relevant_seed = extract_relevant_seed(aligned_conseq, aligned_seed)
                d = Levenshtein.distance(relevant_seed, relevant_conseq)
                if seed_name == name:
                    seed_dist = d
                elif other_dist is None or d < other_dist:
                    other_seed = seed_name
                    other_dist = d

            if seed_dist > other_dist:
                # Consensus is farther from starting seed than another seed: drop it?
                drifted_seeds.append((read_counts[name], name))
            if distance_report is not None:
                distance_report[name] = dict(seed_dist=seed_dist,
                                             other_dist=other_dist,
                                             other_seed=other_seed)
        distance_report = None  # Only update during first iteration.
        if drifted_seeds:
            drifted_seeds.sort()
            dropped_seed = drifted_seeds[0][1]
            del relevant_conseqs[dropped_seed]
        is_filtering = len(drifted_seeds) > 1


def update_counts(rname,
                  qual1,
                  qual2,
                  mseq,
                  merged_inserts,
                  pos_nucs,
                  debug_reports=None):
    """ Update the counts for each position within a merged read.

    @param rname: the reference name this read mapped to
    @param qual1: the quality scores for the forward read
    @param qual2: the quality scores for the reverse read
    @param mseq: the merged sequence of the forward and reverse reads
    @param merged_inserts: {pos: seq}
    @param pos_nucs: {pos: {nuc: count}}
    @param debug_reports: {(rname, pos): {nuc+qual: count}} a dictionary with
        keys for all of the regions and positions that you want a report for.
    """
    is_started = False
    for pos, nuc in enumerate(mseq, 1):
        if not is_started:
            if nuc == '-':
                continue
            is_started = True
        if nuc != 'n':
            nuc_counts = pos_nucs[pos]
            if nuc == 'N':
                nuc_counts[nuc] = -1
            elif nuc == '-':
                nuc_counts[nuc] = -2
            else:
                ins = merged_inserts.get(pos)
                if ins and len(ins) % 3 == 0:
                    nuc_counts[nuc + ins] += 1
                else:
                    nuc_counts[nuc] += 1
                if debug_reports:
                    counts = debug_reports.get((rname, pos))
                    if counts is not None:
                        q = qual1[pos-1] if pos <= len(qual1) else qual2[pos-1]
                        counts[nuc + q] += 1


def counts_to_conseqs(refmap, seeds=None):
    conseqs = {}
    for refname, pos_nucs in refmap.items():
        if not any((any(n > 0 for n in counts.values())
                    for counts in pos_nucs.values())):
            # Nothing mapped, so no consensus.
            continue
        conseq = ''
        deletion = ''
        seed = seeds and seeds.get(refname)
        end = max(pos_nucs.keys())+1
        if seed:
            end = max(end, len(seed)+1)
        for pos in range(1, end):
            nuc_counts = pos_nucs.get(pos)
            most_common = nuc_counts and find_top_token(nuc_counts) or None
            if most_common is None:
                if seed is None:
                    conseq += 'N'
                else:
                    conseq += seed[pos-1]
            elif most_common == '-':
                if seed is None:
                    deletion += '-'
                else:
                    conseq += seed[pos-1]
            else:
                if deletion:
                    if len(deletion) % 3 != 0:
                        conseq += deletion
                    deletion = ''
                conseq += most_common
        conseqs[refname] = conseq
    return conseqs


def build_conseqs(samfilename,
                  seeds=None,
                  is_filtered=False,
                  filter_coverage=1,
                  distance_report=None,
                  original_seeds=None,
                  debug_reports=None):
    """ Build the new consensus sequences from the mapping results.

    @param samfilename: the mapping results in SAM format
    @param seeds: {name: sequence} If this is set,
        any positions without coverage will be set to the base from the seed
        reference. If there are no reads mapped to a reference, it will not
        be included as a new consensus.
    @param is_filtered: if True, then any consensus that has migrated so far
        from its seed that it is closer to a different seed, will not be
        included as a new consensus.
    @param filter_coverage: when filtering on consensus distance, only include
        portions with at least this depth of coverage
    @param distance_report: empty dictionary or None. Dictionary will return:
        {rname: {'seed_dist': seed_dist, 'other_dist': other_dist,
        'other_seed': other_seed}}
    @param original_seeds: {name: sequence} Original seed references used in
        the distance report.
    @param debug_reports: {(rname, pos): None} a dictionary with keys for all
        of the regions and positions that you want a report for. The value
        will be set to a string describing the counts and qualities at that
        position.
    @return: {reference_name: consensus_sequence}
    """
    with open(samfilename) as samfile:
        conseqs = sam_to_conseqs(samfile,
                                 CONSENSUS_Q_CUTOFF,
                                 seeds=seeds,
                                 is_filtered=is_filtered,
                                 filter_coverage=filter_coverage,
                                 distance_report=distance_report,
                                 original_seeds=original_seeds,
                                 debug_reports=debug_reports)

    return conseqs


def write_remap_counts(remap_counts_writer, counts, title, distance_report=None):
    distance_report = distance_report or {}
    for refname in sorted(counts.keys()):
        row = distance_report.get(refname, {})
        row['type'] = title + ' ' + refname
        row['count'] = counts[refname]
        remap_counts_writer.writerow(row)


@cached("remap",
        parameters=["count_threshold", "rdgopen", "rfgopen", "gzip", "debug_file_prefix"],
        outputs=["remap_csv", "remap_counts_csv", "remap_conseq_csv", "unmapped1", "unmapped2"])
def remap(fastq1: Path,
          fastq2: Path,
          prelim_csv: Path,
          remap_csv: Path,
          remap_counts_csv: Path,
          remap_conseq_csv: Path,
          unmapped1: Path,
          unmapped2: Path,
          count_threshold=10,
          rdgopen=READ_GAP_OPEN,
          rfgopen=REF_GAP_OPEN,
          gzip=False,
          debug_file_prefix: Optional[str] = None,
          ) -> None:
    """
    Iterative re-map reads from raw paired FASTQ files to a reference sequence set that
    is being updated as the consensus of the reads that were mapped to the last set.
    @param fastq1: Path to input R1 FASTQ
    @param fastq2: Path to input R2 FASTQ
    @param prelim_csv: Path to input CSV output from prelim_csv()
    @param remap_csv:  Path to output CSV, contents of bowtie2 SAM output
    @param remap_counts_csv:  Path to output CSV, counts of reads mapped to regions
    @param remap_conseq_csv:  Path to output CSV, sample- and region-specific consensus sequences
                                generated while remapping reads
    @param unmapped1:  Path to output FASTQ containing R1 reads that did not map to any region
    @param unmapped2:  Path to output FASTQ containing R2 reads that did not map to any region
    @param count_threshold:  minimum number of reads that map to a region for it to be remapped
    @param rdgopen: read gap open penalty
    @param rfgopen: reference gap open penalty
    @param gzip: True if the FASTQ files are gzipped
    @param debug_file_prefix: the prefix for the file path to write debug files.
        If not None, this will be used to write a copy of the reference FASTA
        files and the output SAM files.
    """

    # Get work_path, stderr, and callback from dynamic scope
    work_path = WorkDir.get()
    stderr = Stderr.get()
    callback = RemapCallback.get()

    reffile = work_path / 'temp.fasta'
    samfile = work_path / 'temp.sam'

    bowtie2, bowtie2_build = find_bowtie2()

    # check that the inputs exist
    fastq1_str = check_fastq(str(fastq1), gzip)
    fastq2_str = check_fastq(str(fastq2), gzip)

    # retrieve reference sequences used for preliminary mapping
    projects = project_config.ProjectConfig.loadDefault()
    seeds = projects.getAllReferences()

    # record the raw read count
    raw_count = line_counter.count(fastq1_str, gzip=gzip) // 2  # 4 lines per record in FASTQ, paired

    # Open all output files that will be used throughout the function
    with open(remap_counts_csv, 'w') as remap_counts_file, \
         open(unmapped1, 'w') as unmapped1_file, \
         open(unmapped2, 'w') as unmapped2_file:

        remap_counts_writer = csv.DictWriter(
            remap_counts_file,
            'type count filtered_count seed_dist other_dist other_seed'.split(),
            lineterminator=os.linesep)
        remap_counts_writer.writeheader()
        remap_counts_writer.writerow(dict(type='raw', count=raw_count))

        # convert preliminary CSV to SAM, count reads
        with open(samfile, 'w') as sam_f, open(prelim_csv) as prelim_f:
            # transfer filtered counts to map counts for remap loop
            map_counts = convert_prelim(prelim_f,
                                        sam_f,
                                        remap_counts_writer,
                                        count_threshold,
                                        projects)

        # debug_reports = {('SARS-CoV-2-seed', pos): None
        #                  for pos in range(20290, 20311)}
        debug_reports = None
        # regenerate consensus sequences based on preliminary map
        prelim_conseqs = build_conseqs(str(samfile),
                                       seeds=seeds,
                                       debug_reports=debug_reports)
        print_debug_reports(debug_reports, 'Preliminary mapping report:')

        # exclude references with low counts (post filtering)
        conseqs = {rname: prelim_conseqs[rname]
                   for rname in map_counts
                   if rname in prelim_conseqs}

        # start remapping loop
        n_remaps = 0
        new_counts = Counter()
        unmapped_count = raw_count

        while conseqs:
            # reset unmapped files with each iteration
            unmapped1_file.seek(0)
            unmapped1_file.truncate()
            unmapped2_file.seek(0)
            unmapped2_file.truncate()

            if debug_file_prefix is None:
                next_debug_prefix = None
            else:
                next_debug_prefix = '{}_remap{}'.format(debug_file_prefix,
                                                        n_remaps+1)
            unmapped_count = map_to_reference(fastq1_str,
                                              fastq2_str,
                                              conseqs,
                                              str(reffile),
                                              str(samfile),
                                              unmapped1_file,
                                              unmapped2_file,
                                              bowtie2,
                                              bowtie2_build,
                                              raw_count,
                                              rdgopen,
                                              rfgopen,
                                              new_counts,
                                              stderr,
                                              callback,
                                              debug_file_prefix=next_debug_prefix)


            old_seed_names = set(conseqs.keys())
            # regenerate consensus sequences
            distance_report = {}
            conseqs = build_conseqs(str(samfile),
                                    seeds=conseqs,
                                    is_filtered=True,
                                    filter_coverage=count_threshold//2,  # pairs
                                    distance_report=distance_report,
                                    original_seeds=seeds,
                                    debug_reports=debug_reports)
            new_seed_names = set(conseqs.keys())
            n_remaps += 1
            print_debug_reports(debug_reports, f'Remap {n_remaps} report:')

            write_remap_counts(remap_counts_writer,
                               new_counts,
                               title='remap-{}'.format(n_remaps),
                               distance_report=distance_report)

            if new_seed_names == old_seed_names:
                # stopping criterion 1 - none of the regions gained reads
                if all((count <= map_counts[refname])
                       for refname, count in new_counts.items()):
                    break

                # stopping criterion 2 - a sufficient fraction of raw data has been mapped
                mapping_efficiency = sum(new_counts.values()) / float(raw_count)
                if mapping_efficiency > MIN_MAPPING_EFFICIENCY:
                    break

                if n_remaps >= MAX_REMAPS:
                    break

            # deep copy of mapping counts
            map_counts = dict(new_counts)

        # finished iterative phase (out of while loop but still in with block)
        # generate SAM CSV output
        with open(remap_csv, 'w') as remap_csv_file:
            remap_writer = csv.DictWriter(remap_csv_file, SAM_FIELDS, lineterminator=os.linesep)
            remap_writer.writeheader()
            if new_counts:
                splitter = MixedReferenceSplitter(str(work_path))
                split_counts = Counter()
                # At least one read was mapped, so samfile has relevant data
                with open(samfile) as f:
                    for fields in splitter.split(f):
                        remap_writer.writerow(dict(zip(SAM_FIELDS, fields)))
                for rname, (split_file1, split_file2) in splitter.splits.items():
                    refseqs = {rname: conseqs[rname]}
                    unmapped_count += map_to_reference(split_file1.name,
                                                       split_file2.name,
                                                       refseqs,
                                                       str(reffile),
                                                       str(samfile),
                                                       unmapped1_file,
                                                       unmapped2_file,
                                                       bowtie2,
                                                       bowtie2_build,
                                                       raw_count,
                                                       rdgopen,
                                                       rfgopen,
                                                       split_counts,
                                                       stderr,
                                                       callback)
                    new_counts.update(split_counts)
                    with open(samfile, 'r') as f:
                        for fields in splitter.walk(f):
                            remap_writer.writerow(dict(zip(SAM_FIELDS, fields)))

        # write consensus sequences and counts
        with open(remap_conseq_csv, 'w') as remap_conseq_file:
            remap_conseq_file.write('region,sequence\n')  # record consensus sequences for later use
            for refname in new_counts.keys():
                # NOTE this is the consensus sequence to which the reads were mapped, NOT the
                # current consensus!
                conseq = conseqs.get(refname) or projects.getReference(refname)
                remap_conseq_file.write('%s,%s\n' % (refname, conseq))
        write_remap_counts(remap_counts_writer,
                           new_counts,
                           title='remap-final')

        # report number of unmapped reads
        remap_counts_writer.writerow(dict(type='unmapped',
                                          count=unmapped_count))


def print_debug_reports(debug_reports, title):
    if debug_reports is None:
        return
    print(title)
    for (rname, pos), summary in debug_reports.items():
        print(pos, summary)


def map_to_contigs(fastq1,
                   fastq2,
                   contigs_csv,
                   remap_csv,
                   remap_counts_csv,
                   remap_conseq_csv,
                   unmapped1,
                   unmapped2,
                   work_path='',
                   callback=None,
                   rdgopen=READ_GAP_OPEN,
                   rfgopen=REF_GAP_OPEN,
                   stderr=sys.stderr,
                   gzip=False,
                   debug_file_prefix=None,
                   excluded_seeds=None):
    """
    Map reads from raw paired FASTQ files to de novo contigs.
    @param fastq1: input R1 FASTQ
    @param fastq2: input R2 FASTQ
    @param contigs_csv: input CSV output from denovo()
    @param remap_csv:  output CSV, contents of bowtie2 SAM output
    @param remap_counts_csv:  output CSV, counts of reads mapped to regions
    @param remap_conseq_csv:  output CSV, contig sequences used to map reads
    @param unmapped1:  output FASTQ containing R1 reads that did not map to any region
    @param unmapped2:  output FASTQ containing R2 reads that did not map to any region
    @param work_path:  optional path to store working files
    @param callback: a function to report progress with three optional
        parameters - callback(message, progress, max_progress)
    @param rdgopen: read gap open penalty
    @param rfgopen: reference gap open penalty
    @param stderr: an open file object to receive stderr from the bowtie2 calls
    @param gzip: True if the FASTQ files are gzipped
    @param debug_file_prefix: the prefix for the file path to write debug files.
        If not None, this will be used to write a copy of the reference FASTA
        files and the output SAM files.
    @param excluded_seeds: a list of seed names to exclude from mapped reads
    """

    reffile = os.path.join(work_path, 'temp.fasta')
    samfile = os.path.join(work_path, 'temp.sam')

    bowtie2, bowtie2_build = find_bowtie2()

    # check that the inputs exist
    fastq1 = check_fastq(fastq1, gzip)
    fastq2 = check_fastq(fastq2, gzip)

    # record the raw read count
    raw_count = line_counter.count(fastq1, gzip=gzip) // 2  # 4 lines per record in FASTQ, paired

    remap_counts_writer = csv.DictWriter(
        remap_counts_csv,
        'type count filtered_count seed_dist other_dist other_seed'.split(),
        lineterminator=os.linesep)
    remap_counts_writer.writeheader()
    remap_counts_writer.writerow(dict(type='raw', count=raw_count))

    conseqs = read_contigs(contigs_csv, excluded_seeds)

    # start remapping loop
    new_counts = Counter()

    if not conseqs:
        unmapped_count = raw_count
    else:
        unmapped_count = map_to_reference(fastq1,
                                          fastq2,
                                          conseqs,
                                          reffile,
                                          samfile,
                                          unmapped1,
                                          unmapped2,
                                          bowtie2,
                                          bowtie2_build,
                                          raw_count,
                                          rdgopen,
                                          rfgopen,
                                          new_counts,
                                          stderr,
                                          callback,
                                          debug_file_prefix=debug_file_prefix)
        write_remap_counts(remap_counts_writer, new_counts, title='remap')

    # generate SAM CSV output
    remap_writer = csv.DictWriter(remap_csv, SAM_FIELDS, lineterminator=os.linesep)
    remap_writer.writeheader()
    if new_counts:
        splitter = MixedReferenceSplitter(work_path)
        split_counts = Counter()
        # At least one read was mapped, so samfile has relevant data
        with open(samfile) as f:
            for fields in splitter.split(f):
                write_remap_row(remap_writer, fields)
        for rname, (split_file1, split_file2) in splitter.splits.items():
            refseqs = {rname: conseqs[rname]}
            unmapped_count += map_to_reference(split_file1.name,
                                               split_file2.name,
                                               refseqs,
                                               reffile,
                                               samfile,
                                               unmapped1,
                                               unmapped2,
                                               bowtie2,
                                               bowtie2_build,
                                               raw_count,
                                               rdgopen,
                                               rfgopen,
                                               split_counts,
                                               stderr,
                                               callback)
            new_counts.update(split_counts)
            with open(samfile, 'r') as f:
                for fields in splitter.walk(f):
                    write_remap_row(remap_writer, fields)

    write_contig_sequences(conseqs, remap_conseq_csv)
    write_remap_counts(remap_counts_writer,
                       new_counts,
                       title='remap-final')

    # report number of unmapped reads
    remap_counts_writer.writerow(dict(type='unmapped',
                                      count=unmapped_count))


def find_bowtie2():
    bowtie2 = Bowtie2()
    bowtie2_build = Bowtie2Build()
    bowtie2_build.set_logger(logger)
    return bowtie2, bowtie2_build


def write_contig_sequences(conseqs, remap_conseq_csv):
    writer = DictWriter(remap_conseq_csv, ['region', 'sequence'], lineterminator=os.linesep)
    writer.writeheader()
    for contig_name, contig_seq in conseqs.items():
        # NOTE this is the contig sequence to which the reads were mapped, NOT the
        # current consensus!
        if is_reported_region(contig_name):
            writer.writerow(dict(region=contig_name,
                                 sequence=contig_seq))


def write_remap_row(remap_writer, fields):
    row = dict(zip(SAM_FIELDS, fields))
    region = row['rname']
    if is_reported_region(region):
        remap_writer.writerow(row)


def is_reported_region(region):
    return not region.endswith(EXCLUDED_CONTIG_SUFFIX)


def read_contigs(contigs_csv, excluded_seeds=None):
    gap_open_penalty = 15
    gap_extend_penalty = 3
    use_terminal_gap_penalty = 1
    contig_groups = defaultdict(list)  # {group_ref_name: [seq, index, index...]}
    conseqs = {}
    projects = ProjectConfig.loadDefault()
    with contigs_csv:
        contigs_reader = DictReader(contigs_csv)
        for i, row in reversed(list(enumerate(contigs_reader, 1))):
            contig_seq = row['contig']
            match_fraction = float(row['match'])
            is_match = 0.25 <= match_fraction
            is_reversed = match_fraction < 0
            if not (ARE_CONTIGS_MERGED and is_match):
                contig_name = get_contig_name(i,
                                              row['ref'],
                                              is_match,
                                              is_reversed,
                                              excluded_seeds)
                conseqs[contig_name] = contig_seq
                continue
            group_ref_name = row['group_ref']
            contig_group = contig_groups[group_ref_name]
            if not contig_group:
                contig_group.append(projects.getReference(group_ref_name))
            contig_group.append(str(i))
            group_seq = contig_group[0]
            agroup, acontig, score = align_it(group_seq,
                                              contig_seq,
                                              gap_open_penalty,
                                              gap_extend_penalty,
                                              use_terminal_gap_penalty)
            match = re.match('-*([^-](.*[^-])?)', acontig)
            start = match.start(1)
            end = match.end(1)
            merged_seq = agroup[:start] + contig_seq + agroup[end:]
            left_trim = len(agroup) - len(agroup.lstrip('-'))
            right_trim = len(agroup) - len(agroup.rstrip('-'))
            contig_group[0] = merged_seq[left_trim:-right_trim or None]

    is_match = True
    is_reversed = False
    for group_ref_name, contig_group in contig_groups.items():
        (group_seq, *contig_nums) = contig_group
        prefix = '_'.join(reversed(contig_nums))
        contig_name = get_contig_name(prefix,
                                      group_ref_name,
                                      is_match,
                                      is_reversed,
                                      excluded_seeds)
        conseqs[contig_name] = group_seq
    return conseqs


def get_contig_name(prefix, ref_name, is_match, is_reversed, excluded_seeds=None):
    name = '{}-{}'.format(prefix, ref_name)
    if excluded_seeds and ref_name in excluded_seeds:
        name += '-' + EXCLUDED_CONTIG_SUFFIX
    elif is_reversed:
        name += '-' + REVERSED_CONTIG_SUFFIX
    elif not is_match:
        name += '-' + PARTIAL_CONTIG_SUFFIX
    return name


def convert_prelim(prelim_csv,
                   target,
                   remap_counts_writer,
                   count_threshold,
                   projects):
    """ Convert prelim.csv to a SAM file.

    Also count the reads that mapped to each seed, and find the best seed
    for each seed group.
    :param prelim_csv: open CSV file to read from
    :param target: open text file to write SAM version to
    :param remap_counts_writer: open CSV writer for counts
    :param count_threshold: minimum read count to be returned
    :param projects: project definitions
    :return dict: { refname: count }
    """
    conseqs = projects.getAllReferences()
    # write SAM header
    target.write('@HD\tVN:1.0\tSO:unsorted\n')
    for rname in sorted(conseqs.keys()):
        refseq = conseqs[rname]
        target.write('@SQ\tSN:%s\tLN:%d\n' % (rname, len(refseq)))
    target.write('@PG\tID:bowtie2\tPN:bowtie2\tVN:2.2.3\tCL:""\n')
    # iterate through prelim CSV and record counts, transfer rows to SAM
    ref_counts = defaultdict(lambda: [0, 0])  # {rname: [filtered_count, count]}
    reader = csv.DictReader(prelim_csv)
    for row in reader:
        is_unmapped = is_unmapped_read(row['flag'])
        refname = row['rname']
        if is_unmapped:
            refname = '*'
        counts = ref_counts[refname]
        counts[1] += 1  # full count

        # write SAM row
        target.write('\t'.join([row[field] for field in SAM_FIELDS]) + '\n')

        if is_unmapped:
            continue
        if is_short_read(row, max_primer_length=50):
            # exclude short reads
            continue

        counts[0] += 1  # filtered count

    refgroups = {}  # { group_name: (refname, count) }
    for refname, (filtered_count, count) in sorted(ref_counts.items()):
        # report preliminary counts to file
        remap_counts_writer.writerow(
            dict(type='prelim %s' % refname,
                 count=count,
                 filtered_count=filtered_count))
        if refname == '*':
            continue
        refgroup = projects.getSeedGroup(refname)
        _best_ref, best_count = refgroups.get(refgroup,
                                              (None, count_threshold - 1))
        if filtered_count > best_count:
            refgroups[refgroup] = (refname, filtered_count)

    seed_counts = {best_ref: best_count
                   for best_ref, best_count in refgroups.values()}
    return seed_counts


def map_to_reference(fastq1,
                     fastq2,
                     refseqs,
                     reffile,
                     samfile,
                     unmapped1,
                     unmapped2,
                     bowtie2,
                     bowtie2_build,
                     raw_count,
                     rdgopen,
                     rfgopen,
                     new_counts,
                     stderr,
                     callback,
                     debug_file_prefix=None):
    """ Map a pair of FASTQ files to a set of reference sequences.

    @param fastq1: FASTQ file with the forward reads
    @param fastq2: FASTQ file with the reverse reads
    @param refseqs: reference sequences to map against
    @param reffile: file path to use for the reference sequences
    @param samfile: file path to use for the SAM file output from mapping
    @param unmapped1: an open file to write unmapped reads to in FASTQ format
    @param unmapped2: an open file to write unmapped reads to in FASTQ format
    @param bowtie2: a wrapper for calls to bowtie2
    @param bowtie2_build: a wrapper for calls to bowtie2-build
    @param raw_count: the number of reads in fastq1
    @param rdgopen: read gap open penalty
    @param rfgopen: reference gap open penalty
    @param new_counts: a Counter to track how many reads are mapped to each
        reference
    @param stderr: an open file object to receive stderr from the bowtie2 calls
    @param callback: a function to report progress with three optional
        parameters - callback(message, progress, max_progress)
    @param debug_file_prefix: the prefix for the file path to write debug files.
        If not None, this will be used to write a copy of the reference FASTA
        file and the output SAM file.
    """
    # generate reference file from current set of consensus sequences
    outfile = open(reffile, 'w')
    for region, conseq in refseqs.items():
        outfile.write('>%s\n%s\n' % (region, conseq))
    outfile.close()

    # regenerate bowtie2 index files
    bowtie2_build.build(reffile, reffile)

    read_gap_open_penalty = rdgopen
    ref_gap_open_penalty = rfgopen

    # stream output from bowtie2
    bowtie_args = ['--wrapper', 'micall-0',
                   '--quiet',
                   '-x', reffile,
                   '--rdg', "{},{}".format(read_gap_open_penalty,
                                           READ_GAP_EXTEND),
                   '--rfg', "{},{}".format(ref_gap_open_penalty,
                                           REF_GAP_EXTEND),
                   '-1', fastq1,
                   '-2', fastq2,
                   '--no-hd',  # no header lines (start with @)
                   '--local',
                   '-X', '1200',
                   '-p', '1']

    new_counts.clear()
    unmapped_count = 0

    with open(samfile, 'w') as f:
        # write SAM header
        f.write('@HD\tVN:1.0\tSO:unsorted\n')
        for rname, refseq in refseqs.items():
            f.write('@SQ\tSN:%s\tLN:%d\n' % (rname, len(refseq)))
        f.write('@PG\tID:bowtie2\tPN:bowtie2\tVN:2.2.3\tCL:""\n')

        # capture stdout stream to count reads before writing to file
        for i, line in enumerate(bowtie2.yield_output(bowtie_args, stderr=stderr)):
            if callback and i % 1000 == 0:
                callback(progress=i)  # progress monitoring in GUI

            f.write(line)

            items = line.split('\t')
            qname, bitflag, rname, _, _, _, _, _, _, seq, qual = items[:11]

            if is_unmapped_read(bitflag):
                # did not map to any reference
                unmapped_file = unmapped1 if is_first_read(bitflag) else unmapped2
                unmapped_file.write('@%s\n%s\n+\n%s\n' % (qname, seq, qual))
                unmapped_count += 1
                continue

            new_counts[rname] += 1
        if callback:
            callback(progress=raw_count)
    if debug_file_prefix is not None:
        shutil.copy(reffile, debug_file_prefix + '_debug_ref.fasta')
        shutil.copy(samfile, debug_file_prefix + '_debug.sam')
    return unmapped_count


class MixedReferenceSplitter(object):
    def __init__(self, work_path):
        self.work_path = work_path
        self.splits = {}

    # noinspection PyMethodMayBeStatic
    def close_split_file(self, split_file):
        split_file.close()

    @staticmethod
    def walk(sam_lines):
        for line in sam_lines:
            if line.startswith('@'):
                continue
            yield line.strip('\n').split('\t')

    def split(self, sam_lines):
        unmatched = {}
        for fields in self.walk(sam_lines):
            if fields[6] == '=':
                yield fields[:11]
            else:
                flags = int(fields[1])
                if flags & (SAM_FLAG_IS_UNMAPPED | SAM_FLAG_IS_MATE_UNMAPPED):
                    pass
                else:
                    qname = fields[0]
                    match = unmatched.pop(qname, None)
                    if match is None:
                        unmatched[qname] = fields
                    else:
                        mapq = fields[4]
                        match_mapq = match[4]
                        if mapq > match_mapq:
                            rname = fields[2]
                        elif mapq < match_mapq:
                            rname = match[2]
                        else:
                            alignment_score = self.get_alignment_score(fields)
                            match_alignment_score = self.get_alignment_score(match)
                            if alignment_score > match_alignment_score:
                                rname = fields[2]
                            else:
                                rname = match[2]
                        ref_splits = self.splits.get(rname, None)
                        if ref_splits is None:
                            ref_splits = (self.create_split_file(rname, 1),
                                          self.create_split_file(rname, 2))
                            self.splits[rname] = ref_splits
                        fastq1, fastq2 = ref_splits
                        if flags & SAM_FLAG_IS_FIRST_SEGMENT:
                            fwd_read, rev_read = fields, match
                        else:
                            fwd_read, rev_read = match, fields
                        self.write_fastq(fwd_read, fastq1)
                        self.write_fastq(rev_read, fastq2, is_reversed=True)
        for fastq1, fastq2 in self.splits.values():
            self.close_split_file(fastq1)
            self.close_split_file(fastq2)

    @staticmethod
    def get_alignment_score(fields):
        for field in fields[11:]:
            if field.startswith('AS:i:'):
                return int(field[5:])

    def get_split_file_name(self, refname, direction):
        suffix = '_R1.fastq' if direction == 1 else '_R2.fastq'
        return os.path.join(self.work_path, refname + suffix)

    def create_split_file(self, refname, direction):
        split_file_name = self.get_split_file_name(refname, direction)
        return open(split_file_name, 'w')

    @staticmethod
    def write_fastq(fields, fastq, is_reversed=False):
        qname = fields[0]
        seq = fields[9]
        quality = fields[10]
        if is_reversed:
            seq = reverse_and_complement(seq)
            quality = ''.join(reversed(quality))
        fastq.write('@{}\n{}\n+\n{}\n'.format(qname, seq, quality))


def matchmaker(samfile, include_singles=False):
    """
    An iterator that returns pairs of reads sharing a common qname from a SAM file.
    Note that unpaired reads will be left in the cached_rows dictionary and
    discarded.
    @param samfile: open file handle to a SAM file
    @param include_singles: True if unpaired reads should be returned, paired
        with a None value: ([qname, flag, rname, ...], None)
    @return: yields a tuple for each read pair with fields split by tab chars:
        ([qname, flag, rname, ...], [qname, flag, rname, ...])
    """
    ref_names = set()
    cached_rows = {}
    for line in samfile:
        row = line.strip('\n').split('\t')

        if line.startswith('@'):
            if row[0] == '@SQ':
                for field in row[1:]:
                    field_name, value = field.split(':', 1)
                    if field_name == 'SN':
                        ref_names.add(value)
            continue

        qname = row[0]
        ref_name = row[2]
        if ref_name in ref_names:
            old_row = cached_rows.pop(qname, None)
            if old_row is None:
                cached_rows[qname] = row
            else:
                # current row should be the second read of the pair
                yield old_row, row
    if include_singles:
        for row in cached_rows.values():
            yield row, None


def find_top_token(base_counts):
    top_count = top_token = None
    for token, count in base_counts.most_common():
        if top_count is None:
            top_token = token
            top_count = count
        elif count < top_count:
            break
        if token < top_token:
            top_token = token
    if top_token == 'N':
        top_token = None
    return top_token


def main():
    parser = argparse.ArgumentParser(
        description='Iterative remapping of bowtie2 by reference.')

    parser.add_argument('fastq1', help='<input> FASTQ containing forward reads')
    parser.add_argument('fastq2', help='<input> FASTQ containing reverse reads')
    parser.add_argument('contigs_csv',
                        type=argparse.FileType('r'),
                        help='<input> CSV containing assembled contigs')
    parser.add_argument('remap_csv',
                        type=argparse.FileType('w'),
                        help='<output> CSV containing remap output (modified SAM)')
    parser.add_argument('remap_counts_csv',
                        type=argparse.FileType('w'),
                        help='<output> CSV containing numbers of mapped reads')
    parser.add_argument('remap_conseq_csv',
                        type=argparse.FileType('w'),
                        help='<output> CSV containing mapping consensus sequences')
    parser.add_argument('unmapped1',
                        type=argparse.FileType('w'),
                        help='<output> FASTQ R1 of reads that failed to map to any region')
    parser.add_argument('unmapped2',
                        type=argparse.FileType('w'),
                        help='<output> FASTQ R2 of reads that failed to map to any region')
    parser.add_argument("--gzip", help="<optional> FASTQ files are compressed",
                        action='store_true')

    args = parser.parse_args()
    work_path = os.path.dirname(args.remap_csv.name)
    map_to_contigs(
        fastq1=args.fastq1,
        fastq2=args.fastq2,
        contigs_csv=args.contigs_csv,
        remap_csv=args.remap_csv,
        remap_counts_csv=args.remap_counts_csv,
        remap_conseq_csv=args.remap_conseq_csv,
        unmapped1=args.unmapped1,
        unmapped2=args.unmapped2,
        work_path=work_path,
        gzip=args.gzip)  # defaults to False


if __name__ == '__main__':
    main()
