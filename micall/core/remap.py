#! /usr/bin/env python3.4

"""
Takes preliminary SAM as CSV input.  Iterative re-mapping of reads from
original FASTQ files.
Also report the number of reads mapped before and after processing.
"""

import argparse
from collections import Counter, defaultdict
import csv
from functools import partial
import logging
import multiprocessing
import os
import re
import shutil
import sys
from threading import Semaphore

# noinspection PyUnresolvedReferences
from gotoh import align_it
import Levenshtein

from micall.core import miseq_logging, project_config
from micall.core.sam2aln import apply_cigar, merge_pairs, merge_inserts
from micall.core.prelim_map import BOWTIE_THREADS, BOWTIE_BUILD_PATH, \
    BOWTIE_PATH, BOWTIE_VERSION, READ_GAP_OPEN, READ_GAP_EXTEND, REF_GAP_OPEN, \
    REF_GAP_EXTEND
from micall.utils.externals import Bowtie2, Bowtie2Build, LineCounter
from micall.utils.translation import reverse_and_complement

CONSENSUS_Q_CUTOFF = 20         # Min Q for base to contribute to conseq (pileup2conseq)
MIN_MAPPING_EFFICIENCY = 0.95   # Fraction of fastq reads mapped needed
MAX_REMAPS = 3                  # Number of remapping attempts if mapping efficiency unsatisfied
READ_MERGE_CHUNK_SIZE = 100
READ_MERGE_BUFFER_SIZE = BOWTIE_THREADS * READ_MERGE_CHUNK_SIZE * 10
SAM_FLAG_IS_UNMAPPED = 0x4
SAM_FLAG_IS_MATE_UNMAPPED = 0x8
SAM_FLAG_IS_FIRST_SEGMENT = 0x40

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

logger = miseq_logging.init_logging_console_only(logging.DEBUG)
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
                   worker_pool=None,
                   filter_coverage=1,
                   distance_report=None):
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
    @param worker_pool: a pool to do some distributed processing
    @param filter_coverage: when filtering on consensus distance, only include
        portions with at least this depth of coverage
    @param distance_report: empty dictionary or None. Dictionary will return:
        {rname: {'seed_dist': seed_dist, 'other_dist': other_dist,
        'other_seed': other_seed}}
    @return: {reference_name: consensus_sequence}
    """

    if debug_reports:
        for key in debug_reports.keys():
            debug_reports[key] = Counter()

    # refmap structure: {refname: {pos: Counter({nuc: count})}}
    refmap = {}

    semaphore = Semaphore(READ_MERGE_BUFFER_SIZE) if worker_pool is not None else None
    pairs = matchmaker(samfile, include_singles=True, semaphore=semaphore)
    if worker_pool is None:
        merged_reads = map(
            partial(merge_reads, quality_cutoff),
            pairs)
    else:
        merged_reads = worker_pool.imap_unordered(
            partial(merge_reads, quality_cutoff),
            pairs,
            chunksize=READ_MERGE_CHUNK_SIZE)
    read_counts = Counter()
    for merged_read in merged_reads:
        if semaphore is not None:
            semaphore.release()
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
    relevant_conseqs = None
    is_filtering = seeds and is_filtered

    gap_open_penalty = 15
    gap_extend_penalty = 3
    use_terminal_gap_penalty = 1
    while is_filtering and len(new_conseqs) > 1:
        drifted_seeds = []  # [(count, name)]
        if relevant_conseqs is None:
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
        for name in sorted(new_conseqs.keys()):
            relevant_conseq = relevant_conseqs[name]
            if not relevant_conseq:
                # None of the coverage was acceptable.
                drifted_seeds.append((read_counts[name], name))
                continue

            other_seed = other_dist = None
            for seed_name in sorted(new_conseqs.keys()):
                seed_ref = seeds[seed_name]
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
            del new_conseqs[dropped_seed]
        is_filtering = len(drifted_seeds) > 1
    return new_conseqs


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
                deletion += '-'
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
                  worker_pool=None,
                  filter_coverage=1,
                  distance_report=None):
    """ Build the new consensus sequences from the mapping results.

    @param samfilename: the mapping results in SAM format
    @param seeds: {name: sequence} If this is set,
        any positions without coverage will be set to the base from the seed
        reference. If there are no reads mapped to a reference, it will not
        be included as a new consensus.
    @param is_filtered: if True, then any consensus that has migrated so far
        from its seed that it is closer to a different seed, will not be
        included as a new consensus.
    @param worker_pool: a pool to do some distributed processing
    @param filter_coverage: when filtering on consensus distance, only include
        portions with at least this depth of coverage
    @param distance_report: empty dictionary or None. Dictionary will return:
        {rname: {'seed_dist': seed_dist, 'other_dist': other_dist,
        'other_seed': other_seed}}
    @return: {reference_name: consensus_sequence}
    """
    with open(samfilename, 'rU') as samfile:
        conseqs = sam_to_conseqs(samfile,
                                 CONSENSUS_Q_CUTOFF,
                                 seeds=seeds,
                                 is_filtered=is_filtered,
                                 worker_pool=worker_pool,
                                 filter_coverage=filter_coverage,
                                 distance_report=distance_report)

    return conseqs


def write_remap_counts(remap_counts_writer, counts, title, distance_report=None):
    distance_report = distance_report or {}
    for refname in sorted(counts.keys()):
        row = distance_report.get(refname, {})
        row['type'] = title + ' ' + refname
        row['count'] = counts[refname]
        remap_counts_writer.writerow(row)


def remap(fastq1,
          fastq2,
          prelim_csv,
          remap_csv,
          remap_counts_csv,
          remap_conseq_csv,
          unmapped1,
          unmapped2,
          work_path='',
          nthreads=BOWTIE_THREADS,
          callback=None,
          count_threshold=10,
          rdgopen=READ_GAP_OPEN,
          rfgopen=REF_GAP_OPEN,
          stderr=sys.stderr,
          gzip=False,
          debug_file_prefix=None):
    """
    Iterative re-map reads from raw paired FASTQ files to a reference sequence set that
    is being updated as the consensus of the reads that were mapped to the last set.
    @param fastq1: input R1 FASTQ
    @param fastq2: input R2 FASTQ
    @param prelim_csv: input CSV output from prelim_csv()
    @param remap_csv:  output CSV, contents of bowtie2 SAM output
    @param remap_counts_csv:  output CSV, counts of reads mapped to regions
    @param remap_conseq_csv:  output CSV, sample- and region-specific consensus sequences
                                generated while remapping reads
    @param unmapped1:  output FASTQ containing R1 reads that did not map to any region
    @param unmapped2:  output FASTQ containing R2 reads that did not map to any region
    @param work_path:  optional path to store working files
    @param nthreads:  optional setting to modify the number of threads used by bowtie2
    @param callback: a function to report progress with three optional
        parameters - callback(message, progress, max_progress)
    @param count_threshold:  minimum number of reads that map to a region for it to be remapped
    @param rdgopen: read gap open penalty
    @param rfgopen: reference gap open penalty
    @param stderr: an open file object to receive stderr from the bowtie2 calls
    @param gzip: True if the FASTQ files are gzipped
    @param debug_file_prefix: the prefix for the file path to write debug files.
        If not None, this will be used to write a copy of the reference FASTA
        files and the output SAM files.
    """

    reffile = os.path.join(work_path, 'temp.fasta')
    samfile = os.path.join(work_path, 'temp.sam')

    try:
        bowtie2 = Bowtie2(BOWTIE_VERSION, BOWTIE_PATH)
        bowtie2_build = Bowtie2Build(BOWTIE_VERSION,
                                     BOWTIE_BUILD_PATH,
                                     logger)
    except RuntimeError:
        bowtie2 = Bowtie2(BOWTIE_VERSION, BOWTIE_PATH + '-' + BOWTIE_VERSION)
        bowtie2_build = Bowtie2Build(BOWTIE_VERSION,
                                     BOWTIE_BUILD_PATH + '-' + BOWTIE_VERSION,
                                     logger)

    # check that the inputs exist
    if not os.path.exists(fastq1):
        logger.error('No FASTQ found at %s', fastq1)
        sys.exit(1)

    if not os.path.exists(fastq2):
        logger.error('No FASTQ found at %s', fastq2)
        sys.exit(1)

    # append .gz extension if necessary
    if gzip:
        if not fastq1.endswith('.gz'):
            try:
                os.symlink(fastq1, fastq1+'.gz')
            except OSError:
                # symbolic link already exists
                pass
            fastq1 += '.gz'

        if not fastq2.endswith('.gz'):
            try:
                os.symlink(fastq2, fastq2+'.gz')
            except OSError:
                # symbolic link already exists
                pass
            fastq2 += '.gz'

    worker_pool = multiprocessing.Pool(processes=nthreads) if nthreads > 1 else None

    # retrieve reference sequences used for preliminary mapping
    projects = project_config.ProjectConfig.loadDefault()
    seeds = projects.getAllReferences()

    # record the raw read count
    raw_count = line_counter.count(fastq1, gzip=gzip) // 2  # 4 lines per record in FASTQ, paired

    remap_counts_writer = csv.DictWriter(
        remap_counts_csv,
        'type count filtered_count seed_dist other_dist other_seed'.split(),
        lineterminator=os.linesep)
    remap_counts_writer.writeheader()
    remap_counts_writer.writerow(dict(type='raw', count=raw_count))

    # convert preliminary CSV to SAM, count reads
    with open(samfile, 'w') as f:
        refgroups = convert_prelim(prelim_csv,
                                   f,
                                   remap_counts_writer,
                                   count_threshold,
                                   projects)

    seed_counts = {best_ref: best_count
                   for best_ref, best_count in refgroups.values()}
    # regenerate consensus sequences based on preliminary map
    conseqs = build_conseqs(samfile, seeds=seeds, worker_pool=worker_pool)

    # exclude references with low counts (post filtering)
    new_conseqs = {}
    map_counts = {}
    for rname, conseq in conseqs.items():
        count = seed_counts.get(rname, None)
        if count is not None:
            map_counts[rname] = count  # transfer filtered counts to map counts for remap loop
            new_conseqs[rname] = conseq
    conseqs = new_conseqs

    # start remapping loop
    n_remaps = 0
    new_counts = Counter()
    unmapped_count = raw_count
    while conseqs:
        if callback:
            callback(message='... remap iteration %d' % n_remaps, progress=0)

        # reset unmapped files with each iteration
        unmapped1.seek(0)
        unmapped1.truncate()
        unmapped2.seek(0)
        unmapped2.truncate()

        if debug_file_prefix is None:
            next_debug_prefix = None
        else:
            next_debug_prefix = '{}_remap{}'.format(debug_file_prefix,
                                                    n_remaps+1)
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
                                          nthreads,
                                          new_counts,
                                          stderr,
                                          callback,
                                          debug_file_prefix=next_debug_prefix)

        old_seed_names = set(conseqs.keys())
        # regenerate consensus sequences
        distance_report = {}
        conseqs = build_conseqs(samfile,
                                seeds=seeds,
                                is_filtered=True,
                                worker_pool=worker_pool,
                                filter_coverage=count_threshold//2,  # pairs
                                distance_report=distance_report)
        new_seed_names = set(conseqs.keys())
        n_remaps += 1
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

    # finished iterative phase
    if worker_pool is not None:
        worker_pool.close()

    # generate SAM CSV output
    remap_writer = csv.DictWriter(remap_csv, SAM_FIELDS, lineterminator=os.linesep)
    remap_writer.writeheader()
    if new_counts:
        splitter = MixedReferenceSplitter(work_path)
        split_counts = Counter()
        # At least one read was mapped, so samfile has relevant data
        with open(samfile, 'rU') as f:
            for fields in splitter.split(f):
                remap_writer.writerow(dict(zip(SAM_FIELDS, fields)))
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
                                               nthreads,
                                               split_counts,
                                               stderr,
                                               callback)
            new_counts.update(split_counts)
            with open(samfile, 'rU') as f:
                for fields in splitter.walk(f):
                    remap_writer.writerow(dict(zip(SAM_FIELDS, fields)))

    # write consensus sequences and counts
    remap_conseq_csv.write('region,sequence\n')  # record consensus sequences for later use
    for refname in new_counts.keys():
        # NOTE this is the consensus sequence to which the reads were mapped, NOT the
        # current consensus!
        conseq = conseqs.get(refname) or projects.getReference(refname)
        remap_conseq_csv.write('%s,%s\n' % (refname, conseq))
    write_remap_counts(remap_counts_writer,
                       new_counts,
                       title='remap-final')

    # report number of unmapped reads
    remap_counts_writer.writerow(dict(type='unmapped',
                                      count=unmapped_count))


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
    :return dict: { group_name: (refname, count) }
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
        counts = ref_counts[row['rname']]
        counts[1] += 1  # full count

        # write SAM row
        target.write('\t'.join([row[field] for field in SAM_FIELDS]) + '\n')

        if is_unmapped_read(row['flag']):
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
    return refgroups


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
                     nthreads,
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
    @param nthreads:  optional setting to modify the number of threads used by bowtie2
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
                   '-p', str(nthreads)]

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
                    yield fields[:11]
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


def matchmaker(samfile, include_singles=False, semaphore=None):
    """
    An iterator that returns pairs of reads sharing a common qname from a SAM file.
    Note that unpaired reads will be left in the cached_rows dictionary and
    discarded.
    @param samfile: open file handle to a SAM file
    @param include_singles: True if unpaired reads should be returned, paired
        with a None value: ([qname, flag, rname, ...], None)
    @param semaphore: acquired before each pair is generated
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
                if semaphore is not None:
                    semaphore.acquire()
                # current row should be the second read of the pair
                yield old_row, row
    if include_singles:
        for row in cached_rows.values():
            if semaphore is not None:
                semaphore.acquire()
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
    parser.add_argument('prelim_csv',
                        type=argparse.FileType('rU'),
                        help='<input> CSV containing preliminary map output (modified SAM)')
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
    remap(fastq1=args.fastq1,
          fastq2=args.fastq2,
          prelim_csv=args.prelim_csv,
          remap_csv=args.remap_csv,
          remap_counts_csv=args.remap_counts_csv,
          remap_conseq_csv=args.remap_conseq_csv,
          unmapped1=args.unmapped1,
          unmapped2=args.unmapped2,
          work_path=work_path,
          gzip=args.gzip)  # defaults to False

if __name__ == '__main__':
    main()
