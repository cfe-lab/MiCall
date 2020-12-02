#! /usr/bin/env python3.6

"""
Kive-style MiSeq pipeline, post-processing 1
Takes aligned CSV as input (produced by sam2aln).
Re-aligns sequence variants to lab standard references (e.g., HXB2).
Reports nucleotide and amino acid frequencies by reference coordinates.
Outputs consensus sequences in same coordinate system.
This assumes a consistent reading frame across the entire region.
"""

import argparse
import re
import typing
from collections import Counter, defaultdict, OrderedDict
import csv
from csv import DictWriter
from enum import IntEnum
from itertools import groupby
from operator import itemgetter, attrgetter
import os
from pathlib import Path

import gotoh
import yaml
from mappy import Aligner

from micall.core.project_config import ProjectConfig, G2P_SEED_NAME
from micall.core.remap import PARTIAL_CONTIG_SUFFIX, REVERSED_CONTIG_SUFFIX
from micall.utils.alignment_wrapper import align_nucs
from micall.utils.big_counter import BigCounter
from micall.utils.spring_beads import Wire, Bead
from micall.utils.translation import translate, ambig_dict

AMINO_ALPHABET = 'ACDEFGHIKLMNPQRSTVWY*'
CONSEQ_MIXTURE_CUTOFFS = [0.01, 0.02, 0.05, 0.1, 0.2, 0.25]
GAP_OPEN_COORD = 40
GAP_EXTEND_COORD = 10
CONSENSUS_MIN_COVERAGE = 100
MAX_CUTOFF = 'MAX'
FIRST_CUTOFF = 'FIRST'


CigarActions = IntEnum(
    'CigarActions',
    'MATCH INSERT DELETE SKIPPED SOFT_CLIPPED HARD_CLIPPED',
    start=0)


def parse_args():
    parser = argparse.ArgumentParser(
        description='Post-processing of short-read alignments.')

    parser.add_argument('aligned_csv',
                        type=argparse.FileType(),
                        help='input CSV with aligned reads')
    parser.add_argument('remap_conseq_csv',
                        type=argparse.FileType(),
                        help='input CSV with consensus sequences from remap step')
    parser.add_argument('--clipping_csv',
                        type=argparse.FileType(),
                        help='input CSV with count of soft-clipped reads at each position')
    parser.add_argument('--conseq_ins_csv',
                        type=argparse.FileType(),
                        help='input CSV with insertions relative to consensus sequence')
    parser.add_argument('--contigs_csv',
                        type=argparse.FileType(),
                        help='input CSV with assembled contigs')
    parser.add_argument('--g2p_aligned_csv',
                        type=argparse.FileType(),
                        help='CSV of aligned reads from the G2P process')
    parser.add_argument('--nuc_csv',
                        type=argparse.FileType('w'),
                        default=os.devnull,
                        help='CSV containing nucleotide frequencies')
    parser.add_argument('--nuc_detail_csv',
                        type=argparse.FileType('w'),
                        help='CSV containing nucleotide frequencies for each '
                             'contig')
    parser.add_argument('--amino_csv',
                        type=argparse.FileType('w'),
                        default=os.devnull,
                        help='CSV containing amino frequencies')
    parser.add_argument('--amino_detail_csv',
                        type=argparse.FileType('w'),
                        help='CSV containing amino frequencies for each contig')
    parser.add_argument('--coord_ins_csv',
                        type=argparse.FileType('w'),
                        default=os.devnull,
                        help='CSV containing insertions relative to coordinate reference')
    parser.add_argument('--conseq_csv',
                        type=argparse.FileType('w'),
                        default=os.devnull,
                        help='CSV containing consensus sequences')
    parser.add_argument('--conseq_all_csv',
                        type=argparse.FileType('w'),
                        help='CSV containing consensus sequences (ignoring inadequate '
                             'coverage)')
    parser.add_argument('--failed_align_csv',
                        type=argparse.FileType('w'),
                        default=os.devnull,
                        help='CSV containing any consensus that failed to align')
    parser.add_argument('--conseq_region_csv',
                        type=argparse.FileType('w'),
                        default=os.devnull,
                        help='CSV containing consensus sequences, split by region')
    parser.add_argument('--genome_coverage_csv',
                        type=argparse.FileType('w'),
                        help='CSV of coverage levels in full-genome coordinates')
    parser.add_argument('--coverage_summary_csv',
                        type=argparse.FileType('w'),
                        help='CSV of coverage levels for the whole sample')
    parser.add_argument('--minimap_hits_csv',
                        type=argparse.FileType('w'),
                        help='CSV of minimap2 match locations')
    return parser.parse_args()


def trim_contig_name(contig_name):
    prefix = contig_name.split('-')[0]
    if not re.fullmatch(r'[\d_]+', prefix):
        seed_name = contig_name
    else:
        _, seed_name = contig_name.split('-', 1)
    return seed_name


class SequenceReport(object):
    """ Hold the data for several reports related to a sample's genetic sequence.

    To use a report object, read a group of aligned reads that mapped to a
    single region, and then write out all the reports for that region.
    """
    def __init__(self,
                 insert_writer,
                 projects,
                 conseq_mixture_cutoffs,
                 clipping_counts=None,
                 conseq_insertion_counts=None,
                 landmarks_yaml=None):
        """ Create an object instance.

        :param insert_writer: InsertionWriter object that will track reads and
            write out any insertions relative to the coordinate reference.
        :param projects: ProjectConfig object
        :param conseq_mixture_cutoffs: a list of cutoff fractions used to
            determine what portion a variant must exceed before it will be
            included as a mixture in the consensus.
        """
        self.repeated_pos = None
        self.consensus_min_coverage = 0
        self.callback_progress = 0
        self.callback_next = self.callback_chunk_size = self.callback_max = None
        self.insert_writer = insert_writer
        self.projects = projects
        self.conseq_mixture_cutoffs = list(conseq_mixture_cutoffs)
        self.conseq_mixture_cutoffs.insert(0, MAX_CUTOFF)
        self.callback = None
        self.seed_aminos = self.reports = self.reading_frames = None
        self.inserts = self.consensus = self.seed = None
        self.contigs = None  # [(ref, group_ref, seq)]
        self.consensus_by_reading_frame = None  # {frame: seq}
        if landmarks_yaml is None:
            landmarks_path = (Path(__file__).parent.parent / 'data' /
                              'landmark_references.yaml')
            landmarks_yaml = landmarks_path.read_text()
        self.landmarks = yaml.safe_load(landmarks_yaml)
        self.coordinate_refs = self.remap_conseqs = None

        # {seed: {coord_region: [ReportAmino]}}
        self.combined_reports = defaultdict(lambda: defaultdict(list))

        # {seed_name: {pos: count}}
        self.clipping_counts = clipping_counts or defaultdict(Counter)

        # {seed_name: {pos: count}
        self.conseq_insertion_counts = (conseq_insertion_counts or
                                        defaultdict(Counter))
        self.nuc_writer = self.nuc_detail_writer = self.conseq_writer = None
        self.amino_writer = self.amino_detail_writer = None
        self.genome_coverage_writer = self.minimap_hits_writer = None
        self.conseq_region_writer = self.fail_writer = None
        self.conseq_all_writer = None

    @property
    def has_detail_counts(self):
        return (self.amino_detail_writer is not None or
                self.nuc_detail_writer is not None)

    def enable_callback(self, callback, file_size):
        """ Enable callbacks to update progress while counting reads.

        @param callback: a function to report progress with three optional
            parameters - callback(message, progress, max_progress)
        @param file_size: the size of the aligned reads file
        """

        self.callback = callback
        self.callback_max = file_size
        self.callback_chunk_size = file_size // 100
        self.callback_next = self.callback_chunk_size
        self.callback(message='... extracting statistics from alignments',
                      progress=0,
                      max_progress=self.callback_max)

    def _count_reads(self, aligned_reads):
        """
        Parses contents of aligned CSV.

        :param aligned_reads: a sequence of Dicts from csv.DictReader
        """

        for i, row in enumerate(aligned_reads):
            if i == 0:
                # these will be the same for all rows, so just assign from the first
                self.detail_seed = row['refname']
                self.seed = trim_contig_name(self.detail_seed)
                self.qcut = row['qcut']
            nuc_seq = row['seq']
            offset = int(row['offset'])
            count = int(row['count'])
            if self.callback:
                row_size = sum(map(len, row.values()))
                self.callback_progress += row_size
                if self.callback_progress >= self.callback_next:
                    self.callback(progress=self.callback_progress)
                    self.callback_next += self.callback_chunk_size

            # first run, prepare containers
            if not self.seed_aminos:
                self.insert_writer.start_group(self.seed,
                                               self.qcut)
                for reading_frame in range(3):
                    self.seed_aminos[reading_frame] = []

            # record this read to calculate insertions later
            self.insert_writer.add_nuc_read('-'*offset + nuc_seq, count)

            # cycle through reading frames
            for reading_frame, frame_seed_aminos in self.seed_aminos.items():
                offset_nuc_seq = ' ' * (reading_frame + offset) + nuc_seq
                # pad to a codon boundary
                offset_nuc_seq += ' ' * ((3 - (len(offset_nuc_seq) % 3)) % 3)
                start = offset - (offset % 3)
                for nuc_pos in range(start, len(offset_nuc_seq), 3):
                    codon_index = nuc_pos // 3
                    # append SeedAmino objects to this list if necessary
                    while len(frame_seed_aminos) <= codon_index:
                        frame_seed_aminos.append(SeedAmino(
                            len(frame_seed_aminos)*3-reading_frame))

                    # update amino acid counts
                    codon = offset_nuc_seq[nuc_pos:nuc_pos + 3]
                    seed_amino = frame_seed_aminos[codon_index]
                    seed_amino.count_aminos(codon, count)
        if self.callback:
            self.callback(progress=self.callback_max)

    def _pair_align(self, reference, query, gap_open=15, gap_extend=5, use_terminal_gap_penalty=1):
        """ Align a query sequence of amino acids to a reference sequence.

        @return: (aligned_ref, aligned_query, score)
        """
        # noinspection PyUnresolvedReferences
        aligned_ref, aligned_query, score = gotoh.align_it_aa(
            reference,
            query,
            gap_open,
            gap_extend,
            use_terminal_gap_penalty)
        return aligned_ref, aligned_query, score

    def _map_to_coordinate_ref(self, coordinate_name, coordinate_ref):
        """
        Align consensus of aligned reads in each reading frame to the
        coordinate reference to determine best frame.
        Generate a map of consensus indices to the coordinate reference.

        :param coordinate_name: Name of coordinate reference.
        :param coordinate_ref: Amino acid sequence of coordinate reference.
        :return: Dict object containing coordinate map.
        """

        # Start max_score with the minimum score we can consider a valid
        # alignment. Anything worse, we won't bother
        consensus_length = len([amino for amino in self.seed_aminos[0] if amino.counts])
        max_score = min(consensus_length, len(coordinate_ref))

        best_alignment = None
        for reading_frame in self.seed_aminos:
            frame_seed_aminos = self.get_seed_aminos(reading_frame)
            consensus = ''.join([seed_amino1.get_consensus()
                                 for seed_amino1 in frame_seed_aminos])
            if reading_frame == 0:
                # best guess before aligning - if alignments fail, this will be non-empty
                self.consensus[coordinate_name] = consensus

            # map to reference coordinates by aligning consensus
            a_coord, a_consensus, score = self._pair_align(
                coordinate_ref,
                consensus,
                gap_open=GAP_OPEN_COORD,
                gap_extend=GAP_EXTEND_COORD)
            if score < max_score:
                continue
            max_score = score
            best_alignment = (reading_frame,
                              consensus,
                              a_coord,
                              a_consensus,
                              frame_seed_aminos)

        report_aminos = []
        if best_alignment is not None:
            (reading_frame,
             consensus,
             a_coord,
             a_consensus,
             frame_seed_aminos) = best_alignment

            ref_offset = 0
            prev_nuc_index = None
            for seed_amino in frame_seed_aminos:
                if seed_amino.consensus_nuc_index is None:
                    pass
                elif prev_nuc_index is None:
                    pass
                elif seed_amino.consensus_nuc_index - prev_nuc_index == 2:
                    ref_offset = -1
                    seed_amino.nucleotides_to_skip = 1
                seed_amino.ref_offset = ref_offset
                prev_nuc_index = seed_amino.consensus_nuc_index

            coord2conseq = self.map_sequences(coordinate_ref,
                                              consensus,
                                              a_coord,
                                              a_consensus)
            self.reading_frames[coordinate_name] = reading_frame
            self.consensus[coordinate_name] = consensus
            coordinate_inserts = {seed_amino.consensus_nuc_index
                                  for seed_amino in frame_seed_aminos}
            prev_conseq_index = None
            prev_consensus_nuc_index = None
            for coord_index in range(len(coordinate_ref)):
                conseq_index = coord2conseq.get(coord_index)
                if conseq_index is None:
                    seed_amino = SeedAmino(None)
                    if (prev_conseq_index is not None and
                            prev_conseq_index+1 < len(frame_seed_aminos)):
                        prev_seed_amino = frame_seed_aminos[prev_conseq_index]
                        prev_count = sum(prev_seed_amino.counts.values())
                        prev_count += prev_seed_amino.deletions
                        next_seed_amino = frame_seed_aminos[prev_conseq_index+1]
                        next_count = sum(next_seed_amino.counts.values())
                        next_count += next_seed_amino.deletions
                        min_count = min(prev_count, next_count)
                        seed_amino.deletions = min_count
                        for nuc in seed_amino.nucleotides:
                            nuc.count_nucleotides('-', min_count)
                else:
                    seed_amino = frame_seed_aminos[conseq_index]
                    if prev_conseq_index is None:
                        coordinate_inserts = {i
                                              for i in coordinate_inserts
                                              if i >= seed_amino.consensus_nuc_index}
                    prev_conseq_index = conseq_index

                report_aminos.append(ReportAmino(seed_amino, coord_index + 1))
                if seed_amino.consensus_nuc_index is not None:
                    coordinate_inserts.remove(seed_amino.consensus_nuc_index)
                    prev_consensus_nuc_index = seed_amino.consensus_nuc_index
            if prev_consensus_nuc_index is None:
                coordinate_inserts.clear()
            else:
                coordinate_inserts = {i
                                      for i in coordinate_inserts
                                      if i <= prev_consensus_nuc_index}
            self.inserts[coordinate_name] = coordinate_inserts

        self.reports[coordinate_name] = report_aminos

    def get_seed_aminos(self, reading_frame: int) -> typing.List['SeedAmino']:
        if self.repeated_pos is None:
            return self.seed_aminos[reading_frame]

        seed_aminos = []
        has_passed_repeat = False
        pos = 1
        while True:
            try:
                seed_amino: SeedAmino = self.seed_aminos[reading_frame][pos-1]
            except IndexError:
                break
            past_repeat = seed_amino.consensus_nuc_index - self.repeated_pos
            if 0 <= past_repeat <= 2 and not has_passed_repeat:
                has_passed_repeat = True
                new_reading_frame = (reading_frame + 1) % 3
                if new_reading_frame < reading_frame:
                    pos -= 1
                reading_frame = new_reading_frame
                try:
                    seed_amino = self.seed_aminos[reading_frame][pos-1]
                except IndexError:
                    break
            seed_aminos.append(seed_amino)
            pos += 1
        return seed_aminos

    def map_sequences(self, from_seq, to_seq, from_aligned=None, to_aligned=None):
        if from_aligned is None or to_aligned is None:
            from_aligned, to_aligned, _score = self._pair_align(
                from_seq,
                to_seq,
                gap_open=GAP_OPEN_COORD,
                gap_extend=GAP_EXTEND_COORD)
        seq_map = {}
        from_index = to_index = 0
        for from_aa, to_aa in zip(from_aligned, to_aligned):
            if (to_index < len(to_seq) and
                    to_aa == to_seq[to_index]):
                seq_map[from_index] = to_index
                to_index += 1
            if (from_index < len(from_seq) and
                    from_aa == from_seq[from_index]):
                from_index += 1
        return seq_map

    def process_reads(self,
                      aligned_csv,
                      coverage_summary=None,
                      included_regions: typing.Optional[typing.Set] = None,
                      excluded_regions: typing.Set = frozenset()):
        """ Parse CSV file containing aligned reads.

        Grouped by reference and quality cutoff.

        @param aligned_csv: an open aligned.CSV file where each line
            corresponds to a single aligned read.
        @param coverage_summary: an open CSV file where the coverage summary
            will be written, or None to ignore.
        @param included_regions: coordinate regions that should be reported,
            all other regions should be excluded, or None to ignore
        @param excluded_regions: coordinate regions that should not be reported.
        """
        aligned_reader = csv.DictReader(aligned_csv)
        for _, aligned_reads in groupby(aligned_reader,
                                        lambda row: (row['refname'], row['qcut'])):
            self.read(aligned_reads,
                      included_regions=included_regions,
                      excluded_regions=excluded_regions)

            if self.insert_writer is not None:
                self.write_insertions(self.insert_writer)
            if self.nuc_detail_writer is not None:
                self.write_nuc_detail_counts(self.nuc_detail_writer)
            elif self.nuc_writer is not None:
                self.write_nuc_counts(self.nuc_writer)
            if self.conseq_writer is not None:
                self.write_consensus(self.conseq_writer)
            if self.conseq_all_writer is not None:
                self.write_consensus_all(self.conseq_all_writer)
            if self.conseq_region_writer is not None:
                self.write_consensus_regions(self.conseq_region_writer)
            if self.genome_coverage_writer is not None:
                self.write_genome_coverage_counts()
            if self.amino_detail_writer is not None:
                self.write_amino_detail_counts()
            elif self.amino_writer is not None:
                self.write_amino_counts(self.amino_writer,
                                        coverage_summary=coverage_summary)
            if self.fail_writer is not None:
                self.write_failure(self.fail_writer)
            if self.has_detail_counts:
                self.combine_reports()
        if self.nuc_detail_writer is not None and self.nuc_writer is not None:
            self.write_nuc_counts(self.nuc_writer)
        if self.amino_detail_writer is not None and self.amino_writer is not None:
            self.write_amino_counts(self.amino_writer,
                                    coverage_summary=coverage_summary)

    def read(self,
             aligned_reads,
             included_regions: typing.Optional[typing.Set] = None,
             excluded_regions: typing.Set = frozenset()):
        """
        Reset all the counters, and read a new section of aligned reads.
        A section must have the same region, and qcut on all lines.

        @param aligned_reads: an iterator of dicts generated by csv.DictReader
            Each dict corresponds to a row from an aligned.CSV file and
            corresponds to a single aligned read.
        @param included_regions: coordinate regions that should be reported,
            all other regions should be excluded, or None to ignore
        @param excluded_regions: coordinate regions that should not be reported.
        """
        aligned_reads = self.align_deletions(aligned_reads)

        self.seed_aminos = {}  # {reading_frame: [SeedAmino(consensus_nuc_index)]}
        self.reports = {}  # {coord_name: [ReportAmino()]}
        self.reading_frames = {}  # {coord_name: reading_frame}
        self.inserts = {}  # {coord_name: set([consensus_index])}
        self.consensus = {}  # {coord_name: consensus_amino_seq}

        # populates these dictionaries, generates amino acid counts
        self._count_reads(aligned_reads)

        if (self.seed is None or
                self.seed.endswith(PARTIAL_CONTIG_SUFFIX) or
                self.seed.endswith(REVERSED_CONTIG_SUFFIX)):
            return

        if not self.seed_aminos:
            # no valid reads were aligned to this region and counted, skip next step
            self.coordinate_refs = {}
        else:
            self.coordinate_refs = self.projects.getCoordinateReferences(self.seed)
            if not self.coordinate_refs:
                # No coordinate references defined for this region.
                # Pad the amino acid count list until it has the same length
                # as the seed reference, then skip next step.
                seed_ref = self.projects.getReference(self.seed)
                while len(self.seed_aminos[0])*3 < len(seed_ref):
                    self.seed_aminos[0].append(SeedAmino(3*len(self.seed_aminos[0])))

        self.find_repeated_pos()

        # iterate over coordinate references defined for this region
        for coordinate_name, coordinate_ref in self.coordinate_refs.items():
            if included_regions is not None and \
                    coordinate_name not in included_regions:
                continue
            if coordinate_name in excluded_regions:
                continue
            self._map_to_coordinate_ref(coordinate_name, coordinate_ref)

    def find_repeated_pos(self):
        self.repeated_pos = None
        repeated_ref_name = 'SARS-CoV-2-seed'
        repeated_ref_pos = 13468  # 1-based position of duplicated base
        if self.seed != repeated_ref_name:
            return

        seed_start = repeated_ref_pos-100
        seed_end = seed_start + 200
        seed_ref = self.projects.getReference(self.seed)[seed_start:seed_end]
        consensus = ''.join(nuc.get_consensus('MAX', no_coverage='x')
                            for amino in self.seed_aminos[0]
                            for nuc in amino.nucleotides)
        consensus_stripped = consensus.replace('x', '')
        aref, aseq, score = align_nucs(seed_ref, consensus_stripped)
        offset = len(aseq) - len(aseq.lstrip('-'))
        aref = aref[offset:]
        aseq = aseq[offset:]
        ref_pos = offset+seed_start
        conseq_pos = 0
        for ref_nuc, conseq_nuc in zip(aref, aseq):
            while conseq_pos < len(consensus) and consensus[conseq_pos] == 'x':
                conseq_pos += 1
            if conseq_pos >= len(consensus):
                break

            if ref_nuc == seed_ref[ref_pos-seed_start]:
                ref_pos += 1
            if conseq_nuc == consensus[conseq_pos]:
                conseq_pos += 1
            if ref_pos == repeated_ref_pos:
                self.repeated_pos = conseq_pos
                break
            if ref_pos >= len(seed_ref) + seed_start:
                break

    def read_clipping(self, clipping_csv):
        for row in csv.DictReader(clipping_csv):
            pos = int(row['pos'])
            count = int(row['count'])
            self.clipping_counts[row['refname']][pos] += count

    def read_insertions(self, conseq_ins_csv):
        reader = csv.DictReader(conseq_ins_csv)

        # {ref: {pos: set([qname])}}
        insertion_names = defaultdict(lambda: defaultdict(set))

        for row in reader:
            ref_name = row['refname']
            pos = int(row['pos'])
            pos_names = insertion_names[ref_name][pos]
            pos_names.add(row['qname'])
        for ref_name, ref_positions in insertion_names.items():
            self.conseq_insertion_counts[ref_name] = Counter(
                {pos: len(names) for pos, names in ref_positions.items()})

    @staticmethod
    def _create_amino_writer(amino_file):
        columns = ['seed',
                   'region',
                   'q-cutoff',
                   'query.nuc.pos',
                   'refseq.aa.pos']
        columns.extend(AMINO_ALPHABET)
        columns.extend(
            ('X', 'partial', 'del', 'ins', 'clip', 'v3_overlap', 'coverage'))
        return csv.DictWriter(amino_file,
                              columns,
                              lineterminator=os.linesep)

    def write_amino_header(self, amino_file):
        self.amino_writer = self._create_amino_writer(amino_file)
        self.amino_writer.writeheader()

    def write_amino_counts(self, amino_writer=None, coverage_summary=None):
        """ Write amino counts file.

        Must have already called write_nuc_counts() to calculate max_clip_count
        for each ReportAmino.
        """
        amino_writer = amino_writer or self.amino_writer
        if not self.combined_reports:
            self.write_amino_report(amino_writer,
                                    self.reports,
                                    self.seed,
                                    coverage_summary)
        else:
            for seed, reports in self.combined_reports.items():
                self.write_amino_report(amino_writer,
                                        reports,
                                        seed,
                                        coverage_summary)

    def write_amino_detail_header(self, amino_detail_file):
        self.amino_detail_writer = self._create_amino_writer(amino_detail_file)
        self.amino_detail_writer.writeheader()

    def write_amino_detail_counts(self, amino_detail_writer=None):
        amino_detail_writer = amino_detail_writer or self.amino_detail_writer
        self.write_amino_report(amino_detail_writer, self.reports, self.detail_seed)

    def combine_reports(self):
        if self.contigs is None:
            group_ref = self.seed
        elif self.detail_seed == G2P_SEED_NAME:
            group_ref = G2P_SEED_NAME
        else:
            contig_num, _ = self.detail_seed.split('-', 1)
            group_ref = self.contigs[int(contig_num)-1][1]
        old_reports = self.combined_reports[group_ref]
        for region, report in self.reports.items():
            old_report = old_reports[region]
            for i, report_amino in enumerate(report):
                while len(old_report) <= i:
                    old_report.append(ReportAmino(SeedAmino(None),
                                                  len(old_report)+1))
                old_report[i].seed_amino.add(report_amino.seed_amino)
        self.reports.clear()

    def write_amino_report(self, amino_writer, reports, seed, coverage_summary=None):
        if not reports:
            return
        regions = sorted(reports.keys())
        for region in regions:
            coverage_sum = 0.0
            pos_count = 0
            for report_amino in reports[region]:
                seed_amino = report_amino.seed_amino
                query_pos = (str(seed_amino.consensus_nuc_index + 1)
                             if seed_amino.consensus_nuc_index is not None
                             else '')
                max_clip_count = 0
                total_insertion_count = 0
                for seed_nuc in seed_amino.nucleotides:
                    max_clip_count = max(seed_nuc.clip_count, max_clip_count)
                    total_insertion_count += seed_nuc.insertion_count
                row = {'seed': seed,
                       'region': region,
                       'q-cutoff': self.qcut,
                       'query.nuc.pos': query_pos,
                       'refseq.aa.pos': report_amino.position,
                       'X': seed_amino.low_quality,
                       'partial': seed_amino.partial,
                       'del': seed_amino.deletions,
                       'ins': total_insertion_count,
                       'clip': max_clip_count,
                       'v3_overlap': seed_amino.v3_overlap,
                       'coverage': seed_amino.deletions}
                for letter in AMINO_ALPHABET:
                    letter_count = seed_amino.counts[letter]
                    row[letter] = letter_count
                    coverage_sum += letter_count
                    row['coverage'] += letter_count
                amino_writer.writerow(row)
                pos_count += 1
            if coverage_summary is not None and pos_count > 0:
                region_coverage = coverage_sum / pos_count
                old_coverage = coverage_summary.get('avg_coverage', -1)
                if region_coverage > old_coverage:
                    coverage_summary['avg_coverage'] = region_coverage
                    coverage_summary['coverage_region'] = region
                    coverage_summary['region_width'] = pos_count

    def write_minimap_hits_header(self, minimap_hits_file):
        self.minimap_hits_writer = csv.DictWriter(minimap_hits_file,
                                                  ['contig',
                                                   'ref_name',
                                                   'start',
                                                   'end',
                                                   'ref_start',
                                                   'ref_end'],
                                                  lineterminator=os.linesep)
        self.minimap_hits_writer.writeheader()

    def write_genome_coverage_header(self, genome_coverage_file):
        columns = ['contig',
                   'coordinates',
                   'query_nuc_pos',
                   'refseq_nuc_pos',
                   'dels',
                   'coverage',
                   'link']
        self.genome_coverage_writer = csv.DictWriter(genome_coverage_file,
                                                     columns,
                                                     lineterminator=os.linesep)
        self.genome_coverage_writer.writeheader()

    def write_genome_coverage_counts(self):
        seed_nucs = [(seed_nuc.get_consensus(MAX_CUTOFF) or '?', seed_nuc)
                     for seed_amino in self.seed_aminos[0]
                     for seed_nuc in seed_amino.nucleotides]
        aligned_consensus = ''.join(nuc for nuc, coverage in seed_nucs)
        consensus = aligned_consensus.lstrip('?')
        consensus_offset = len(aligned_consensus) - len(consensus)
        consensus = consensus.rstrip('?')

        # Fill in gaps from remap conseq.
        seed_seq = self.remap_conseqs and self.remap_conseqs.get(self.detail_seed)
        if seed_seq:
            consensus = ''.join(nuc if nuc != '?' else seed_seq[i+consensus_offset]
                                for i, nuc in enumerate(consensus))
        is_partial = (self.seed.endswith(PARTIAL_CONTIG_SUFFIX) or
                      self.seed.endswith(REVERSED_CONTIG_SUFFIX))
        if is_partial:
            coordinate_name = None
        else:
            for seed_landmarks in self.landmarks:
                if re.fullmatch(seed_landmarks['seed_pattern'], self.seed):
                    coordinate_name = seed_landmarks['coordinates']
                    break
            else:
                coordinate_name = None

        self.write_sequence_coverage_counts(self.detail_seed,
                                            coordinate_name,
                                            consensus,
                                            consensus_offset,
                                            seed_nucs)
        if is_partial or not self.contigs:
            return
        seed_prefix = self.detail_seed.split('-')[0]
        try:
            contig_nums = [int(item) for item in seed_prefix.split('_')]
        except ValueError:
            return
        for contig_num in contig_nums:
            contig_ref, group_ref, contig_seq = self.contigs[contig_num-1]
            contig_name = f'contig-{contig_num}-{contig_ref}'
            self.write_sequence_coverage_counts(contig_name,
                                                coordinate_name,
                                                contig_seq)

    def write_sequence_coverage_counts(self,
                                       contig_name: str,
                                       coordinate_name: typing.Optional[str],
                                       consensus: str,
                                       consensus_offset=0,
                                       seed_nucs=None):
        if coordinate_name is None:
            alignments = []
        else:
            coordinate_seq = self.projects.getReference(coordinate_name)
            aligner = Aligner(seq=coordinate_seq, preset='map-ont')
            alignments = sorted(aligner.map(consensus), key=attrgetter('q_st'))
            alignments = [alignment
                          for alignment in alignments
                          if alignment.is_primary]
        wire = Wire()
        bead_end = source_end = 0
        for alignment in alignments:
            curr_start = alignment.q_st
            if source_end < curr_start:
                bead_size = curr_start - source_end
                wire.append(Bead(start=bead_end, end=bead_end+bead_size))
                bead_end += bead_size
                skipped_source = 0
            else:
                skipped_source = source_end - curr_start
            bead_size = alignment.r_en - alignment.r_st - skipped_source
            wire.append(Bead(bead_end,
                             bead_end+bead_size,
                             alignment.r_st+skipped_source,
                             alignment.r_en,
                             alignment,
                             skipped_source))
            if self.minimap_hits_writer is not None:
                ref_start = alignment.r_st+1
                ref_end = alignment.r_en
                if alignment.strand < 0:
                    ref_start, ref_end = ref_end, ref_start
                self.minimap_hits_writer.writerow(dict(contig=contig_name,
                                                       ref_name=coordinate_name,
                                                       start=alignment.q_st+1,
                                                       end=alignment.q_en,
                                                       ref_start=ref_start,
                                                       ref_end=ref_end))
            bead_end += bead_size
            source_end = alignment.q_en
        expected_end = len(consensus)
        if source_end < expected_end:
            wire.append(Bead(start=source_end, end=expected_end))
        wire.align()
        seq_pos = skipped_source = 0
        for bead in wire:
            ref_pos = bead.display_start - bead.skipped
            skipped_source += bead.skipped
            seq_pos -= bead.skipped
            if bead.alignment is not None:
                cigar = bead.alignment.cigar
            else:
                cigar = [(bead.end - bead.start, None)]
            for length, action in cigar:
                if action == CigarActions.DELETE:
                    ref_pos += length
                    continue
                for _ in range(length):
                    if action == CigarActions.INSERT:
                        seq_pos += 1
                        ref_pos_display = None
                        link = 'I'
                    elif action == CigarActions.MATCH or action is None:
                        seq_pos += 1
                        ref_pos += 1
                        ref_pos_display = ref_pos
                        link = 'U' if action is None else 'M'
                    else:
                        name = CigarActions(action).name
                        raise ValueError(f'Unexpected CIGAR action: {name}.')
                    offset_seq_pos = seq_pos + consensus_offset
                    if not seed_nucs:
                        coverage = dels = None
                    else:
                        seed_nuc: SeedNucleotide = seed_nucs[offset_seq_pos - 1][1]
                        coverage = seed_nuc.get_coverage()
                        if coverage == 0:
                            continue
                        dels = seed_nuc.counts['-']
                    if 0 < skipped_source:
                        skipped_source -= 1
                    else:
                        row = dict(contig=contig_name,
                                   coordinates=coordinate_name,
                                   query_nuc_pos=offset_seq_pos,
                                   refseq_nuc_pos=ref_pos_display,
                                   dels=dels,
                                   coverage=coverage,
                                   link=link)
                        self.genome_coverage_writer.writerow(row)

    @staticmethod
    def _create_nuc_writer(nuc_file):
        return csv.DictWriter(nuc_file,
                              ['seed',
                               'region',
                               'q-cutoff',
                               'query.nuc.pos',
                               'refseq.nuc.pos',
                               'A',
                               'C',
                               'G',
                               'T',
                               'N',
                               'del',
                               'ins',
                               'clip',
                               'v3_overlap',
                               'coverage'],
                              lineterminator=os.linesep)

    def write_nuc_header(self, nuc_file):
        self.nuc_writer = self._create_nuc_writer(nuc_file)
        self.nuc_writer.writeheader()

    def write_nuc_detail_header(self, nuc_detail_file):
        self.nuc_detail_writer = self._create_nuc_writer(nuc_detail_file)
        self.nuc_detail_writer.writeheader()

    def write_counts(self,
                     seed: str,
                     region: str,
                     seed_amino: 'SeedAmino',
                     report_amino: typing.Optional['ReportAmino'],
                     nuc_writer: DictWriter):
        """ Write rows of nucleotide counts for a single codon.

        :param seed: the seed reference name
        :param region: the coordinate reference name
        :param seed_amino: the SeedAmino object for this position
        :param ReportAmino? report_amino: the ReportAmino object for this position
        :param nuc_writer: the CSV writer to write the row into
        """
        ref_offset = seed_amino.ref_offset
        for i, seed_nuc in enumerate(seed_amino.nucleotides):
            if i < seed_amino.nucleotides_to_skip:
                continue
            if seed_amino.consensus_nuc_index is None:
                query_pos_txt = ''
            else:
                query_pos = i + seed_amino.consensus_nuc_index + 1
                query_pos_txt = str(query_pos)
            ref_pos = (str(i + 3*report_amino.position - 2 + ref_offset)
                       if report_amino is not None
                       else '')
            row = {'seed': seed,
                   'region': region,
                   'q-cutoff': self.qcut,
                   'query.nuc.pos': query_pos_txt,
                   'refseq.nuc.pos': ref_pos,
                   'del': seed_nuc.counts['-'],
                   'ins': seed_nuc.insertion_count,
                   'clip': seed_nuc.clip_count,
                   'v3_overlap': seed_nuc.v3_overlap,
                   'coverage': seed_nuc.get_coverage()}
            for base in 'ACTGN':
                nuc_count = seed_nuc.counts[base]
                row[base] = nuc_count
            nuc_writer.writerow(row)

    def merge_extra_counts(self):
        for region, report_aminos in self.reports.items():
            first_amino_index = None
            last_amino_index = None
            last_consensus_nuc_index = None
            for i, report_amino in enumerate(report_aminos):
                seed_amino = report_amino.seed_amino
                if seed_amino.consensus_nuc_index is not None:
                    if first_amino_index is None:
                        first_amino_index = i
                    last_amino_index = i
                    last_consensus_nuc_index = seed_amino.consensus_nuc_index
            if first_amino_index is not None:
                if not self.has_detail_counts:
                    seed_name = self.seed
                else:
                    seed_name = self.detail_seed
                seed_clipping = self.clipping_counts[seed_name]
                seed_insertion_counts = self.conseq_insertion_counts[seed_name]
                for i, report_amino in enumerate(report_aminos):
                    seed_amino = report_amino.seed_amino
                    if i < first_amino_index:
                        consensus_nuc_index = (i - first_amino_index) * 3
                    elif i > last_amino_index:
                        consensus_nuc_index = (i - last_amino_index) * 3 + \
                                              last_consensus_nuc_index
                    else:
                        consensus_nuc_index = seed_amino.consensus_nuc_index
                        if consensus_nuc_index is None:
                            continue
                    for j, seed_nuc in enumerate(seed_amino.nucleotides):
                        query_pos = j + consensus_nuc_index + 1
                        seed_nuc.clip_count = seed_clipping[query_pos]
                        seed_nuc.insertion_count = seed_insertion_counts[query_pos]
                        if j == 2:
                            insertion_counts = self.insert_writer.insert_pos_counts[
                                (self.seed, region)]
                            seed_nuc.insertion_count += insertion_counts[report_amino.position]

    def write_nuc_detail_counts(self, nuc_detail_writer=None):
        nuc_detail_writer = nuc_detail_writer or self.nuc_detail_writer
        self.write_nuc_report(nuc_detail_writer, self.reports, self.detail_seed)

    def write_nuc_counts(self, nuc_writer=None):
        nuc_writer = nuc_writer or self.nuc_writer

        if self.combined_reports:
            for seed, reports in self.combined_reports.items():
                self.write_nuc_report(nuc_writer, reports, seed)
        elif self.reports:
            self.write_nuc_report(nuc_writer, self.reports, self.seed)

    def write_nuc_report(self, nuc_writer, reports, seed):
        self.merge_extra_counts()
        if not self.coordinate_refs:
            for seed_amino in self.seed_aminos[0]:
                self.write_counts(seed, seed, seed_amino, None, nuc_writer)
        else:
            for region, report_aminos in sorted(reports.items()):
                start_ref_offset = None
                for report_amino in report_aminos:
                    # If the region starts with an offset, we reverse it.
                    if start_ref_offset is None:
                        start_ref_offset = report_amino.seed_amino.ref_offset
                    report_amino.seed_amino.ref_offset -= start_ref_offset
                    self.write_counts(seed,
                                      region,
                                      report_amino.seed_amino,
                                      report_amino,
                                      nuc_writer)
                    report_amino.seed_amino.ref_offset += start_ref_offset

    @staticmethod
    def _create_consensus_writer(
            conseq_file,
            include_seed=False,
            include_seed_region_offsets=False,
    ):
        columns = ["region", "q-cutoff", "consensus-percent-cutoff"]
        offsets = ["offset"]
        if include_seed_region_offsets:
            offsets = ["seed-offset", "region-offset"]
        columns = columns + offsets
        columns.append("sequence")
        if include_seed:
            columns = ["seed"] + columns
        return csv.DictWriter(conseq_file, columns, lineterminator=os.linesep)

    def get_consensus_rows(
            self,
            seed_amino_entries,
            ignore_coverage=False,
    ):
        mixture_cutoffs = ([MAX_CUTOFF] if ignore_coverage
                           else self.conseq_mixture_cutoffs)
        min_coverage = 1 if ignore_coverage else self.consensus_min_coverage
        for mixture_cutoff in mixture_cutoffs:
            consensus = ''
            seed_offset = None
            region_offset = None
            for seed_pos, region_pos, seed_amino in seed_amino_entries:
                for nuc_index, seed_nuc in enumerate(seed_amino.nucleotides):
                    nuc_coverage = seed_nuc.get_coverage()
                    if nuc_coverage < min_coverage:
                        if seed_offset is not None:
                            consensus += 'x'
                    else:
                        nuc_consensus = seed_nuc.get_consensus(mixture_cutoff)
                        if seed_offset is None and nuc_consensus:
                            seed_offset = seed_pos + nuc_index
                            if region_pos is not None:
                                region_offset = region_pos + nuc_index
                        consensus += nuc_consensus
                if seed_offset is None:
                    # Still haven't started, so reset the consensus.
                    consensus = ''
            if seed_offset is not None:
                consensus = consensus.rstrip('x')
                yield {
                    'q-cutoff': self.qcut,
                    'consensus-percent-cutoff': format_cutoff(mixture_cutoff),
                    'seed-offset': seed_offset,
                    'region-offset': region_offset,
                    'sequence': consensus
                }

    def _write_consensus_helper(
            self,
            amino_entries,
            csv_writer,
            row_metadata,
            ignore_coverage=False,
            include_seed_region_offsets=False,
    ):
        for row in self.get_consensus_rows(
                amino_entries,
                ignore_coverage=ignore_coverage
        ):
            row.update(row_metadata)
            if not include_seed_region_offsets:
                row.pop("region-offset")
                row["offset"] = row.pop("seed-offset")
            csv_writer.writerow(row)

    def write_consensus_header(self, conseq_file):
        self.conseq_writer = self._create_consensus_writer(conseq_file)
        self.conseq_writer.writeheader()

    def write_consensus(self, conseq_writer=None):
        conseq_writer = conseq_writer or self.conseq_writer
        seed_amino_entries = [
            (seed_amino.consensus_nuc_index, None, seed_amino)
            for seed_amino in self.seed_aminos[0]
        ]
        self._write_consensus_helper(
            seed_amino_entries,
            conseq_writer,
            {"region": self.detail_seed},
        )

    def write_consensus_all_header(self, conseq_all_file):
        self.conseq_all_writer = self._create_consensus_writer(
            conseq_all_file,
            include_seed=True,
            include_seed_region_offsets=True,
        )
        self.conseq_all_writer.writeheader()

    def write_consensus_all(self, conseq_all_writer=None):
        conseq_all_writer = conseq_all_writer or self.conseq_all_writer
        seed_amino_entries = [
            (seed_amino.consensus_nuc_index, None, seed_amino)
            for seed_amino in self.seed_aminos[0]
        ]
        self._write_consensus_helper(
            seed_amino_entries,
            conseq_all_writer,
            {
                "seed": self.detail_seed,
                "region": None,
            },
            ignore_coverage=True,
            include_seed_region_offsets=True,
        )

        regions = sorted(self.reports.keys())
        for region in regions:
            region_amino_entries = [
                (
                    report_amino.seed_amino.consensus_nuc_index,
                    3 * (report_amino.position - 1),
                    report_amino.seed_amino
                )
                for report_amino in self.reports[region]
            ]
            self._write_consensus_helper(
                region_amino_entries,
                conseq_all_writer,
                {
                    "seed": self.seed,
                    "region": region,
                },
                ignore_coverage=True,
                include_seed_region_offsets=True,
            )

    def write_consensus_regions_header(self, conseq_region_file):
        self.conseq_region_writer = self._create_consensus_writer(
            conseq_region_file,
            include_seed=True
        )
        self.conseq_region_writer.writeheader()

    def write_consensus_regions(self, conseq_region_writer=None):
        conseq_region_writer = conseq_region_writer or self.conseq_region_writer
        regions = sorted(self.reports.keys())
        for region in regions:
            region_amino_entries = [
                (3 * (report_amino.position-1), None, report_amino.seed_amino)
                for report_amino in self.reports[region]
            ]
            self._write_consensus_helper(
                region_amino_entries,
                conseq_region_writer,
                {
                    "seed": self.seed,
                    "region": region,
                },
            )

    @staticmethod
    def _create_failure_writer(fail_file):
        return csv.DictWriter(fail_file,
                              ['seed',
                               'region',
                               'qcut',
                               'queryseq',
                               'refseq'],
                              lineterminator=os.linesep)

    def write_failure_header(self, fail_file):
        self.fail_writer = self._create_failure_writer(fail_file)
        self.fail_writer.writeheader()

    def write_failure(self, fail_writer=None):
        fail_writer = fail_writer or self.fail_writer
        if any(self.reports.values()):
            return  # Successfully aligned at least one coordinate reference.
        for region, report_aminos in sorted(self.reports.items()):
            coordinate_ref = self.projects.getReference(region)
            fail_writer.writerow(dict(seed=self.seed,
                                      region=region,
                                      qcut=self.qcut,
                                      queryseq=self.consensus[region],
                                      refseq=coordinate_ref))

    def write_insertions(self, insert_writer=None):
        insert_writer = insert_writer or self.insert_writer
        for coordinate_name, coordinate_inserts in self.inserts.items():
            insert_writer.write(coordinate_inserts,
                                coordinate_name,
                                self.reports[coordinate_name])

    def read_remap_conseqs(self, remap_conseq_csv):
        # noinspection PyTypeChecker
        self.remap_conseqs = dict(map(itemgetter('region', 'sequence'),
                                      csv.DictReader(remap_conseq_csv)))

    def read_contigs(self, contigs_csv):
        self.contigs = list(map(itemgetter('ref', 'group_ref', 'contig'),
                                csv.DictReader(contigs_csv)))

    @staticmethod
    def group_deletions(positions):
        """ Group deletion positions into groups that are within 13 bases.

        They also have to have at least 13 bases without deletions in between
        the groups.
        """
        if not positions:
            return
        gap_size = 13
        working_positions = positions[:]
        # Avoid special check for end condition by adding an extra group.
        working_positions.append(positions[-1] + gap_size)
        next_yield = 0
        for i, pos in enumerate(working_positions):
            if i == 0:
                continue
            start = positions[next_yield]
            end = positions[i-1]
            if pos - end >= gap_size:
                if end - start < gap_size and (i - next_yield) % 3 == 0:
                    yield positions[next_yield:i]
                next_yield = i

    def align_deletions(self, aligned_reads):
        """ Align codon deletions to the codon boundaries.

        :param aligned_reads: group of rows from the CSV reader of the aligned
            reads.
        :return: a generator of rows from the CSV reader with codon deletions
            aligned to the codon boundaries.
        """
        if self.remap_conseqs is None:
            yield from aligned_reads
        reading_frames = None
        for row in aligned_reads:
            if reading_frames is None:
                reading_frames = self.load_reading_frames(row['refname'])
            seq = row['seq']
            seq_offset = int(row['offset'])
            for old_positions in self.group_deletions(
                    [match.start() for match in re.finditer('-', seq)]):
                chars = list(seq)
                anchor_index = len(old_positions)//2
                anchor_pos = old_positions[anchor_index]
                start_pos = anchor_pos-len(old_positions)//2
                end_pos = start_pos + len(old_positions)
                nuc_pos = seq_offset + start_pos
                reading_frame = reading_frames[nuc_pos]
                offset = (nuc_pos + reading_frame) % 3
                if offset == 1:
                    start_pos -= 1
                    end_pos -= 1
                elif offset == 2:
                    start_pos += 1
                    end_pos += 1
                new_positions = list(range(start_pos, end_pos))
                for pos in reversed(old_positions):
                    del chars[pos]
                for pos in new_positions:
                    chars.insert(pos, '-')
                seq = ''.join(chars)
                row['seq'] = seq
            yield row

    def load_reading_frames(self, seed_name):
        """ Calculate reading frames along a consensus sequence.

        :param seed_name: the name of the seed to look up
        :return: {pos: frame} zero-based position and reading frame for each
            position. Frame 1 needs one nucleotide inserted at start.
        """
        result = Counter()
        conseq = self.remap_conseqs[seed_name]
        coord_refs = self.projects.getCoordinateReferences(seed_name)
        if not coord_refs:
            _, seed_name = seed_name.split('-', 1)
            coord_refs = self.projects.getCoordinateReferences(seed_name)
        for coord_ref in coord_refs.values():
            best_alignment = (-1000000, '', '', 0)
            for frame_index in range(3):
                conseq_aminos = translate('-'*frame_index + conseq)
                aconseq, acoord, score = self._pair_align(conseq_aminos,
                                                          coord_ref,
                                                          GAP_OPEN_COORD,
                                                          GAP_EXTEND_COORD)
                best_alignment = max(best_alignment, (score, aconseq, acoord, frame_index))
            score, aconseq, acoord, frame_index = best_alignment
            if frame_index == 0:
                continue  # defaults to 0, no need to record
            conseq_codon_index = -1
            coord_codon_index = -1
            for conseq_amino, coord_amino in zip(aconseq, acoord):
                if conseq_amino != '-':
                    conseq_codon_index += 1
                if coord_amino == '-':
                    continue
                coord_codon_index += 1
                
                nuc_pos = conseq_codon_index * 3 - frame_index
                for i in range(3):
                    result[nuc_pos+i] = frame_index
        return result


class SeedAmino(object):
    """
    Records the frequencies of amino acids at a given position of the
    aligned reads as determined by the consensus sequence.
    """
    def __init__(self, consensus_nuc_index, counts=None):
        self.v3_overlap = 0
        self.consensus_nuc_index = consensus_nuc_index
        self.counts = counts or Counter()  # {amino: count}
        self.codon_counts = Counter()  # {codon_nucs: count}
        self.nucleotides = [SeedNucleotide() for _ in range(3)]
        self.low_quality = 0
        self.partial = 0
        self.deletions = 0
        self.read_count = 0
        self.ref_offset = 0
        self.nucleotides_to_skip = 0

    def __repr__(self):
        if self.counts:
            return 'SeedAmino({!r}, {!r})'.format(self.consensus_nuc_index,
                                                  dict(self.counts))
        return 'SeedAmino({})'.format(self.consensus_nuc_index)

    def count_aminos(self, codon_seq, count):
        """ Record a set of reads at this position in the seed reference.
        @param codon_seq: a string of three nucleotides that were read at this
                          position, may be padded with spaces at the start
                          or end of a sequence, or dashes for deletions
        @param count: the number of times they were read
        """
        self.read_count += count
        self.codon_counts[codon_seq] += count
        if 'N' in codon_seq:
            self.low_quality += count
        elif '---' == codon_seq:
            self.deletions += count
        elif '-' in codon_seq:
            self.partial += count  # Partial deletion
        elif ' ' not in codon_seq and 'n' not in codon_seq:
            amino = translate(codon_seq.upper())
            self.counts[amino] += count
        elif 'nnn' == codon_seq:
            # Don't count the gap between forward and reverse reads in a pair.
            self.read_count -= count
        for i, nuc in enumerate(codon_seq):
            if nuc != ' ':
                seed_nucleotide = self.nucleotides[i]
                seed_nucleotide.count_nucleotides(nuc, count)

    def add(self, other: 'SeedAmino'):
        self.counts += other.counts
        self.partial += other.partial
        self.deletions += other.deletions
        self.read_count += other.read_count
        self.low_quality += other.low_quality
        self.nucleotides_to_skip = other.nucleotides_to_skip
        self.ref_offset = other.ref_offset
        for nuc, other_nuc in zip(self.nucleotides, other.nucleotides):
            nuc.add(other_nuc)

    def get_report(self) -> str:
        """ Build a report string with the counts of each amino acid.

        Report how many times each amino acid was seen in count_aminos().
        @return: comma-separated list of counts in the same order as the
        AMINO_ALPHABET list
        """
        return ','.join([str(self.counts[amino])
                         for amino in AMINO_ALPHABET])

    def apply_repeat(self, repeated_nuc: int) -> 'SeedAmino':
        new_amino = SeedAmino(self.consensus_nuc_index)
        for codon, count in self.codon_counts.items():
            new_codon = codon[:repeated_nuc + 1] + codon[repeated_nuc:2]
            new_amino.count_aminos(new_codon, count)
        return new_amino

    def get_consensus(self) -> str:
        """ Find the amino acid that was seen most often in count_aminos().

        If there is a tie, just pick one of the tied amino acids.
        @return: the letter of the most common amino acid
        """
        consensus = self.counts.most_common(1)
        if consensus:
            return consensus[0][0]
        if self.read_count:
            return '?'
        return '-'

    def count_overlap(self, other):
        for nuc1, nuc2 in zip(self.nucleotides, other.nucleotides):
            nuc1.count_overlap(nuc2)
            self.v3_overlap = max(self.v3_overlap, nuc1.v3_overlap)


class SeedNucleotide(object):
    """
    Records the frequencies of nucleotides at a given position of the
    aligned reads as determined by the consensus sequence.
    """
    COUNTED_NUCS = 'ACTG-'

    def __init__(self, counts=None):
        self.v3_overlap = self.clip_count = self.insertion_count = 0
        self.counts = counts or Counter()

    def __repr__(self):
        return 'SeedNucleotide({!r})'.format(dict(self.counts))

    def count_nucleotides(self, nuc_seq, count=1):
        """ Record a set of reads at this position in the seed reference.
        @param nuc_seq: a single nucleotide letter that was read at this
        position
        @param count: the number of times it was read
        """
        if nuc_seq == 'n':
            "Represents gap between forward and reverse read, ignore."
        else:
            self.counts[nuc_seq] += count

    def add(self, other):
        self.counts += other.counts
        self.clip_count += other.clip_count
        self.insertion_count += other.insertion_count

    def get_report(self):
        """ Build a report string with the counts of each nucleotide.

        Report how many times each nucleotide was seen in count_nucleotides().
        @return: comma-separated list of counts for A, C, G, and T.
        """
        return ','.join(map(str, [self.counts[nuc] for nuc in 'ACGT']))

    def get_coverage(self):
        return sum(self.counts[nuc] for nuc in self.COUNTED_NUCS)

    def get_consensus(self, mixture_cutoff, no_coverage=''):
        """ Choose consensus nucleotide or mixture from the counts.

        @param mixture_cutoff: the minimum fraction of reads
            that a nucleotide must be found in for it to be considered,
            or MAX_CUTOFF to consider only the most common nucleotide.
        @param no_coverage: what to return when there are no reads mapped to
            this position.
        @return: The letter for the consensus nucleotide or mixture.
            Nucleotide mixtures are encoded by IUPAC symbols, and the most common
            nucleotide can be a mixture if there is a tie.
        """
        if not self.counts:
            return no_coverage

        coverage = self.get_coverage()
        if mixture_cutoff not in (MAX_CUTOFF, FIRST_CUTOFF):
            min_count = coverage * mixture_cutoff
        else:
            min_count = 0
        mixture = []
        for nuc, count in self.counts.most_common():
            if count < min_count:
                break
            if nuc in self.COUNTED_NUCS:
                mixture.append(nuc)
                if mixture_cutoff in (MAX_CUTOFF, FIRST_CUTOFF):
                    # Catch any ties before breaking out.
                    min_count = count

        has_deletion = '-' in mixture
        if has_deletion:
            mixture.remove('-')
        if len(mixture) > 1:
            mixture.sort()
            if mixture_cutoff == FIRST_CUTOFF:
                consensus = mixture[0]
            else:
                consensus = ambig_dict[''.join(mixture)]
        elif len(mixture) == 1:
            # no ambiguity
            consensus = mixture[0]
        else:
            # Nothing left to go in the mixture.
            consensus = '-' if has_deletion else 'N'
        if has_deletion:
            consensus = consensus.lower()
        return consensus

    def count_overlap(self, other):
        for nuc in 'ACGT':
            self.v3_overlap += other.counts[nuc]


class ReportAmino(object):
    def __init__(self, seed_amino, position):
        """ Create a new instance.

        @param seed_amino: Counts for the
        """
        self.seed_amino = seed_amino
        self.position = position
        self.max_clip_count = 0
        self.insertion_count = 0

    def __repr__(self):
        return 'ReportAmino({!r}, {})'.format(self.seed_amino, self.position)


class InsertionWriter(object):
    def __init__(self, insert_file):
        """ Initialize a writer object.

        @param insert_file: an open file that the data will be written to
        """
        self.insert_file = insert_file
        self.insert_writer = csv.DictWriter(insert_file,
                                            ['seed',
                                             'region',
                                             'qcut',
                                             'left',
                                             'insert',
                                             'count',
                                             'before'],
                                            lineterminator=os.linesep)
        self.insert_writer.writeheader()

        # {(seed, region): {pos: insert_count}}
        self.insert_pos_counts = defaultdict(Counter)
        self.seed = self.qcut = None
        self.insert_file_name = getattr(insert_file, 'name', None)
        if self.insert_file_name is None or self.insert_file_name == os.devnull:
            self.nuc_seqs = Counter()
            self.nuc_seqs_context = None
        else:
            dirname = os.path.dirname(self.insert_file_name)
            file_prefix = os.path.join(os.path.abspath(dirname),
                                       'nuc_read_counts')
            self.nuc_seqs = self.nuc_seqs_context = BigCounter(
                file_prefix=file_prefix)

    def __enter__(self):
        """ Context manager makes sure that BigCounter gets cleaned up. """
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        if self.nuc_seqs_context is not None:
            self.nuc_seqs_context.__exit__(exc_type, exc_val, exc_tb)

    def start_group(self, seed, qcut):
        """ Start a new group of reads.

        @param seed: the name of the region these reads mapped to
        @param qcut: the quality cut off used for these reads
        """
        self.seed = seed
        self.qcut = qcut
        self.nuc_seqs.clear()

    def add_nuc_read(self, offset_sequence, count):
        """ Add a read to the group.

        @param offset_sequence: the nucleotide sequence of the read that has
            had dashes added to offset it into the consensus sequence
            coordinates
        @param count: the number of times this sequence was read
        """
        self.nuc_seqs[offset_sequence] += count

    def write(self, inserts, region, report_aminos=None):
        """ Write any insert ranges to the file.

        Sequence data comes from the reads that were added to the current group.
        @param inserts: indexes of positions in the reads that should be
            reported as insertions.
        @param region: the name of the coordinate region the current group was
            mapped to
        @param report_aminos: a list of ReportAmino objects that represent the
            sequence that successfully mapped to the coordinate reference.
        """
        if len(inserts) == 0:
            return

        report_aminos = report_aminos or []

        region_insert_pos_counts = self.insert_pos_counts[(self.seed, region)]
        inserts = list(inserts)
        inserts.sort()

        # convert insertion coordinates into contiguous ranges
        insert_ranges = []
        for insert in inserts:
            if not insert_ranges or insert != insert_ranges[-1][1]:
                # just starting or we hit a gap
                insert_ranges.append([insert, insert + 3])
            else:
                insert_ranges[-1][1] += 3

        # enumerate insertions by popping out all AA sub-string variants
        insert_counts = OrderedDict()  # {left: {insert_seq: count}}
        insert_targets = {}  # {left: inserted_before_pos}
        for left, right in insert_ranges:
            for report_amino in report_aminos:
                seed_amino = report_amino.seed_amino
                if seed_amino.consensus_nuc_index == right:
                    insert_targets[left] = report_amino.position
                    break
            current_counts = Counter()
            insert_counts[left] = current_counts
            for nuc_seq, count in self.nuc_seqs.items():
                insert_nuc_seq = nuc_seq[left:right]
                is_valid = (insert_nuc_seq and
                            'n' not in insert_nuc_seq and
                            '-' not in insert_nuc_seq)
                if is_valid:
                    insert_amino_seq = translate(insert_nuc_seq)
                    if insert_amino_seq:
                        current_counts[insert_amino_seq] += count

        # record insertions to CSV
        for left, counts in insert_counts.items():
            for insert_seq, count in counts.most_common():
                insert_before = insert_targets.get(left)
                # Only care about insertions in the middle of the sequence,
                # so ignore any that come before or after the reference.
                # Also report if we're in test mode (no report_aminos).
                if not report_aminos or insert_before not in (1, None):
                    row = dict(seed=self.seed,
                               region=region,
                               qcut=self.qcut,
                               left=left + 1,
                               insert=insert_seq,
                               count=count,
                               before=insert_before)
                    self.insert_writer.writerow(row)
                    if insert_before is not None:
                        region_insert_pos_counts[insert_before-1] += count


def format_cutoff(cutoff):
    """ Format the cutoff fraction as a string to use as a name. """

    if cutoff == MAX_CUTOFF:
        return cutoff
    return '{:0.3f}'.format(cutoff)


def aln2counts(aligned_csv,
               nuc_csv,
               amino_csv,
               coord_ins_csv,
               conseq_csv,
               failed_align_csv,
               callback=None,
               coverage_summary_csv=None,
               clipping_csv=None,
               conseq_ins_csv=None,
               g2p_aligned_csv=None,
               remap_conseq_csv=None,
               conseq_region_csv=None,
               amino_detail_csv=None,
               genome_coverage_csv=None,
               nuc_detail_csv=None,
               contigs_csv=None,
               conseq_all_csv=None,
               minimap_hits_csv=None):
    """
    Analyze aligned reads for nucleotide and amino acid frequencies.
    Generate consensus sequences.
    @param aligned_csv:         Open file handle containing aligned reads (from sam2aln)
    @param nuc_csv:             Open file handle to write nucleotide frequencies.
    @param amino_csv:           Open file handle to write amino acid frequencies.
    @param coord_ins_csv:       Open file handle to write insertions relative to coordinate reference.
    @param conseq_csv:          Open file handle to write consensus sequences.
    @param failed_align_csv:    Open file handle to write sample consensus sequences that failed to
                                align to the coordinate reference.
    @param callback: a function to report progress with three optional
        parameters - callback(message, progress, max_progress)
    @param coverage_summary_csv: Open file handle to write coverage depth.
    @param clipping_csv: Open file handle containing soft clipping counts
    @param conseq_ins_csv: Open file handle containing insertions relative to
        consensus sequence
    @param g2p_aligned_csv: Open file handle containing aligned reads (from
        fastq_g2p)
    @param remap_conseq_csv: Open file handle containing consensus sequences
        from the remap step.
    @param conseq_region_csv: Open file handle to write consensus sequences
        split into regions.
    @param amino_detail_csv: Open file handle to write amino acid frequencies
        for individual contigs.
    @param nuc_detail_csv: Open file handle to write nucleotide frequencies
        for individual contigs.
    @param genome_coverage_csv: Open file handle to write coverage for individual
        contigs.
    @param contigs_csv: Open file handle to read contig sequences.
    @param conseq_all_csv: Open file handle to write consensus sequences *ignoring
        inadequate coverage*.
    @param minimap_hits_csv: Open file handle to write minimap2 match locations.
    """
    # load project information
    projects = ProjectConfig.loadDefault()

    landmarks_path = (Path(__file__).parent.parent / 'data' /
                      'landmark_references.yaml')
    landmarks_yaml = landmarks_path.read_text()

    # initialize reporter classes
    with InsertionWriter(coord_ins_csv) as insert_writer:
        report = SequenceReport(insert_writer,
                                projects,
                                CONSEQ_MIXTURE_CUTOFFS,
                                landmarks_yaml=landmarks_yaml)
        report.consensus_min_coverage = CONSENSUS_MIN_COVERAGE
        report.write_amino_header(amino_csv)
        report.write_consensus_header(conseq_csv)
        report.write_consensus_regions_header(conseq_region_csv)
        report.write_failure_header(failed_align_csv)
        report.write_nuc_header(nuc_csv)
        if coverage_summary_csv is None:
            coverage_summary = coverage_writer = None
        else:
            coverage_writer = csv.DictWriter(coverage_summary_csv,
                                             ['avg_coverage',
                                              'coverage_region',
                                              'region_width'],
                                             lineterminator=os.linesep)
            coverage_writer.writeheader()
            coverage_summary = {}
        if nuc_detail_csv is not None:
            report.write_nuc_detail_header(nuc_detail_csv)
        if amino_detail_csv is not None:
            report.write_amino_detail_header(amino_detail_csv)

        if callback:
            aligned_filename = getattr(aligned_csv, 'name', None)
            if aligned_filename:
                file_size = os.stat(aligned_filename).st_size
                report.enable_callback(callback, file_size)

        if conseq_all_csv is not None:
            report.write_consensus_all_header(conseq_all_csv)
        if clipping_csv is not None:
            report.read_clipping(clipping_csv)
        if conseq_ins_csv is not None:
            report.read_insertions(conseq_ins_csv)
        if remap_conseq_csv is not None:
            report.read_remap_conseqs(remap_conseq_csv)
        if contigs_csv is not None:
            report.read_contigs(contigs_csv)
        if genome_coverage_csv is not None:
            report.write_genome_coverage_header(genome_coverage_csv)
        if minimap_hits_csv is not None:
            report.write_minimap_hits_header(minimap_hits_csv)

        report.process_reads(aligned_csv,
                             coverage_summary,
                             excluded_regions={'V3LOOP'})
        if g2p_aligned_csv is not None:
            report.nuc_detail_writer = report.amino_detail_writer = None
            report.combined_reports.clear()
            if report.remap_conseqs is not None:
                report.remap_conseqs[G2P_SEED_NAME] = projects.getReference(
                    G2P_SEED_NAME)
            report.process_reads(g2p_aligned_csv,
                                 coverage_summary,
                                 included_regions={'V3LOOP'})

        if coverage_summary_csv is not None:
            if coverage_summary:
                coverage_writer.writerow(coverage_summary)


def main():
    args = parse_args()
    aln2counts(args.aligned_csv,
               args.nuc_csv,
               args.amino_csv,
               args.coord_ins_csv,
               args.conseq_csv,
               args.failed_align_csv,
               coverage_summary_csv=args.coverage_summary_csv,
               clipping_csv=args.clipping_csv,
               conseq_ins_csv=args.conseq_ins_csv,
               g2p_aligned_csv=args.g2p_aligned_csv,
               remap_conseq_csv=args.remap_conseq_csv,
               conseq_region_csv=args.conseq_region_csv,
               amino_detail_csv=args.amino_detail_csv,
               nuc_detail_csv=args.nuc_detail_csv,
               genome_coverage_csv=args.genome_coverage_csv,
               contigs_csv=args.contigs_csv,
               conseq_all_csv=args.conseq_all_csv,
               minimap_hits_csv=args.minimap_hits_csv)


if __name__ == '__main__':
    main()
