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
from collections import Counter, defaultdict
import csv
from csv import DictWriter
from itertools import groupby, chain
from operator import itemgetter
import os
from pathlib import Path
import logging

import gotoh
import yaml

from micall.core.project_config import ProjectConfig, G2P_SEED_NAME
from micall.core.remap import PARTIAL_CONTIG_SUFFIX, REVERSED_CONTIG_SUFFIX
from micall.data.landmark_reader import LandmarkReader
from micall.utils.big_counter import BigCounter
from micall.utils.consensus_aligner import CigarActions, ConsensusAligner
from micall.utils.report_amino import ReportAmino, MAX_CUTOFF, SeedAmino, AMINO_ALPHABET, ReportNucleotide, \
    SeedNucleotide
from micall.utils.spring_beads import Wire, Bead
from micall.utils.translation import translate

logger = logging.getLogger(__name__)

CONSEQ_MIXTURE_CUTOFFS = [0.01, 0.02, 0.05, 0.1, 0.2, 0.25]
GAP_OPEN_COORD = 40
GAP_EXTEND_COORD = 10
CONSENSUS_MIN_COVERAGE = 100


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
    parser.add_argument('--insertions_csv',
                        type=argparse.FileType('w'),
                        default=os.devnull,
                        help='CSV containing insertions relative to coordinate reference (new)')
    parser.add_argument('--conseq_csv',
                        type=argparse.FileType('w'),
                        default=os.devnull,
                        help='CSV containing consensus sequences')
    parser.add_argument('--conseq_all_csv',
                        type=argparse.FileType('w'),
                        help='CSV containing consensus sequences (ignoring inadequate '
                             'coverage)')
    parser.add_argument('--conseq_stitched_csv',
                        type=argparse.FileType('w'),
                        help='CSV containing stitched whole genome consensus sequences')
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


def get_insertion_info(right, report_aminos, report_nucleotides):
    insert_behind = None
    insertion_coverage = 0
    for i, report_amino in enumerate(report_aminos):
        seed_amino = report_amino.seed_amino
        if seed_amino.consensus_nuc_index == right:
            insert_behind = (report_amino.position - 1) * 3
            coverage_left = report_aminos[i - 1].seed_amino.nucleotides[2].get_coverage() if i >= 0 else 0
            coverage_right = report_amino.seed_amino.nucleotides[0].get_coverage()
            insertion_coverage = max(coverage_left, coverage_right)
            break
    if len(report_aminos) == 0:
        for i, report_nuc in enumerate(report_nucleotides):
            seed_nuc = report_nuc.seed_nucleotide
            if seed_nuc.consensus_index == right:
                insert_behind = report_nuc.position - 1
                coverage_left = report_nucleotides[i - 1].seed_nucleotide.get_coverage() if i >= 0 else 0
                coverage_right = report_nuc.seed_nucleotide.get_coverage()
                insertion_coverage = max(coverage_left, coverage_right)
                break
    return insert_behind, insertion_coverage


def combine_region_nucleotides(nuc_dict, region_nucleotides, region_start):
    assert region_start is not None
    for nuc_index, nucleotide in enumerate(region_nucleotides):
        position = region_start + nuc_index - 1
        if position not in nuc_dict.keys():
            nuc_dict[position] = nucleotide.seed_nucleotide
        else:
            if len(nuc_dict[position].counts) == 0 and len(nucleotide.seed_nucleotide.counts) != 0:
                logger.debug(f"Zero count in dict at position {position + 1}. Inserting non-zero counts.")
                nuc_dict[position] = nucleotide.seed_nucleotide
            elif len(nuc_dict[position].counts) != 0 and len(nucleotide.seed_nucleotide.counts) == 0:
                logger.debug(f"Zero count in nucleotide at position {position + 1}. Inserting non-zero counts.")
            else:
                if nuc_dict[position].counts != nucleotide.seed_nucleotide.counts:
                    logger.debug(f"Counts don't match up. Position {position + 1}")
                    logger.debug(f"Counts in dict: {nuc_dict[position].counts}")
                    logger.debug(f"Counts in nucleotide: {nucleotide.seed_nucleotide.counts}")
                    logger.debug("Continuing with dict counts.")
    return nuc_dict


def combine_region_insertions(insertions_dict, region_insertions, region_start):
    assert region_start is not None
    if region_insertions is None:
        return insertions_dict
    new_insertions = {}
    for position in region_insertions.keys():
        ref_position = position + region_start - 1
        if ref_position not in insertions_dict:
            shifted_positions = (ref_position-3, ref_position-2, ref_position-1, ref_position+1, ref_position+2,
                                 ref_position+3)
            if any(shift_pos in insertions_dict for shift_pos in shifted_positions):
                logger.debug(f"Disagreement or shift in insertions between regions. Position {ref_position + 1}")
            else:
                new_insertions[ref_position] = region_insertions[position]
        else:
            if insertions_dict[ref_position].items() != region_insertions[position].items():
                logger.debug(f"Insertion counts don't match up. Position {ref_position + 1}")
                logger.debug("Continuing with dict counts.")
    insertions_dict.update(new_insertions)
    return insertions_dict


def insert_insertions(insertions, consensus_nucs):
    for position in insertions.keys():
        neighbour_left = consensus_nucs.get(position)
        neighbour_right = consensus_nucs.get(position + 1)
        coverage_left = neighbour_left.get_coverage() if neighbour_left is not None else 0
        coverage_right = neighbour_right.get_coverage() if neighbour_left is not None else 0
        coverage_nuc = max(coverage_left, coverage_right)
        num_insertions = len(insertions[position])
        for i in range(num_insertions):
            # use non-integer positions for insertion nucleotides
            insertion_position = position + (i + 1) / (num_insertions + 1)
            coverage_insertion = insertions[position][i].get_coverage()
            coverage_no_insertion = coverage_nuc - coverage_insertion
            if coverage_no_insertion > 0:
                insertions[position][i].counts['-'] = coverage_no_insertion
            consensus_nucs[insertion_position] = insertions[position][i]
    return consensus_nucs


def aggregate_insertions(insertions_counter, coverage_nuc=0, consensus_pos=None):
    aggregated_insertions = defaultdict()

    if len(insertions_counter) == 0:
        return aggregated_insertions

    sorted_counts = sorted(insertions_counter.items(), key=lambda item: -len(item[0]))
    length = len(sorted_counts[0][0])
    for i in range(length):
        insertion_nuc = SeedNucleotide()
        insertion_nuc.consensus_index = consensus_pos
        for insertion in sorted_counts:
            if len(insertion[0]) > i:
                insertion_nuc.count_nucleotides(insertion[0][i], insertion[1])
        coverage_insertion = insertion_nuc.get_coverage(count_nucs='ACGT')
        coverage_no_insertion = coverage_nuc - coverage_insertion
        if coverage_no_insertion > 0:
            insertion_nuc.count_nucleotides('-', coverage_no_insertion)
        aggregated_insertions[i] = insertion_nuc

    return aggregated_insertions


class ConsensusBuilder(object):
    """ Helper class to build the consensus for different mixture cutoffs."""
    def __init__(self,
                 mixture_cutoffs,
                 consensus_min_coverage):
        self.conseq_mixture_cutoffs = mixture_cutoffs
        self.consensus_min_coverage = consensus_min_coverage
        self.qcut = None

    def get_consensus_rows(
            self,
            seed_entries,
            ignore_coverage=False,
            is_nucleotide=False,
            discard_deletions=True,
    ):
        mixture_cutoffs = ([MAX_CUTOFF] if ignore_coverage
                           else self.conseq_mixture_cutoffs)
        min_coverage = 1 if ignore_coverage else self.consensus_min_coverage
        for mixture_cutoff in mixture_cutoffs:
            consensus = ''
            seed_offset = None
            if not is_nucleotide:
                for seed_pos, seed_amino in seed_entries:
                    for nuc_index, seed_nuc in enumerate(seed_amino.nucleotides):
                        nuc_coverage = seed_nuc.get_coverage()
                        if nuc_coverage < min_coverage:
                            if seed_offset is not None:
                                consensus += 'x'
                        else:
                            nuc_consensus = seed_nuc.get_consensus(mixture_cutoff, discard_deletions=discard_deletions)
                            if seed_offset is None and nuc_consensus:
                                seed_offset = seed_pos + nuc_index
                            consensus += nuc_consensus
                    if seed_offset is None:
                        # Still haven't started, so reset the consensus.
                        consensus = ''
            else:
                for nuc_index, seed_nuc in seed_entries:
                    nuc_coverage = seed_nuc.get_coverage()
                    if nuc_coverage < min_coverage:
                        if seed_offset is not None:
                            consensus += 'x'
                    else:
                        nuc_consensus = seed_nuc.get_consensus(mixture_cutoff, discard_deletions=discard_deletions)
                        if seed_offset is None and nuc_consensus:
                            seed_offset = nuc_index
                        consensus += nuc_consensus
            if seed_offset is not None:
                consensus = consensus.rstrip('x')
                yield {
                    'q-cutoff': self.qcut,
                    'consensus-percent-cutoff': format_cutoff(mixture_cutoff),
                    'seed-offset': seed_offset,
                    'sequence': consensus
                }


class SequenceReport(object):
    """ Hold the data for several reports related to a sample's genetic sequence.

    To use a report object, read a group of aligned reads that mapped to a
    single region, and then write out all the reports for that region.
    """
    def __init__(self,
                 insert_writer: 'InsertionWriter',
                 projects: ProjectConfig,
                 conseq_mixture_cutoffs: typing.List[float],
                 clipping_counts=None,
                 conseq_insertion_counts=None,
                 landmarks_yaml=None,
                 consensus_min_coverage=0):
        """ Create an object instance.

        :param insert_writer: will track reads and write out any insertions
            relative to the coordinate reference.
        :param projects: details of reference sequences and project codes
        :param conseq_mixture_cutoffs: a list of cutoff fractions used to
            determine what portion a variant must exceed before it will be
            included as a mixture in the consensus.
        """
        self.consensus_min_coverage = consensus_min_coverage
        self.callback_progress = 0
        self.callback_next = self.callback_chunk_size = self.callback_max = None
        self.insert_writer = insert_writer
        self.projects = projects
        self.conseq_mixture_cutoffs = list(conseq_mixture_cutoffs)
        self.conseq_mixture_cutoffs.insert(0, MAX_CUTOFF)
        self.callback = None
        self.seed_aminos = self.reading_frames = None
        self.report_nucleotides: typing.Dict[
            str,  # coordinate_name
            typing.List[ReportNucleotide]] = defaultdict(list)
        self.reports: typing.Dict[
            str,  # coordinate_name
            typing.List[ReportAmino]] = defaultdict(list)
        self.inserts = self.consensus = self.seed = self.insert_nucs = None
        self.contigs = None  # [(ref, group_ref, seq)]
        self.consensus_by_reading_frame = None  # {frame: seq}
        self.consensus_aligner = ConsensusAligner(projects)
        if landmarks_yaml is None:
            landmarks_path = (Path(__file__).parent.parent / 'data' /
                              'landmark_references.yaml')
            landmarks_yaml = landmarks_path.read_text()
        self.landmarks = yaml.safe_load(landmarks_yaml)
        self.coordinate_refs = self.remap_conseqs = None

        self.consensus_builder = ConsensusBuilder(self.conseq_mixture_cutoffs,
                                                  consensus_min_coverage)

        self.combined_reports: typing.Dict[
            str,  # seed
            typing.Dict[str,  # coord_region
                        typing.List[ReportAmino]]] = defaultdict(
            lambda: defaultdict(list))
        self.combined_report_nucleotides: typing.Dict[
            str,  # seed
            typing.Dict[str,  # coord_region
                        typing.List[ReportNucleotide]]] = defaultdict(
            lambda: defaultdict(list))
        self.combined_insertions = defaultdict(lambda:  # reference
                                               defaultdict(lambda:  # region
                                                           defaultdict(lambda:  # ref position
                                                                       defaultdict(lambda:  # insertion position
                                                                                   SeedNucleotide()))))

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
        self.conseq_stitched_writer = None

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
                self.consensus_builder.qcut = row['qcut']
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

    def _map_to_coordinate_ref(self, coordinate_name):
        """ Extract the coordinate region from the aligned reads.

        Consensus has already been aligned against the full-genome nucleotide
        sequence. Populates self.reports with ReportAmino objects and
        self.report_nucleotides with ReportNucleotide objects.

        :param coordinate_name: Name of coordinate reference.
        """

        is_amino_coordinate = self.projects.isAmino(coordinate_name)
        landmark_reader = LandmarkReader(self.landmarks)
        region_info = landmark_reader.get_gene(
            self.consensus_aligner.coordinate_name,
            coordinate_name,
            drop_stop_codon=False)
        report_nucleotides = self.report_nucleotides[coordinate_name]
        if is_amino_coordinate:
            report_aminos = self.reports[coordinate_name]
            amino_ref = self.projects.getReference(coordinate_name)
        else:
            report_aminos = amino_ref = None
        repeated_ref_pos = region_info.get('duplicated_pos')
        skipped_ref_pos = region_info.get('skipped_pos')
        self.consensus_aligner.report_region(region_info['start'],
                                             region_info['end'],
                                             report_nucleotides,
                                             report_aminos,
                                             repeated_ref_pos,
                                             skipped_ref_pos,
                                             amino_ref)
        self.consensus[coordinate_name] = self.consensus_aligner.amino_consensus
        self.inserts[coordinate_name] = self.consensus_aligner.inserts

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
        if self.conseq_stitched_writer is not None:
            self.write_whole_genome_consensus_from_nuc(self.conseq_stitched_writer)

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
        self.reports.clear()  # {coord_name: [ReportAmino()]}
        self.report_nucleotides.clear()  # {coord_name: [ReportNucleotide()]
        self.reading_frames = {}  # {coord_name: reading_frame}
        self.inserts = {}  # {coord_name: set([consensus_index])}
        self.insert_nucs = {}  # {coord_name: {position: insertion counts}}
        self.consensus = {}  # {coord_name: consensus_amino_seq}
        if self.consensus_aligner.consensus is not None:
            self.consensus_aligner = ConsensusAligner(self.projects)

        # populates these dictionaries, generates amino acid counts
        self._count_reads(aligned_reads)
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

        self.consensus_aligner.start_contig(coordinate_name,
                                            reading_frames=self.seed_aminos)

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

        # iterate over coordinate references defined for this region
        for coordinate_name in self.coordinate_refs:
            if included_regions is not None and \
                    coordinate_name not in included_regions:
                continue
            if coordinate_name in excluded_regions:
                continue
            self._map_to_coordinate_ref(coordinate_name)

    def read_clipping(self, clipping_csv):
        for row in csv.DictReader(clipping_csv):
            pos = int(row['pos'])
            count = int(row['count'])
            self.clipping_counts[row['refname']][pos] += count

    def read_insertions(self, conseq_ins_csv):
        reader = csv.DictReader(conseq_ins_csv)

        # {ref: {pos: set([qname])}}
        insertion_names = defaultdict(lambda: defaultdict(set))
        insertion_nucs = defaultdict(lambda: defaultdict(Counter))

        for row in reader:
            ref_name = row['refname']
            pos = int(row['pos'])
            pos_insertions = insertion_nucs[ref_name][pos]
            pos_names = insertion_names[ref_name][pos]
            if row['qname'] not in pos_names:
                pos_insertions.update([row['insert']])
            pos_names.add(row['qname'])
        for ref_name, ref_positions in insertion_names.items():
            self.conseq_insertion_counts[ref_name] = Counter(
                {pos: len(names) for pos, names in ref_positions.items()})
        for ref_name in insertion_nucs.keys():
            for ref_position in insertion_nucs[ref_name].keys():
                self.insert_writer.conseq_insertions[ref_name][ref_position] = \
                    aggregate_insertions(insertion_nucs[ref_name][ref_position], consensus_pos=ref_position - 1)

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
                                                  len(old_report) + 1))
                old_report[i].seed_amino.add(report_amino.seed_amino)
        old_report_nucs = self.combined_report_nucleotides[group_ref]
        for region, report in self.report_nucleotides.items():
            old_report_nuc = old_report_nucs[region]
            for i, report_nuc in enumerate(report):
                while len(old_report_nuc) <= i:
                    old_report_nuc.append(ReportNucleotide(len(old_report_nuc)+1))
                old_report_nuc[i].seed_nucleotide.add(report_nuc.seed_nucleotide)
        old_insertions = self.combined_insertions[group_ref]
        for region in self.insert_writer.ref_insertions:
            if self.insert_writer.ref_insertions[region] is not None:
                for position in self.insert_writer.ref_insertions[region]:
                    for insertion_position in self.insert_writer.ref_insertions[region][position]:
                        nuc = old_insertions[region][position].get(insertion_position)
                        if nuc is not None:
                            nuc.add(self.insert_writer.ref_insertions[region][position][insertion_position])
                        else:
                            old_insertions[region][position][insertion_position] =\
                                self.insert_writer.ref_insertions[region][position][insertion_position]

        self.reports.clear()
        self.report_nucleotides.clear()
        self.insert_writer.ref_insertions.clear()

    def write_amino_report(self,
                           amino_writer: DictWriter,
                           reports: typing.Dict[str, typing.List['ReportAmino']],
                           seed: str,
                           coverage_summary: dict = None):
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
                pos_count += 1
                for field_name in ('coverage',
                                   'clip',
                                   'v3_overlap',
                                   'ins',
                                   'del',
                                   'partial',
                                   'X'):
                    if row[field_name]:
                        break
                else:
                    # Nothing useful, don't write this row.
                    continue
                amino_writer.writerow(row)
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
        self.write_sequence_coverage_counts(
            self.detail_seed,
            self.consensus_aligner.coordinate_name,
            self.consensus_aligner.consensus,
            self.consensus_aligner.consensus_offset,
            self.consensus_aligner.seed_nucs)
        if not self.consensus_aligner.coordinate_name or not self.contigs:
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
                                                self.consensus_aligner.coordinate_name,
                                                contig_seq)

    def write_sequence_coverage_counts(self,
                                       contig_name: str,
                                       coordinate_name: typing.Optional[str],
                                       consensus: str,
                                       consensus_offset=0,
                                       seed_nucs=None):
        if not coordinate_name:
            coordinate_name = None
            alignments = []
        else:
            consensus_aligner = ConsensusAligner(self.projects)
            consensus_aligner.start_contig(coordinate_name, consensus)
            alignments = consensus_aligner.alignments
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
                action = CigarActions.MATCH if coordinate_name is None else None
                cigar = [(bead.end - bead.start, action)]
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
                        seed_nuc = seed_nucs[offset_seq_pos - 1]
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
                               'genome.pos',
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
                     seed_nuc: SeedNucleotide,
                     report_nuc: ReportNucleotide,
                     nuc_writer: DictWriter,
                     genome_start_pos=1):
        """ Write rows of nucleotide counts for a single codon.

        :param seed: the seed reference name
        :param region: the coordinate reference name
        :param seed_nuc: the SeedNucleotide object for this position
        :param report_nuc: the ReportNucleotide object for this position
        :param nuc_writer: the CSV writer to write the row into
        :param genome_start_pos: one-based start position of region
        """
        if seed_nuc.consensus_index is None:
            query_pos_txt = ''
        else:
            query_pos_txt = str(seed_nuc.consensus_index + 1)
        ref_pos = (str(report_nuc.position)
                   if report_nuc.position is not None
                   else '')
        genome_pos = (str(report_nuc.position+genome_start_pos - 1)
                      if report_nuc.position is not None
                      else '')
        row = {'seed': seed,
               'region': region,
               'q-cutoff': self.qcut,
               'query.nuc.pos': query_pos_txt,
               'refseq.nuc.pos': ref_pos,
               'genome.pos': genome_pos,
               'del': seed_nuc.counts['-'],
               'ins': seed_nuc.insertion_count,
               'clip': seed_nuc.clip_count,
               'v3_overlap': seed_nuc.v3_overlap,
               'coverage': seed_nuc.get_coverage()}
        for base in 'ACTGN':
            nuc_count = seed_nuc.counts[base]
            row[base] = nuc_count
        for field_name in ('coverage',
                           'clip',
                           'N',
                           'ins',
                           'del',
                           'v3_overlap'):
            if row[field_name]:
                nuc_writer.writerow(row)
                break

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
                report_nucleotides = self.report_nucleotides[region]
                report_nuc_index = 0
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
                        if len(report_nucleotides) <= report_nuc_index:
                            break
                        report_nuc = report_nucleotides[report_nuc_index]
                        query_pos = j + consensus_nuc_index + 1
                        seed_nuc.clip_count = seed_clipping[query_pos]
                        seed_nuc.insertion_count = seed_insertion_counts[query_pos]
                        if j == 2:
                            insertion_counts = self.insert_writer.insert_pos_counts[
                                (self.seed, region)]
                            seed_nuc.insertion_count += (
                                insertion_counts[report_amino.position])
                            report_nuc.seed_nucleotide.insertion_count += (
                                insertion_counts[report_amino.position])
                        if (report_nuc.seed_nucleotide.consensus_index in
                                (None, seed_nuc.consensus_index)):
                            report_nuc.seed_nucleotide.clip_count = \
                                seed_clipping[query_pos]
                            report_nuc.seed_nucleotide.insertion_count += \
                                seed_insertion_counts[query_pos]
                            report_nuc_index += 1

    def write_nuc_detail_counts(self, nuc_detail_writer=None):
        nuc_detail_writer = nuc_detail_writer or self.nuc_detail_writer
        self.write_nuc_report(nuc_detail_writer,
                              self.report_nucleotides,
                              self.detail_seed)

    def write_nuc_counts(self, nuc_writer=None):
        nuc_writer = nuc_writer or self.nuc_writer

        if self.combined_report_nucleotides:
            for seed, reports in self.combined_report_nucleotides.items():
                self.write_nuc_report(nuc_writer, reports, seed)
        elif self.reports:
            self.write_nuc_report(nuc_writer, self.report_nucleotides, self.seed)

    def write_nuc_report(self,
                         nuc_writer: DictWriter,
                         reports: typing.Dict[str, typing.List[ReportNucleotide]],
                         seed: str):
        self.merge_extra_counts()
        landmark_reader = LandmarkReader(self.landmarks)
        report_nucleotides: typing.List[ReportNucleotide]
        for region, report_nucleotides in sorted(reports.items()):
            try:
                gene = landmark_reader.get_gene(self.seed, region)
                genome_start_pos = gene['start']
            except ValueError:
                genome_start_pos = 1
            for report_nucleotide in report_nucleotides:
                self.write_counts(seed,
                                  region,
                                  report_nucleotide.seed_nucleotide,
                                  report_nucleotide,
                                  nuc_writer,
                                  genome_start_pos)

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

    def _write_consensus_helper(
            self,
            amino_entries,
            csv_writer,
            row_metadata,
            ignore_coverage=False,
            is_nucleotide=False,
            discard_deletions=True,
    ):
        for row in self.consensus_builder.get_consensus_rows(amino_entries,
                                                             ignore_coverage=ignore_coverage,
                                                             is_nucleotide=is_nucleotide,
                                                             discard_deletions=discard_deletions):
            row.update(row_metadata)
            if 'region-offset' not in row_metadata:
                row["offset"] = row.pop("seed-offset")
            csv_writer.writerow(row)

    def write_consensus_header(self, conseq_file):
        self.conseq_writer = self._create_consensus_writer(conseq_file)
        self.conseq_writer.writeheader()

    def write_consensus(self, conseq_writer=None):
        conseq_writer = conseq_writer or self.conseq_writer
        seed_amino_entries = [
            (seed_amino.consensus_nuc_index, seed_amino)
            for seed_amino in self.seed_aminos[0]
        ]
        self._write_consensus_helper(
            seed_amino_entries,
            conseq_writer,
            {"region": self.detail_seed},
        )

    def write_whole_genome_consensus_from_nuc(self, conseq_stitched_writer=None):
        conseq_stitched_writer = conseq_stitched_writer or self.conseq_stitched_writer
        landmark_reader = LandmarkReader(self.landmarks)
        for entry in self.combined_report_nucleotides:
            nuc_dict = {}
            insertions_dict = {}
            for seed_landmarks in self.landmarks:
                if re.fullmatch(seed_landmarks['seed_pattern'], entry):
                    coordinate_name = seed_landmarks['coordinates']
                    break
            else:
                coordinate_name = None
                if entry is not None and entry != '':
                    logger.warning(f'No coordinate reference found for entry: {entry}')
                continue
            sorted_regions: list = sorted(
                self.combined_report_nucleotides[entry].items(),
                key=lambda item: landmark_reader.get_gene(
                    coordinate_name,
                    item[0])['start'])
            regions_dict = dict(sorted_regions)
            for region in regions_dict:
                region_info = landmark_reader.get_gene(coordinate_name,
                                                       region)
                region_start = region_info['start']
                nuc_dict = combine_region_nucleotides(
                    nuc_dict,
                    self.combined_report_nucleotides[entry][region],
                    region_start)
                insertions_dict = combine_region_insertions(insertions_dict,
                                                            self.combined_insertions[entry][region],
                                                            region_start)
            nuc_dict = insert_insertions(insertions_dict, nuc_dict)
            nuc_entries = list(nuc_dict.items())
            nuc_entries.sort(key=lambda elem: elem[0])
            self._write_consensus_helper(
                nuc_entries,
                conseq_stitched_writer,
                {
                    "seed": entry,
                    "region": "whole genome consensus",
                },
                is_nucleotide=True,
            )

    def write_consensus_all_header(self, conseq_all_file):
        self.conseq_all_writer = self._create_consensus_writer(
            conseq_all_file,
            include_seed=True,
            include_seed_region_offsets=True,
        )
        self.conseq_all_writer.writeheader()

    def write_consensus_stitched_header(self, conseq_stitched_file):
        self.conseq_stitched_writer = self._create_consensus_writer(
            conseq_stitched_file,
            include_seed=True,
        )
        self.conseq_stitched_writer.writeheader()

    def write_consensus_all(self, conseq_all_writer=None):
        conseq_all_writer = conseq_all_writer or self.conseq_all_writer
        seed_amino_entries = [
            (seed_amino.consensus_nuc_index, seed_amino)
            for seed_amino in self.seed_aminos[0]
        ]
        self._write_consensus_helper(
            seed_amino_entries,
            conseq_all_writer,
            {
                "seed": self.detail_seed,
                "region": None,
                "region-offset": None,
            },
            ignore_coverage=True,
        )

        regions = sorted(self.reports.keys())
        for region in regions:
            consensus_nuc_indexes = set(chain(*(
                report_amino.seed_amino.all_consensus_nuc_indexes
                for report_amino in self.reports[region])))
            if not consensus_nuc_indexes:
                continue
            min_index = min(consensus_nuc_indexes)
            max_index = max(consensus_nuc_indexes) + 2
            region_amino_start = min(
                report_amino.position
                for report_amino in self.reports[region]
                if min_index in report_amino.seed_amino.all_consensus_nuc_indexes)
            region_nuc_offset = 3*(region_amino_start-1)
            if min_index < 0:
                region_nuc_offset -= min_index
            reading_frame = -min_index % 3
            region_amino_entries = [
                (seed_amino.consensus_nuc_index, seed_amino)
                for seed_amino in self.seed_aminos[reading_frame]
                if min_index <= seed_amino.consensus_nuc_index <= max_index
            ]
            self._write_consensus_helper(
                region_amino_entries,
                conseq_all_writer,
                {
                    "seed": self.seed,
                    "region": region,
                    "region-offset": region_nuc_offset
                },
                ignore_coverage=True,
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
                (3 * (report_amino.position-1), report_amino.seed_amino)
                for report_amino in self.reports[region]
            ]
            self._write_consensus_helper(
                region_amino_entries,
                conseq_region_writer,
                {
                    "seed": self.seed,
                    "region": region,
                },
                discard_deletions=False,
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
        insert_writer.write(self.inserts,
                            self.detail_seed,
                            self.reports,
                            self.report_nucleotides,
                            self.landmarks,
                            self.consensus_builder)

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


class InsertionWriter(object):
    def __init__(self, insert_file):
        """ Initialize a writer object.

        @param insert_file: an open file that the data will be written to
        """
        self.insert_file = insert_file
        self.insert_writer = csv.DictWriter(insert_file,
                                            ["seed",
                                             "mixture_cutoff",
                                             "region",
                                             "ref_region_pos",
                                             "ref_genome_pos",
                                             "query_pos",
                                             "insertion"],
                                            lineterminator=os.linesep
                                            )
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
        # insertions relative to consensus, by consensus position:
        self.conseq_insertions = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: SeedNucleotide())))
        # insertions relative to consensus, by ref position:
        self.ref_insertions = defaultdict(
            lambda: defaultdict(lambda: defaultdict(lambda: SeedNucleotide())))

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

    def write(self, insertions, seed_name, report_aminos_all, report_nucleotides_all, landmarks, consensus_builder):
        """ Write any insert ranges to the file.

        Sequence data comes from the reads that were added to the current group.
        @param insertions: regions and indexes of positions in the reads that should be
            reported as insertions.
        @param seed_name: name of the current seed
        @param report_aminos_all: a list of ReportAmino objects that represent the
            sequence that successfully mapped to the coordinate reference.
        @param report_nucleotides_all: a list of ReportNucleotides objects that
            represent the sequence, in case the region is not translated.
        @param landmarks: landmarks for the seed
        @param consensus_builder: helper function to write insertion consensus
        """
        if len(insertions) == 0:
            return

        for region, inserts in insertions.items():
            self.ref_insertions[region] = defaultdict(lambda: defaultdict(lambda: SeedNucleotide()))

            report_aminos = report_aminos_all[region] or []
            report_nucleotides = report_nucleotides_all[region] or []

            region_insert_pos_counts = self.insert_pos_counts[(self.seed, region)]
            inserts = list(inserts)
            inserts.sort()

            # convert insertion coordinates into contiguous ranges
            insert_ranges = []
            if len(report_aminos) == 0 and len(report_nucleotides) != 0:
                distance_insertions = 1
            else:
                distance_insertions = 3
            for insert in inserts:
                if not insert_ranges or insert != insert_ranges[-1][1]:
                    # just starting or we hit a gap
                    insert_ranges.append([insert, insert + distance_insertions])
                else:
                    insert_ranges[-1][1] += distance_insertions

            # enumerate insertions by popping out all AA sub-string variants
            insert_behind = {}
            insert_nuc_counts = {}
            insertion_coverage = {}  # max coverage left and right of the insertion
            for left, right in insert_ranges:
                current_insert_behind, current_insert_coverage =\
                    get_insertion_info(right, report_aminos, report_nucleotides)
                insert_behind[left] = current_insert_behind
                insertion_coverage[left] = current_insert_coverage
                current_nuc_counts = Counter()
                insert_nuc_counts[left] = current_nuc_counts
                for nuc_seq, count in self.nuc_seqs.items():
                    insert_nuc_seq = nuc_seq[left:right]
                    is_valid = (insert_nuc_seq and
                                'n' not in insert_nuc_seq and
                                '-' not in insert_nuc_seq)
                    if is_valid:
                        current_nuc_counts[insert_nuc_seq] += count

            for left in insert_nuc_counts.keys():
                ref_pos = insert_behind.get(left) - 1
                self.ref_insertions[region][ref_pos] =\
                    aggregate_insertions(insert_nuc_counts[left],
                                         consensus_pos=left-1,
                                         coverage_nuc=insertion_coverage[left])

            for left, counts in insert_nuc_counts.items():
                insert_pos = insert_behind.get(left) // 3
                count = sum(counts.values())
                # Only care about insertions in the middle of the sequence,
                # so ignore any that come before or after the reference.
                # Also report if we're in test mode (no report_aminos).
                if not report_aminos or insert_pos not in (0, None):
                    if insert_pos is not None:
                        region_insert_pos_counts[insert_pos] += count

        self.parse_conseq_insertions(seed_name, report_nucleotides_all)
        self.write_insertions_file(landmarks, consensus_builder)

    def parse_conseq_insertions(self, seed_name, report_nucleotides_all):
        for insertion_position in self.conseq_insertions[seed_name].keys():
            insertions = self.conseq_insertions[seed_name][insertion_position]
            for region in report_nucleotides_all:
                report_aminos = []
                report_nucleotides = report_nucleotides_all[region]
                current_insert_behind, current_insert_coverage = \
                    get_insertion_info(insertion_position, report_aminos, report_nucleotides)
                if current_insert_behind is not None:
                    for position in insertions:
                        insertions[position].count_nucleotides('-', count=current_insert_coverage)
                    if self.ref_insertions[region] is None:
                        self.ref_insertions[region] = defaultdict(lambda: defaultdict(lambda: SeedNucleotide()))
                    self.ref_insertions[region][current_insert_behind - 1] = insertions

    def write_insertions_file(self, landmarks, consensus_builder):
        landmark_reader = LandmarkReader(landmarks)
        for region in self.ref_insertions.keys():
            try:
                region_info = landmark_reader.get_gene(self.seed, region)
                region_start = region_info['start']
            except ValueError:
                region_start = 1
            if self.ref_insertions[region] is not None:
                for ref_pos in self.ref_insertions[region].keys():
                    insertions = list(self.ref_insertions[region][ref_pos].items())
                    for row in consensus_builder.get_consensus_rows(insertions, is_nucleotide=True):
                        insertions_row = dict(insertion=row["sequence"],
                                              seed=self.seed,
                                              mixture_cutoff=row['consensus-percent-cutoff'],
                                              region=region,
                                              ref_region_pos=ref_pos + 1,
                                              ref_genome_pos=ref_pos + region_start,
                                              query_pos=insertions[0][1].consensus_index + 1)
                        self.insert_writer.writerow(insertions_row)


def format_cutoff(cutoff):
    """ Format the cutoff fraction as a string to use as a name. """

    if cutoff == MAX_CUTOFF:
        return cutoff
    return '{:0.3f}'.format(cutoff)


def aln2counts(aligned_csv,
               nuc_csv,
               amino_csv,
               insertions_csv,
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
               conseq_stitched_csv=None,
               minimap_hits_csv=None):
    """
    Analyze aligned reads for nucleotide and amino acid frequencies.
    Generate consensus sequences.
    @param aligned_csv:         Open file handle containing aligned reads (from sam2aln)
    @param nuc_csv:             Open file handle to write nucleotide frequencies.
    @param amino_csv:           Open file handle to write amino acid frequencies.
    @param insertions_csv:      Open file handle to write insertions (new).
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
    @param conseq_stitched_csv: Open file handle to write stitched whole genome
        consensus sequences.
    @param minimap_hits_csv: Open file handle to write minimap2 match locations.
    """
    # load project information
    projects = ProjectConfig.loadDefault()

    landmarks_path = (Path(__file__).parent.parent / 'data' /
                      'landmark_references.yaml')
    landmarks_yaml = landmarks_path.read_text()

    # initialize reporter classes
    with InsertionWriter(insertions_csv) as insert_writer:
        report = SequenceReport(insert_writer,
                                projects,
                                CONSEQ_MIXTURE_CUTOFFS,
                                landmarks_yaml=landmarks_yaml,
                                consensus_min_coverage=CONSENSUS_MIN_COVERAGE)
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
        if conseq_stitched_csv is not None:
            report.write_consensus_stitched_header(conseq_stitched_csv)
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
            report.combined_report_nucleotides.clear()
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
    logging.basicConfig(level=logging.DEBUG)
    args = parse_args()
    aln2counts(args.aligned_csv,
               args.nuc_csv,
               args.amino_csv,
               args.insertions_csv,
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
               conseq_stitched_csv=args.conseq_stitched_csv,
               minimap_hits_csv=args.minimap_hits_csv)


if __name__ == '__main__':
    main()
