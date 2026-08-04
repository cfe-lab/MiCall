from typing import Iterable, Optional, Tuple, Dict, Literal, TypeVar, TextIO, Sequence, List
from collections import defaultdict
import csv
import os
from pathlib import Path
from dataclasses import replace
from math import ceil
from functools import reduce
from itertools import tee, islice, chain
from queue import LifoQueue
from Bio import Seq
import logging
from aligntools import CigarHit, connect_nonoverlapping_cigar_hits, drop_overlapping_cigar_hits, CigarActions

from micall.core.project_config import ProjectConfig
from micall.core.plot_contigs import plot_stitcher_coverage
from micall.utils.contig_stitcher_context import ReferencefullStitcherContext
from micall.utils.contig_stitcher_contigs import GenotypedContig, AlignedContig
from micall.utils.consensus_aligner import align_consensus
from micall.utils.overlap_stitcher import align_queries, calculate_concordance_norm, sort_concordance_indexes
import micall.utils.referencefull_contig_stitcher_events as events


T = TypeVar("T")
logger = logging.getLogger(__name__)


def log(e: events.EventType) -> None:
    ReferencefullStitcherContext.get().emit(e)
    if isinstance(e, events.Warning):
        logger.warning("%s", e)
    else:
        logger.debug("%s", e)


def cut_query(self: GenotypedContig, cut_point: float) -> Tuple[GenotypedContig, GenotypedContig]:
    """ Cuts query sequence in two parts with cut_point between them. """

    cut_point = max(0.0, cut_point)
    left_len = ceil(cut_point)
    total_len = len(self.seq)

    # Distribute reads_count proportionally based on sequence length
    left_reads_count: Optional[int] = None
    right_reads_count: Optional[int] = None
    if self.reads_count is not None and total_len > 0:
        left_reads_count = round(self.reads_count * left_len / total_len)
        right_reads_count = self.reads_count - left_reads_count

    left = replace(self, name=None, seq=self.seq[:left_len], reads_count=left_reads_count)
    right = replace(self, name=None, seq=self.seq[left_len:], reads_count=right_reads_count)
    return left, right


def cut_reference(self: AlignedContig, cut_point: float) -> Tuple[AlignedContig, AlignedContig]:
    """ Cuts this alignment in two parts with cut_point between them. """

    alignment_left, alignment_right = self.alignment.cut_reference(cut_point)
    left = replace(self, name=None, alignment=alignment_left)
    right = replace(self, name=None, alignment=alignment_right)
    log(events.Cut(self, left, right, cut_point))
    return left, right


def lstrip(self: AlignedContig) -> AlignedContig:
    """
    Trims the query sequence of the contig from its beginning up to the start of the
    alignment. The CIGAR alignment is also updated to reflect the trimming.
    """

    alignment = self.alignment.lstrip_reference().lstrip_query()
    q_remainder, query = cut_query(self, alignment.q_st - 0.5)
    alignment = alignment.translate(0, -1 * alignment.q_st)
    result = AlignedContig.make(query, alignment, self.strand)
    log(events.LStrip(self, result))
    return result


def rstrip(self: AlignedContig) -> AlignedContig:
    """
    Trims the query sequence of the contig from its end based on the end of the
    alignment. The CIGAR alignment is also updated to reflect the trimming.
    """

    alignment = self.alignment.rstrip_reference().rstrip_query()
    query, q_remainder = cut_query(self, alignment.q_ei + 0.5)
    result = AlignedContig.make(query, alignment, self.strand)
    log(events.RStrip(self, result))
    return result


def overlap(a: AlignedContig, b: AlignedContig) -> bool:
    def intervals_overlap(x, y):
        return x[0] <= y[1] and x[1] >= y[0]

    if a.group_ref != b.group_ref:
        return False

    return intervals_overlap((a.alignment.r_st, a.alignment.r_ei),
                             (b.alignment.r_st, b.alignment.r_ei))


def munge(self: AlignedContig, other: AlignedContig) -> AlignedContig:
    """
    Combines two adjacent contigs into a single contig by joining their
    query sequences and alignments.
    """

    match_fraction = min(self.match_fraction, other.match_fraction)
    ref_name = max([self, other], key=lambda x: x.alignment.ref_length).ref_name
    query = GenotypedContig(seq=self.seq + other.seq,
                            name=None,
                            ref_name=ref_name,
                            group_ref=self.group_ref,
                            ref_seq=self.ref_seq,
                            match_fraction=match_fraction,
                            reads_count=None)  # Combined contigs lose their reads count

    self_alignment = self.alignment
    other_alignment = \
        other.alignment.translate(
            query_delta=(-1 * other.alignment.q_st + self.alignment.q_ei + 1),
            reference_delta=0)
    alignment = self_alignment.connect(other_alignment)

    ret = AlignedContig.make(query=query, alignment=alignment, strand=self.strand)
    log(events.Munge(self, other, ret))
    return ret


def sliding_window(sequence: Iterable[T]) -> Iterable[Tuple[Optional[T], T, Optional[T]]]:
    """
    Generate a three-element sliding window of a sequence.

    Each element generated contains a tuple with the previous item (None if the first item),
    the current item, and the next item (None if the last item) in the sequence.
    """

    a, b, c = tee(sequence, 3)
    prevs = chain([None], a)
    nexts = chain(islice(c, 1, None), [None])
    return zip(prevs, b, nexts)


def combine_contigs(parts: Sequence[AlignedContig]) -> AlignedContig:
    """
    Combine a list of contigs into a single AlignedContig by trimming and merging overlapping parts.

    Left-trimming and right-trimming occur at any shared overlapping points
    between adjacent parts. munge() is used to combine contiguous parts without overlap.
    """

    stripped_parts = []
    for prev_part, part, next_part in sliding_window(parts):
        if prev_part is not None:
            part = lstrip(part)
        if next_part is not None:
            part = rstrip(part)
        stripped_parts.append(part)

    ret = reduce(munge, stripped_parts)
    log(events.Combine(tuple(stripped_parts), ret))
    return ret


def align_to_reference(contig: GenotypedContig) -> Iterable[GenotypedContig]:
    """
    Align a single Contig to its reference sequence, producing potentially multiple aligned contigs.

    If the reference sequence (ref_seq) is unavailable, the contig is returned unaltered.
    Otherwise, alignments are performed and contigs corresponding to each alignment are yielded.
    """

    if contig.ref_seq is None:
        log(events.NoRef(contig))
        yield contig
        return

    alignments, _algo = align_consensus(contig.ref_seq, contig.seq)
    hits = [x.to_cigar_hit() for x in alignments]
    strands: Tuple[Literal["forward", "reverse"], ...] = tuple(
        "forward" if x.strand == 1 else "reverse" for x in alignments)

    for i, (hit, strand) in enumerate(zip(hits, strands)):
        log(events.InitialHit(contig, i, hit, strand))

    if not hits:
        log(events.ZeroHits(contig))
        yield contig
        return

    if len(set(strands)) > 1:
        log(events.StrandConflict(contig, hits, strands))
        yield contig
        return

    strand = strands[0]
    if strand == "reverse":
        rc = str(Seq.Seq(contig.seq).reverse_complement())
        original_contig = contig
        new_contig = replace(contig, seq=rc)
        contig = new_contig
        hits = [replace(hit, q_st=len(rc)-hit.q_ei-1, q_ei=len(rc)-hit.q_st-1) for hit in hits]

        log(events.ReverseComplement(original_contig, new_contig))
        for i, (hit, strand) in enumerate(zip(hits, strands)):
            log(events.InitialHit(contig, i, hit, strand))

    def quality(x: CigarHit):
        mlen = sum(1 for x in x.cigar.relax().iterate_operations()
                   if x == CigarActions.MATCH)
        return (mlen, x.ref_length)

    filtered = tuple(drop_overlapping_cigar_hits(hits, quality))
    connected = tuple(connect_nonoverlapping_cigar_hits(filtered))
    log(events.HitNumber(contig, tuple(zip(hits, strands)), connected))

    for i, single_hit in enumerate(connected):
        query = replace(contig, name=None)
        part = AlignedContig.make(query, single_hit, strand)
        log(events.ConnectedHit(contig, part, i))
        yield part


def strip_conflicting_mappings(contigs: Iterable[GenotypedContig]) -> Iterable[GenotypedContig]:
    contigs = list(contigs)
    names = {contig.id: contig for contig in contigs}

    def get_indexes(id: int) -> Tuple[int, int]:
        contig = names[id]
        if isinstance(contig, AlignedContig):
            return contig.alignment.q_st, contig.alignment.r_st
        else:
            return -1, -1

    reference_sorted = list(sorted(names.keys(), key=lambda id: get_indexes(id)[1]))
    query_sorted = list(sorted(names.keys(), key=lambda id: get_indexes(id)[0]))

    def is_out_of_order(id: int) -> bool:
        return reference_sorted.index(id) != query_sorted.index(id)

    sorted_by_query = sorted(contigs, key=lambda contig: get_indexes(contig.id))
    for prev_contig, contig, next_contig in sliding_window(sorted_by_query):
        if isinstance(contig, AlignedContig):
            original = contig
            start = prev_contig.alignment.q_ei + 1 if isinstance(prev_contig, AlignedContig) else 0
            end = next_contig.alignment.q_st - 1 if isinstance(next_contig, AlignedContig) else len(contig.seq) - 1

            if prev_contig is not None or is_out_of_order(original.id):
                contig = lstrip(contig)
                log(events.InitialStrip(original, start, original.alignment.q_st - 1))
            if next_contig is not None or is_out_of_order(original.id):
                contig = rstrip(contig)
                log(events.InitialStrip(original, original.alignment.q_ei + 1, end))

        yield contig


def align_all_to_reference(contigs: Iterable[GenotypedContig]) -> Iterable[GenotypedContig]:
    """
    Align multiple contigs to their respective reference sequences.

    Applies align_to_reference to each contig in the given collection,
    flattening the result into a single list.
    """

    groups = map(align_to_reference, contigs)
    groups = map(strip_conflicting_mappings, groups)
    for group in groups:
        yield from group


def find_all_overlapping_contigs(self: AlignedContig, aligned_contigs):
    """
    Yield all contigs from a collection that overlap with a given contig.
    Contigs are considered overlapping if they have overlapping intervals on the same reference genome.
    """

    for other in aligned_contigs:
        if overlap(self, other):
            yield other


def find_overlapping_contig(self: AlignedContig, aligned_contigs):
    """
    Find the single contig in a collection that overlaps the most with a given contig.
    It returns the contig with the maximum overlapped reference length with the given contig (self).
    """

    every = find_all_overlapping_contigs(self, aligned_contigs)
    return max(every, key=lambda other: other.alignment.ref_length if other else 0, default=None)


def concordance_to_cut_points(left_overlap, right_overlap, aligned_left, aligned_right, concordance):
    """ Determine optimal cut points for stitching based on sequence concordance in the overlap region. """

    sorted_concordance_indexes = sort_concordance_indexes(concordance)

    def remove_dashes(s: str):
        return s.replace('-', '')

    for max_concordance_index in sorted_concordance_indexes:
        aligned_left_q_index = len(remove_dashes(aligned_left[:max_concordance_index]))
        aligned_right_q_index = right_overlap.alignment.query_length - \
            len(remove_dashes(aligned_right[max_concordance_index:])) + 1
        aligned_left_r_index = left_overlap.alignment.coordinate_mapping.query_to_ref.left_max(aligned_left_q_index)
        if aligned_left_r_index is None:
            aligned_left_r_index = left_overlap.alignment.r_st - 1
        aligned_right_r_index = right_overlap.alignment.coordinate_mapping.query_to_ref.right_min(aligned_right_q_index)
        if aligned_right_r_index is None:
            aligned_right_r_index = right_overlap.alignment.r_ei + 1
        if aligned_right_r_index > aligned_left_r_index:
            return aligned_left_r_index + 0.5, aligned_right_r_index - 0.5, max_concordance_index

    return left_overlap.alignment.r_st - 1 + 0.5, right_overlap.alignment.r_ei + 1 - 0.5, 0


def stitch_2_contigs(left, right):
    """
    Stitch two contigs together into a single coherent contig.

    The function handles the overlap by cutting both contigs into segments, aligning the
    overlapping segments, and then choosing the optimal stitching points based on sequence
    concordance. Non-overlapping segments are retained as is.
    """

    # Cut in 4 parts.
    left_remainder, left_overlap = cut_reference(left, right.alignment.r_st - 0.5)
    right_overlap, right_remainder = cut_reference(right, left.alignment.r_ei + 0.5)
    left_overlap = lstrip(rstrip(left_overlap))
    right_overlap = lstrip(rstrip(right_overlap))
    left_remainder = rstrip(left_remainder)
    right_remainder = lstrip(right_remainder)
    log(events.StitchCut(left, right, left_overlap, right_overlap, left_remainder, right_remainder))

    # Align overlapping parts, then recombine based on concordance.
    aligned_left, aligned_right = align_queries(left_overlap.seq, right_overlap.seq)
    concordance = calculate_concordance_norm(aligned_left, aligned_right)
    aligned_left_cutpoint, aligned_right_cutpoint, max_concordance_index = \
        concordance_to_cut_points(left_overlap, right_overlap, aligned_left, aligned_right, concordance)
    left_overlap_take, left_overlap_drop = cut_reference(left_overlap, aligned_left_cutpoint)
    right_overlap_drop, right_overlap_take = cut_reference(right_overlap, aligned_right_cutpoint)

    # Log it.
    cut_point_location_scaled = max_concordance_index / (((len(concordance) or 1) - 1) or 1)
    log(events.Overlap(left, right, left_overlap, right_overlap,
                       left_remainder, right_remainder, left_overlap_take,
                       right_overlap_take, tuple(concordance),
                       max_concordance_index, cut_point_location_scaled))

    return combine_contigs([left_remainder, left_overlap_take, right_overlap_take, right_remainder])


def combine_overlaps(contigs: Sequence[AlignedContig]) -> Iterable[AlignedContig]:
    """
    Repeatedly combine all overlapping aligned contigs into an iterable collection of contiguous AlignedContigs.
    It proceeds by iterating through sorted contigs and stitching any overlapping ones until none are left.
    """

    # Going left-to-right through aligned contigs.
    contigs = list(sorted(contigs, key=lambda x: x.alignment.r_st))
    while contigs:
        current = contigs.pop(0)

        # Find overlap. If there isn't one - we are done with the current contig.
        overlapping_contig = find_overlapping_contig(current, contigs)
        if not overlapping_contig:
            log(events.NoOverlap(current))
            yield current
            continue

        # Replace two contigs by their stitched version, then loop with it.
        new_contig = stitch_2_contigs(current, overlapping_contig)
        contigs.remove(overlapping_contig)
        contigs.insert(0, new_contig)
        log(events.Stitch(current, overlapping_contig, new_contig))


def merge_intervals(intervals: Sequence[Tuple[int, int]]) -> Sequence[Tuple[int, int]]:
    """
    Merge overlapping and adjacent intervals.
    Note that intervals are inclusive.

    :param intervals: A list of intervals [start, end] where 'start' and 'end' are integers.
    :return: A list of merged intervals.
    """

    if not intervals:
        return []

    # Sort intervals by their starting values
    sorted_intervals = sorted(intervals, key=lambda x: x[0])

    merged_intervals = [sorted_intervals[0]]
    for current in sorted_intervals[1:]:
        current_start, current_end = current
        last_start, last_end = merged_intervals[-1]
        if current_start <= last_end + 1:
            # Extend the last interval if there is an overlap or if they are adjacent
            merged_intervals[-1] = (min(last_start, current_start), max(last_end, current_end))
        else:
            # Add the current interval if there is no overlap
            merged_intervals.append(current)

    return merged_intervals


def find_covered_contig(contigs: Sequence[AlignedContig]) -> Tuple[Optional[AlignedContig], Sequence[AlignedContig]]:
    """
    Find and return the first contig that is completely covered by other contigs.

    :param contigs: Sequence of all aligned contigs to be considered.
    :return: An AlignedContig if there is one completely covered by others, None otherwise.
    """

    def calculate_cumulative_coverage(others) -> Sequence[Tuple[int, int]]:
        intervals = [(contig.alignment.r_st, contig.alignment.r_ei) for contig in others]
        merged_intervals = merge_intervals(intervals)
        return merged_intervals

    by_reads = all(contig.reads_count is not None for contig in contigs)
    def order(contig: AlignedContig) -> int:
        if by_reads:
            assert contig.reads_count is not None
            return contig.reads_count
        else:
            return contig.alignment.ref_length

    # Iterating from the least significant contig to the most significant one.
    for current in sorted(contigs, key=order):
        current_interval = (current.alignment.r_st, current.alignment.r_ei)

        # Create a map of cumulative coverage for contigs
        overlaping_contigs = [x for x in contigs if x.id != current.id and overlap(current, x)]
        cumulative_coverage = calculate_cumulative_coverage(overlaping_contigs)

        # Check if the current contig is covered by the cumulative coverage intervals
        if any((cover_interval[0] < current_interval[0] and cover_interval[1] >= current_interval[1]
                or cover_interval[0] <= current_interval[0] and cover_interval[1] > current_interval[1])
               for cover_interval in cumulative_coverage):
            # Strictly covered.
            return current, overlaping_contigs
        elif any((cover_interval[0] == current_interval[0] and cover_interval[1] == current_interval[1])
               for cover_interval in cumulative_coverage):
            if any(contig.reads_count is None for contig in overlaping_contigs + [current]):
                log(events.IgnoreCoverage(current, overlaping_contigs))
                continue

            # Unstrict coverage (exact boundary match).
            # In this case we must compare the read counts to decide which contig to keep.
            total_coverage = sum(contig.reads_count or 0 for contig in overlaping_contigs)
            if total_coverage > (current.reads_count or 0):
                return current, overlaping_contigs

    return None, []


def drop_completely_covered(contigs: Sequence[AlignedContig]) -> Sequence[AlignedContig]:
    """ Filter out all contigs that are contained within other contigs. """

    contigs = list(contigs)
    while contigs:
        covered, covering = find_covered_contig(contigs)
        if covered:
            contigs.remove(covered)
            log(events.Drop(covered, tuple(covering)))
        else:
            break

    return contigs


def split_contigs_with_gaps(contigs: Sequence[AlignedContig]) -> Sequence[AlignedContig]:
    """
    Split contigs at large gaps if those gaps are covered by other contigs in the list.

    A gap within a contig is considered large based on a pre-defined threshold. If another contig aligns
    within that gap's range, the contig is split into two around the midpoint of the gap.
    """

    contigs = list(contigs)

    def covered_by(gap, other):
        # Check if any 1 reference coordinate in gap is mapped in `other`.
        gap_coords = gap.coordinate_mapping.ref_to_query.domain
        cover_coords = set(other.alignment.coordinate_mapping.ref_to_query.keys())
        return not gap_coords.isdisjoint(cover_coords)

    def covered(self, gap):
        return any(covered_by(gap, other) for other in contigs if other != self)

    def significant(gap):
        # noinspection PyLongLine
        # The size of the gap is unavoidably, to some point, arbitrary. Here we tried to adjust it to common gaps in HIV, as HIV is the primary test subject in MiCall. A notable feature of HIV-1 reverse transcription is the appearance of periodic deletions of approximately 21 nucleotides. These deletions have been reported to occur in the HIV-1 genome and are thought to be influenced by the structure of the viral RNA. Specifically, the secondary structures and foldings of the RNA can lead to pause sites for the reverse transcriptase, resulting in staggered alignment when the enzyme slips. This misalignment can cause the reverse transcriptase to "jump," leading to deletions in the newly synthesized DNA. The unusually high frequency of about 21-nucleotide deletions is believed to correspond to the pitch of the RNA helix, which reflects the spatial arrangement of the RNA strands. The 21 nucleotide cycle is an average measure and is thought to be associated with the length of one turn of the RNA helix, meaning that when reverse transcriptase slips and reattaches, it often does so one helical turn away from the original site.         # noqa: E501
        return gap.ref_length > 21

    def try_split(self: AlignedContig):
        for gap in self.alignment.deletions():
            if not significant(gap):
                # Really we do not want to split on every little deletion
                # because that would mean that we would need to stitch
                # overlaps around them.
                # And we are likely to lose quality with every stitching operation.
                # By skipping we assert that this gap is aligner's fault.
                log(events.IgnoreGap(self, gap))
                continue

            if covered(self, gap):
                midpoint = gap.r_st + (gap.r_ei - gap.r_st) / 2 + self.alignment.epsilon
                left_part, right_part = cut_reference(self, midpoint)
                left_part = rstrip(left_part)
                right_part = lstrip(right_part)

                contigs.remove(self)
                contigs.append(left_part)
                contigs.append(right_part)
                process_queue.put(right_part)
                log(events.SplitGap(self, gap, left_part, right_part))
                return

    process_queue: LifoQueue = LifoQueue()
    for contig in contigs:
        process_queue.put(contig)

    while not process_queue.empty():
        contig = process_queue.get()
        try_split(contig)

    return contigs


def stitch_contigs(contigs: Iterable[GenotypedContig]) -> Iterable[GenotypedContig]:
    contigs = list(contigs)
    for contig in contigs:
        log(events.Intro(contig))
        contig.register()

    maybe_aligned = list(align_all_to_reference(contigs))

    # Contigs that did not align do not need any more processing
    yield from (x for x in maybe_aligned if not isinstance(x, AlignedContig))
    aligned: Sequence[AlignedContig] = [
        x for x in maybe_aligned if isinstance(x, AlignedContig)]

    aligned = split_contigs_with_gaps(aligned)
    aligned = drop_completely_covered(aligned)
    yield from combine_overlaps(aligned)


GroupRef = Optional[str]


def stitch_consensus(contigs: Iterable[GenotypedContig]) -> Iterable[GenotypedContig]:
    contigs = list(stitch_contigs(contigs))
    consensus_parts: Dict[GroupRef, List[AlignedContig]] = defaultdict(list)

    for contig in contigs:
        if isinstance(contig, AlignedContig):
            consensus_parts[contig.group_ref].append(contig)
        else:
            yield contig

    def combine(group_ref):
        ctgs = sorted(consensus_parts[group_ref], key=lambda x: x.alignment.r_st)
        result = combine_contigs(ctgs)
        log(events.FinalCombine(tuple(ctgs), result))
        return result

    yield from map(combine, consensus_parts)


def write_contigs(output_csv: TextIO, contigs: Iterable[GenotypedContig]):
    writer = csv.DictWriter(output_csv,
                            ['ref', 'match', 'group_ref', 'contig'],
                            lineterminator=os.linesep)
    writer.writeheader()
    for contig in contigs:
        writer.writerow(dict(ref=contig.ref_name,
                             match=contig.match_fraction,
                             group_ref=contig.group_ref,
                             contig=contig.seq))

    output_csv.flush()


def read_remap_counts(remap_counts_csv: TextIO) -> Dict[str, int]:
    """Read remap counts CSV and extract read counts per contig.

    Args:
        remap_counts_csv: Open file handle to remap_counts.csv

    Returns:
        Dictionary mapping contig names (like "1-HCV-1a") to read counts

    Raises:
        ValueError: If count values are not valid integers

    Note:
        - If multiple entries exist for the same contig (e.g., from different
          remap iterations), the last occurrence is used
        - Entries not starting with "remap" are ignored (e.g., "raw", "unmapped")
        - Malformed entries (missing space after remap prefix) are ignored
    """
    counts = {}
    reader = csv.DictReader(remap_counts_csv)

    for row_num, row in enumerate(reader, start=2):  # Start at 2 (header is row 1)
        type_field = row.get('type', '')

        # Extract contig name from type field like "remap 1-HCV-1a" or "remap-1 1-HCV-1a"
        if type_field.startswith('remap'):
            # Remove "remap" or "remap-1", "remap-2", etc. prefix
            parts = type_field.split(' ', 1)
            if len(parts) == 2:
                contig_name = parts[1]
                count_str = row.get('count', '0')

                try:
                    count = int(count_str)
                except ValueError as e:
                    raise ValueError(
                        f"Invalid count value '{count_str}' for contig '{contig_name}' "
                        f"at row {row_num} in remap_counts.csv"
                    ) from e

                counts[contig_name] = count

    return counts


def read_contigs(input_csv: TextIO, contig_read_counts: Dict[str, int]) -> Iterable[GenotypedContig]:
    projects = ProjectConfig.loadDefault()

    for i, row in enumerate(csv.DictReader(input_csv)):
        seq = row['contig']
        ref_name = row['ref']
        group_ref = row['group_ref']
        match_fraction = float(row['match'])

        try:
            ref_seq = projects.getGenotypeReference(group_ref)
        except KeyError:
            try:
                ref_seq = projects.getReference(group_ref)
            except KeyError:
                ref_seq = None

        # Generate contig name to look up read count (format: "1-HCV-1a")
        contig_name = f"{i+1}-{ref_name}"
        reads_count = contig_read_counts.get(contig_name)

        yield GenotypedContig(name=None,
                              seq=seq,
                              ref_name=ref_name,
                              group_ref=group_ref,
                              ref_seq=str(ref_seq) if ref_seq is not None else None,
                              match_fraction=match_fraction,
                              reads_count=reads_count)


def referencefull_contig_stitcher(input_csv: TextIO,
                                  output_csv: Optional[TextIO],
                                  stitcher_plot_path: Optional[Path],
                                  remap_counts_csv: Optional[TextIO] = None,
                                  ) -> int:
    with ReferencefullStitcherContext.fresh() as ctx:
        # Read remap counts if provided
        if remap_counts_csv is not None:
            contig_read_counts = read_remap_counts(remap_counts_csv)
        else:
            contig_read_counts = {}

        # Read contigs with read counts
        contigs = list(read_contigs(input_csv, contig_read_counts))

        if output_csv is not None or stitcher_plot_path is not None:
            contigs = list(stitch_consensus(contigs))

        if output_csv is not None:
            write_contigs(output_csv, contigs)

        if stitcher_plot_path is not None:
            plot_stitcher_coverage(ctx.events, stitcher_plot_path)

        return len(contigs)
