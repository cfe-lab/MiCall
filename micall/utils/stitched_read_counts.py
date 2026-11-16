"""Carry original remap counts to stitched sequences without positional aliasing.

The referenceless stitcher can reorder, merge, or omit its input contigs.
Original remap-count row numbers are therefore NOT stitched row numbers.
"""
import csv
from collections import Counter
from typing import Dict, TextIO

from Bio.Seq import Seq
from micall.core.remap import get_contig_read_count
from micall.utils.referencefull_contig_stitcher import read_remap_counts


def molecule_key(sequence: str) -> str:
    """Identify a complete molecule regardless of its reported strand."""
    dna = sequence.upper()
    reverse = str(Seq(dna).reverse_complement())
    return min(dna, reverse)


def project_stitched_counts(
    unstitched_csv: TextIO,
    stitched_csv: TextIO,
    original_counts_csv: TextIO,
) -> Dict[int, int]:
    """Transfer counts by stitched row only when identity is demonstrated.

    Counts are only exact when an original sequence and its stitched output
    each occur exactly once, including reverse-complement equivalence.
    Distinct original contigs can receive the *same* reads under multi-mapping,
    so summing duplicate contig counts would fabricate molecule support.
    Ambiguous, joined, altered or unmapped sequences retain UNKNOWN support.
    The downstream stitcher already handles a missing count.
    """
    source_rows = list(csv.DictReader(unstitched_csv))
    target_rows = list(csv.DictReader(stitched_csv))
    source_keys = [molecule_key(row["contig"]) for row in source_rows]
    target_keys = [molecule_key(row["contig"]) for row in target_rows]
    original_duplicates = Counter(source_keys)
    stitched_duplicates = Counter(target_keys)
    counts = read_remap_counts(original_counts_csv)
    supported = {}
    for index, (row, sequence) in enumerate(zip(source_rows, source_keys), 1):
        if (not sequence or original_duplicates[sequence] != 1
                or stitched_duplicates[sequence] != 1):
            continue
        count = get_contig_read_count(index, row["ref"],
                                      float(row["match"]), counts)
        if count is not None:
            supported[sequence] = count

    # The caller owns the stitched row identities. Do not serialize a
    # fake remap-count CSV or invent stitched contig names.
    projected = {}
    for index, sequence in enumerate(target_keys, 1):
        if sequence in supported:
            projected[index] = supported[sequence]
    return projected
