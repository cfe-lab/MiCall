#! /usr/bin/env python

"""
This module computes k-mer frequencies (and matches in contexts) from
given FASTA sequences.  It processes each contig in the FASTA file
and, for each k-mer (of sizes 1 up to a maximum), computes the number
of occurrences in the input sequence.
"""

import argparse
import sys
from dataclasses import dataclass
from collections import defaultdict
from pathlib import Path
from typing import Sequence, Iterator, TextIO, TypedDict, MutableMapping
import logging
from micall.utils.contig_stitcher_contigs import Contig
from Bio import SeqIO
from itertools import chain
import csv


logging.basicConfig(
    level=logging.INFO,
    format="%(levelname)s: %(message)s",
)
logger = logging.getLogger(__name__)


class Row(TypedDict, total=True):
    """
    A TypedDict for CSV rows.
    """

    size: int
    matches: int
    occurrences: int
    kmer: str


# FIELDNAMES defines the header for the CSV output.
FIELDNAMES = tuple(Row.__annotations__.keys())


class UserError(ValueError):
    """
    Custom error for user-facing messages.
    """

    def __init__(self, fmt: str, *fmt_args: object):
        self.fmt = fmt
        self.fmt_args = fmt_args
        self.code = 1


@dataclass(frozen=True)
class KMer:
    """
    Represents a k-mer from a sequence along with its flanking context.

    Attributes:
        sequence: The k-mer substring.
        left: All bases to the left of the k-mer within the sequence.
        right: All bases to the right of the k-mer within the sequence.
    """

    sequence: str
    left: str
    right: str


@dataclass(frozen=True)
class KMerWithCounter:
    """
    Combines a KMer with its associated counter that maps matching
    score to occurrence.

    Attributes:
        kmer: The KMer instance.
        counter: A mapping from match count (int) to number of
                 occurrences (int) of context k-mers that have this
                 many matching positions.
    """

    kmer: KMer
    counter: MutableMapping[int, int]


# Pool is defined as a mapping from k-mer string to its
# KMerWithCounter.
Pool = MutableMapping[str, KMerWithCounter]


def read_contigs(input_fasta: TextIO) -> Iterator[Contig]:
    """
    Read contig records from a FASTA file.

    Parameters:
        input_fasta (TextIO): A file-like object in FASTA format.

    Yields:
        Contig: A contig with a name and a sequence string.
    """

    for record in SeqIO.parse(input_fasta, "fasta"):
        yield Contig(name=record.name, seq=str(record.seq), reads_count=None)


def kmers(sequence: str, size: int) -> Iterator[KMer]:
    """
    Generate all k-mers of a specified size from a DNA sequence along
    with their left and right contexts.

    Explanation:
      For each valid window in the sequence (using a sliding window approach),
      yields a KMer object with:
         - sequence: the current k-mer substring,
         - left: substring preceding the k-mer,
         - right: substring following the k-mer.

    Parameters:
        sequence (str): The DNA sequence.
        size (int): The desired k-mer size (window length).

    Yields:
        KMer: A k-mer with its associated context.
    """

    for i in range(len(sequence) - size + 1):
        end = i+size
        seq = sequence[i:end]
        left = sequence[:i]
        right = sequence[end:]
        yield KMer(sequence=seq, left=left, right=right)


def slide_kmer(with_counter: KMerWithCounter) -> None:
    """
    Compute match counts between a given k-mer and all k-mers
    generated from its flanking contexts.

    For the provided KMerWithCounter object, this function slides a
    window over the left and right context sequences to generate
    k-mers of the same size. For each such context k-mer, it compares
    the bases with the original k-mer, tallying the number of matching
    positions. Results are stored in the provided counter mapping
    within the KMerWithCounter.

    Parameters:
        with_counter (KMerWithCounter): Combines a KMer with its
        counter where counts will be updated.
    """

    kmer = with_counter.kmer
    size = len(kmer.sequence)
    counter = with_counter.counter

    for other in chain(kmers(kmer.left, size),
                       kmers(kmer.right, size),
                       ):
        matches = sum(x == y for x, y in zip(kmer.sequence, other.sequence))
        if matches > 0:
            counter[matches] += 1


def process_contig(pool: Pool,
                   contig: Contig,
                   max_kmer: int,
                   ) -> None:
    """
    Process a single contig to compute k-mer statistics for k = 1 to
    max_kmer.

    For each k value:
      1. All k-mers are generated from the contig sequence.
      2. Duplicates (based on k-mer sequence) are removed.
      3. Each unique k-mer is then tested against its flanking context
      using slide_kmer().
      4. Results are accumulated within a shared pool (deduplicating
      k-mers across contigs).

    Parameters:
        pool (Pool): A dictionary that accumulates KMerWithCounter
        objects across contigs.
        contig (Contig): The input contig to be processed.
        max_kmer (int): Maximum k-mer length to analyze.
    """

    for size in range(1, max_kmer + 1):
        all_kmers = kmers(contig.seq, size)
        deduplicated = {x.sequence: x for x in all_kmers}.values()
        for kmer in deduplicated:
            if kmer.sequence in pool:
                existing = pool[kmer.sequence]
                with_counter = KMerWithCounter(kmer, existing.counter)
            else:
                with_counter = KMerWithCounter(kmer, defaultdict(int))
                pool[kmer.sequence] = with_counter
            slide_kmer(with_counter)


def process_all_contigs(pool: Pool,
                        contigs: Sequence[Contig],
                        max_kmer: int,
                        ) -> None:
    for i, contig in enumerate(contigs):
        logger.debug("Processing contig %s (%s/%s).",
                     contig.name, i+1, len(contigs))
        process_contig(pool, contig, max_kmer)


def counters_to_rows(pool: Pool) -> Iterator[Row]:
    for with_counter in pool.values():
        kmer = with_counter.kmer
        counter = with_counter.counter
        size = len(kmer.sequence)
        for matches, occurrences in counter.items():
            row: Row = {"size": size,
                        "matches": matches,
                        "occurrences": occurrences,
                        "kmer": kmer.sequence,
                        }
            yield row


def main_typed(input: Path, output: Path, max_kmer: int) -> None:
    """
    Main processing function: reads input FASTA, processes contigs,
    and writes CSV output.

    Parameters:
        input_path (Path): Path to the FASTA input file.
        output_path (Path): Path to the CSV output file.
        max_kmer (int): Maximum k-mer length to process.
    """

    if not input.exists():
        raise UserError("Input file %r does not exist.", str(input))

    with input.open() as input_file:
        contigs = tuple(read_contigs(input_file))
        logger.debug("Read %s input sequences.", len(contigs))

    pool: Pool = {}
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w") as output_file:
        writer = csv.DictWriter(output_file, fieldnames=FIELDNAMES)
        writer.writeheader()
        process_all_contigs(pool, contigs, max_kmer)
        for row in counters_to_rows(pool):
            writer.writerow(row)


def parse_arguments(argv: Sequence[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Calculate k-mer frequencies from FASTA sequences and output results as CSV.")

    parser.add_argument('input', type=Path,
                        help='Input FASTA file.')
    parser.add_argument('output', type=Path,
                        help='Input CSV file.')
    parser.add_argument('--max', type=int, default=10,
                        help='Largest k-mer to be analyzed.')

    verbosity_group = parser.add_mutually_exclusive_group()
    verbosity_group.add_argument('--verbose', action='store_true',
                                 help='Increase output verbosity.')
    verbosity_group.add_argument('--no-verbose', action='store_true',
                                 help='Normal output verbosity.', default=True)
    verbosity_group.add_argument('--debug', action='store_true',
                                 help='Maximum output verbosity.')
    verbosity_group.add_argument('--quiet', action='store_true',
                                 help='Minimize output verbosity.')

    return parser.parse_args(argv)


def configure_logging(args: argparse.Namespace) -> None:
    if args.quiet:
        logger.setLevel(logging.ERROR)
    elif args.debug:
        logger.setLevel(logging.DEBUG)
    elif args.verbose:
        logger.setLevel(logging.INFO)
    else:
        logger.setLevel(logging.WARNING)


def main(argv: Sequence[str]) -> int:
    args = parse_arguments(argv)
    configure_logging(args)

    try:
        main_typed(args.input, args.output, args.max)
        logger.debug("Done.")
        return 0
    except BrokenPipeError:
        logger.debug("Broken pipe.")
        return 1
    except KeyboardInterrupt:
        logger.debug("Interrupted.")
        return 1
    except UserError as e:
        logger.fatal(e.fmt, *e.fmt_args)
        return e.code


def entry() -> None:
    sys.exit(main(sys.argv[1:]))


if __name__ == "__main__": entry()  # noqa
