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

    qseqid: str
    size: int
    matches: int
    occurrences: int
    kmer: str


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
    Represents a k-mer and its context within a sequence.
    """

    sequence: str
    left: str
    right: str


def read_contigs(input_fasta: TextIO) -> Iterator[Contig]:
    """
    Read contigs from a FASTA file.

    Parameters:
        input_fasta (TextIO): Input file handle in FASTA format.

    Yields:
        Contig: A contig with a name and sequence.
    """

    for record in SeqIO.parse(input_fasta, "fasta"):
        yield Contig(name=record.name, seq=str(record.seq))


def kmers(sequence: str, size: int) -> Iterator[KMer]:
    """
    Generate all k-mers of a given size from a sequence along with
    their left and right context.

    Parameters:
        sequence (str): The input DNA sequence.
        size (int): The k-mer size.

    Yields:
        KMer: The k-mer and its context.
    """

    for i in range(len(sequence) - size + 1):
        end = i+size
        seq = sequence[i:end]
        left = sequence[:i]
        right = sequence[end:]
        yield KMer(sequence=seq, left=left, right=right)


def slide_kmer(qseqid: str, kmer: KMer) -> Iterator[Row]:
    """
    Compare kmer against k-mers from its surrounding context and yield
    match counts.

    For each k-mer extracted from the left and right contexts, count
    how many positions match.

    Parameters:
        qseqid (str): The contig identifier.
        kmer (KMer): The k-mer with its context.

    Yields:
        Row: A dictionary row with match details.
    """

    size = len(kmer.sequence)
    counter: MutableMapping[int, int] = defaultdict(int)

    for other in chain(kmers(kmer.left, size),
                       kmers(kmer.right, size),
                       ):
        matches = sum(x == y for x, y in zip(kmer.sequence, other.sequence))
        if matches > 0:
            counter[matches] += 1

    for matches, occurrences in counter.items():
        ret: Row = {"qseqid": qseqid,
                    "size": size,
                    "matches": matches,
                    "occurrences": occurrences,
                    "kmer": kmer.sequence,
                    }
        yield ret


def process_contig(contig: Contig, max_kmer: int) -> Iterator[Row]:
    """
    Process a single contig, compute k-mer statistics for k=1..max_kmer.

    Parameters:
        contig (Contig): The contig to process.
        max_kmer (int): The maximum k-mer size to consider.

    Yields:
        Row: The computed CSV row data for each unique k-mer.
    """

    qseqid = (contig.name and str(contig.name)) or ''
    for size in range(1, max_kmer + 1):
        all_kmers = kmers(contig.seq, size)
        deduplicated = {x.sequence: x for x in all_kmers}.values()
        for kmer in sorted(deduplicated, key=lambda x: x.sequence):
            yield from slide_kmer(qseqid, kmer)


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

    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w") as output_file:
        writer = csv.DictWriter(output_file, fieldnames=FIELDNAMES)
        writer.writeheader()
        for contig in contigs:
            for row in process_contig(contig, max_kmer):
                writer.writerow(row)


def parse_arguments(argv: Sequence[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Concatenate inputs and write to an output file.")

    parser.add_argument('input', type=Path,
                        help='Input FASTA file.')
    parser.add_argument('output', type=Path,
                        help='Input CSV file.')
    parser.add_argument('--max', type=int, default=10,
                        help='Largest k-mer to be analyzed.')

    return parser.parse_args(argv)


def main(argv: Sequence[str]) -> int:
    try:
        args = parse_arguments(argv)
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
