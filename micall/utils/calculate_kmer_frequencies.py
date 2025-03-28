#! /usr/bin/env python

"""
This script computes kmer frequencies of given sequences.
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
    qseqid: str
    size: int
    matches: int
    occurrences: int
    kmer: str


FIELDNAMES = tuple(Row.__annotations__.keys())


class UserError(ValueError):
    def __init__(self, fmt: str, *fmt_args: object):
        self.fmt = fmt
        self.fmt_args = fmt_args
        self.code = 1


@dataclass(frozen=True)
class KMer:
    sequence: str
    left: str
    right: str


def read_contigs(input_fasta: TextIO) -> Iterator[Contig]:
    for record in SeqIO.parse(input_fasta, "fasta"):
        yield Contig(name=record.name, seq=str(record.seq))


def kmers(sequence: str, size: int) -> Iterator[KMer]:
    for i in range(len(sequence) - size):
        end = i+size
        seq = sequence[i:end]
        left = sequence[:i]
        right = sequence[end:]
        yield KMer(sequence=seq, left=left, right=right)


def slide_kmer(qseqid: str, kmer: KMer) -> Iterator[Row]:
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
    qseqid = (contig.name and str(contig.name)) or ''
    for size in range(1, max_kmer + 1):
        all_kmers = kmers(contig.seq, size)
        deduplicated = {x.sequence: x for x in all_kmers}.values()
        for kmer in sorted(deduplicated, key=lambda x: x.sequence):
            yield from slide_kmer(qseqid, kmer)


def main_typed(input: Path, output: Path, max_kmer: int) -> int:

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

    return 0


def main(argv: Sequence[str]) -> int:
    parser = argparse.ArgumentParser(
        description="Concatenate inputs and write to an output file.")

    parser.add_argument('input', type=Path,
                        help='Input FASTA file.')
    parser.add_argument('output', type=Path,
                        help='Input CSV file.')
    parser.add_argument('--max', type=int, default=10,
                        help='Largest k-mer to be analyzed.')

    args = parser.parse_args(argv)
    return main_typed(args.input, args.output, args.max)


def entry() -> None:
    try:
        rc = main(sys.argv[1:])
        logger.debug("Done.")
    except BrokenPipeError:
        logger.debug("Broken pipe.")
        rc = 1
    except KeyboardInterrupt:
        logger.debug("Interrupted.")
        rc = 1
    except UserError as e:
        logger.fatal(e.fmt, *e.fmt_args)
        rc = e.code

    sys.exit(rc)


if __name__ == "__main__": entry()  # noqa
