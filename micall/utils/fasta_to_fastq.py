#! /usr/bin/env python

"""
This module generates a random FASTQ file of simulated reads derived
from a reference FASTA.  Each read is a random substring (without
mutations) from the reference.  When the random reads are assembled
(given sufficient overlap and coverage), they should produce the
original FASTA.
"""

import argparse
import random
import sys
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord, Seq
from typing import Sequence, Iterator
from pathlib import Path

MAX_QUALITY = 40


def simulate_reads(reference: Seq,
                   n_reads: int,
                   is_reversed: bool,
                   min_length: int,
                   max_length: int,
                   extract_num: int,
                   ) -> Iterator[SeqRecord]:

    """
    Generate random reads records from the given reference sequence.

    Args:
      reference: The reference sequence.
      n_reads: Total number of reads to generate.
      is_reversed: Whether this is reverse complement reads.
      min_length: Minimum length of each read.
      max_length: Maximum length of each read.
      extract_num: Extraction number, a metadata component.

    Returns:
      Bio.SeqRecord objects representing FASTQ reads.

    Each record's sequence is a random substring (of a random length
    between min_length and max_length) of the reference sequence. The
    quality string is a dummy string, consisting of the letter 'I'
    repeated for the length of the read (Phred score 40).
    """

    if is_reversed:
        reference = reference.reverse_complement()

    ref_len = len(reference)
    file_num = 2 if is_reversed else 1

    for i in range(n_reads):
        # Choose a read length uniformly between min_length and max_length.
        read_length = random.randint(min_length, max_length)

        # Choose a random start index ensuring the read fits within
        # the reference.
        if ref_len - read_length <= 0:
            start = 0
        else:
            start = random.randrange(0, ref_len - read_length + 1)

        end = start + read_length
        read_seq_seq = reference[start:end]
        read_seq_str = str(read_seq_seq)
        read_seq = Seq(read_seq_str)

        # Create a dummy quality list (here "max" for each base) and
        # then convert that into a FASTQ quality string.
        # Bio.SeqIO.write when given a SeqRecord with
        # letter_annotations["phred_quality"] writes in FASTQ.
        qualities = [MAX_QUALITY] * read_length

        y_coord = i + 1
        record_id = f"""\
M01234:01:000000000-AAAAA:1:1101:{extract_num}:{y_coord:04d} {file_num}:N:0:1
""".strip()
        description = f"start={start} length={read_length}"

        record = SeqRecord(read_seq, id=record_id, description=description)
        record.letter_annotations["phred_quality"] = qualities

        yield record


def generate_fastq(fasta: Path,
                   fastq: Path,
                   n_reads: int,
                   is_reversed: bool,
                   min_length: int,
                   max_length: int,
                   extract_num: int,
                   ) -> None:

    """
    Read a reference from a FASTA file and write a FASTQ file with
    simulated random reads.

    Args:
      fasta: Path to the input FASTA file containing one or
                   more reference sequences.
      fastq: Path to the output FASTQ file.
      n_reads: Total number of reads to generate.
      is_reversed: Whether this is reverse complement reads.
      min_length: Minimum length of each read.
      max_length: Maximum length of each read.
      extract_num: Extraction number, a metadata component.
    """

    with open(fasta, "r") as fasta_handle, \
         open(fastq, "w") as fastq_handle:
        for record in SeqIO.parse(fasta_handle, "fasta"):
            reference: Seq = record.seq
            simulated_reads = simulate_reads(reference=reference,
                                             n_reads=n_reads,
                                             is_reversed=is_reversed,
                                             min_length=min_length,
                                             max_length=max_length,
                                             extract_num=extract_num,
                                             )
            SeqIO.write(simulated_reads, fastq_handle, format="fastq")


def get_parser() -> argparse.ArgumentParser:
    description = """\
Generate a random FASTQ file of simulated reads from a reference FASTA.
The simulated reads, when assembled,\
 should produce the original FASTA sequence.
"""
    p = argparse.ArgumentParser(description=description)
    p.add_argument("fasta", type=Path,
                   help="Path to the input reference FASTA file.")
    p.add_argument("fastq", type=Path,
                   help="Path to the output FASTQ file.")
    p.add_argument("--nreads", type=int, default=1000,
                   help="Number of reads to generate.")
    p.add_argument("--reversed", action='store_true',
                   help="Generate FASTQ file for the reverse complement.")
    p.add_argument("--min_length", type=int, default=250,
                   help="Minimum length of each simulated read.")
    p.add_argument("--max_length", type=int, default=250,
                   help="Maximum length of each simulated read.")
    p.add_argument("--extract_num", type=int, default=1234,
                   help="Extraction number, a metadata component.")
    return p


def main(argv: Sequence[str]) -> int:
    parser = get_parser()
    args = parser.parse_args(argv)

    generate_fastq(fasta=args.fasta,
                   fastq=args.fastq,
                   n_reads=args.nreads,
                   is_reversed=args.reversed,
                   min_length=args.min_length,
                   max_length=args.max_length,
                   extract_num=args.extract_num,
                   )
    return 0


def entry() -> None:
    sys.exit(main(sys.argv[1:]))


if __name__ == '__main__': entry()  # noqa
