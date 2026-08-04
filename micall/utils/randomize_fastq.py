#! /usr/bin/env python

import argparse
import random
import sys
from pathlib import Path
from typing import Sequence, Tuple, TextIO, Iterator

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

NUCLEOTIDES = "ACGT"
SNP_POSSIBILITIES = {
    base: tuple(nuc for nuc in NUCLEOTIDES if nuc != base)
    for base in NUCLEOTIDES
}


def introduce_errors(seq: str,
                     qualities: Sequence[int],
                     subst_rate: float,
                     ins_rate: float,
                     del_rate: float,
                     ins_quality: int,
                     rng: random.Random,
                     ) -> Tuple[str, Sequence[int]]:

    """
    Introduce substitution, insertion, and deletion errors into a sequence.

    Args:
      seq: Original sequence string.
      qualities: List of phred quality values corresponding to each base.
      subst_rate: Per-base probability of substituting the base.
      ins_rate: Per-base probability (after processing a base)
                of inserting an extra nucleotide.
      del_rate: Per-base probability of deleting the base.
      ins_quality: Quality score for any inserted bases.

    Returns:
      A tuple (new_seq, new_qualities) where new_seq is the
      error-introduced sequence and new_qualities is a list of quality
      scores for the new sequence.
    """

    new_seq_chars = []
    new_quals = []

    for base, q in zip(seq, qualities):
        # First, decide whether to drop (delete) this base.
        if rng.random() < del_rate:
            # Base is dropped (no substitution or insertion for this
            # base).  You could decide to also sometimes insert an
            # extra unwanted base even if deletion occurs.
            continue

        # Otherwise, decide whether to substitute the base.
        if rng.random() < subst_rate:
            # Choose a new nucleotide that is different from the original.
            new_base = rng.choice(SNP_POSSIBILITIES[base])
        else:
            new_base = base

        new_seq_chars.append(new_base)
        new_quals.append(q)

        # Now, decide whether to insert an extra (random) nucleotide
        # *after* this base.
        while rng.random() < ins_rate:
            inserted_base = rng.choice(NUCLEOTIDES)
            new_seq_chars.append(inserted_base)
            new_quals.append(ins_quality)

    return "".join(new_seq_chars), new_quals


def process_records(input_handle: TextIO,
                    subst_rate: float,
                    ins_rate: float,
                    del_rate: float,
                    ins_quality: int,
                    rng: random.Random,
                    ) -> Iterator[SeqRecord]:

    for record in SeqIO.parse(input_handle, "fastq"):
        # Get original sequence and quality digits.
        original_seq = str(record.seq)
        original_quals = record.letter_annotations["phred_quality"]

        # Introduce errors into this read.
        new_seq_str, new_quals = introduce_errors(original_seq,
                                                  original_quals,
                                                  subst_rate,
                                                  ins_rate,
                                                  del_rate,
                                                  ins_quality,
                                                  rng)

        # Create a new SeqRecord with updated sequence and quality.
        new_record = SeqRecord(Seq(new_seq_str),
                               id=record.id,
                               name=record.name,
                               description=record.description)
        new_record.letter_annotations["phred_quality"] = new_quals
        yield new_record


def process_fastq(in_fastq: Path,
                  out_fastq: Path,
                  subst_rate: float,
                  ins_rate: float,
                  del_rate: float,
                  ins_quality: int,
                  rng: random.Random,
                  ) -> None:

    """
    Process an input FASTQ file, introducing errors into each read,
    and write the modified reads to an output FASTQ file.

    Args:
      in_fastq: Path to the input FASTQ file.
      out_fastq: Path to the output FASTQ file.
      subst_rate: Substitution error rate.
      ins_rate: Insertion error rate.
      del_rate: Deletion error rate.
      ins_quality: Base quality score for inserted bases.
    """

    with open(in_fastq, "r") as input_handle, \
         open(out_fastq, "w") as output_handle:
        new_records = process_records(input_handle=input_handle,
                                      subst_rate=subst_rate,
                                      ins_rate=ins_rate,
                                      del_rate=del_rate,
                                      ins_quality=ins_quality,
                                      rng=rng,
                                      )
        SeqIO.write(new_records, output_handle, "fastq")


def get_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Introduce errors into a FASTQ file to generate test data for error-prone reads.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("in_fastq", type=Path,
                        help="Input FASTQ file containing original reads.")
    parser.add_argument("out_fastq", type=Path,
                        help="Output FASTQ file that will contain reads with errors.")
    parser.add_argument("--subst_rate", type=float, default=1/500,
                        help="Per-base substitution error rate.")
    parser.add_argument("--ins_rate", type=float, default=1/5000,
                        help="Per-base insertion error rate.")
    parser.add_argument("--del_rate", type=float, default=1/5000,
                        help="Per-base deletion error rate.")
    parser.add_argument("--ins_quality", type=int, default=20,
                        help="Quality score assigned to any inserted bases.")
    parser.add_argument("--seed", type=int, default=42,
                        help="Random seed for reproducibility.")
    return parser


def main(argv: Sequence[str]) -> int:
    parser = get_parser()
    args = parser.parse_args(argv)
    rng = random.Random(args.seed)
    process_fastq(in_fastq=args.in_fastq,
                  out_fastq=args.out_fastq,
                  subst_rate=args.subst_rate,
                  ins_rate=args.ins_rate,
                  del_rate=args.del_rate,
                  ins_quality=args.ins_quality,
                  rng=rng,
                  )
    return 0


def entry() -> None:
    sys.exit(main(sys.argv[1:]))


if __name__ == '__main__': entry()  # noqa
