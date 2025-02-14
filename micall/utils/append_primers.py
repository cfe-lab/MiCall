#!/usr/bin/env python

"""
Append configurable primer sequences to every FASTA sequence
in an input file.
"""

import argparse
import sys
import logging
from pathlib import Path
from typing import Sequence

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

logging.basicConfig(
    level=logging.INFO,
    format="%(levelname)s: %(message)s",
)
logger = logging.getLogger(__name__)

# Default no-mixture primer sequences.
DEFAULT_FORWARD_PRIMER = "GCGCCCGAACAGGGACCTGAAAGCGAAAG"
DEFAULT_REVERSE_PRIMER = "TAAGCCTCAATAAAGCTTGCCTTGAGTGC"


class UserError(ValueError):
    def __init__(self, fmt: str, *fmt_args: object):
        self.fmt = fmt
        self.fmt_args = fmt_args
        self.code = 1


def append_primers_to_record(record: SeqRecord,
                             fwd_primer: str,
                             rev_primer: str) -> SeqRecord:

    """
    Create a new SeqRecord with primer sequences appended. The new sequence
    is the forward primer + original sequence + reverse primer.
    """

    if fwd_primer in str(record.seq):
        pos = str(record.seq).index(fwd_primer)
        logger.error("Sequence %r contains the forward primer at position %s.",
                     record.name, pos)

    if rev_primer in str(record.seq):
        pos = str(record.seq).index(rev_primer)
        logger.error("Sequence %r contains the reverse primer at position %s.",
                     record.name, pos)

    new_seq = Seq(fwd_primer + str(record.seq) + rev_primer)
    new_record = SeqRecord(
        new_seq,
        id=record.id,
        name=record.name,
        description=record.description
    )

    return new_record


def parse_arguments(argv: Sequence[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Append primer sequences to every "
                    "FASTA sequence in an input file. "
                    "Sequences are appended at both ends."
    )
    parser.add_argument("input_fasta", type=Path,
                        help="Path to the input FASTA file")
    parser.add_argument("output_fasta", type=Path,
                        help="Path to the output FASTA file")
    parser.add_argument("--forward-primer", default=DEFAULT_FORWARD_PRIMER,
                        help="The forward primer to prepend to each sequence.")
    parser.add_argument("--reverse-primer", default=DEFAULT_REVERSE_PRIMER,
                        help="The reverse primer to append to each sequence.")

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


def append_primers(input: Path,
                   output: Path,
                   forward_primer: str,
                   reverse_primer: str,
                   ) -> None:

    """
    Main routine that processes the input FASTA file,
    appends primers to each sequence,
    and writes the results to an output file.
    """

    if not input.exists() or not input.is_file():
        raise UserError("Input FASTA file does not exist: %r.", str(input))

    logger.info("Processing FASTA file: %r.", str(input))
    records = []
    for record in SeqIO.parse(str(input), "fasta"):
        new_record = append_primers_to_record(
            record, forward_primer, reverse_primer)
        records.append(new_record)

    logger.info("Primer appending complete. Writing to output file: %r",
                str(output))
    count = SeqIO.write(records, str(output), "fasta")
    logger.info("Wrote %d modified FASTA records to %r.",
                count, str(output))


def main(argv: Sequence[str]) -> int:
    args = parse_arguments(argv)
    configure_logging(args)
    append_primers(input=args.input_fasta,
                   output=args.output_fasta,
                   forward_primer=args.forward_primer,
                   reverse_primer=args.reverse_primer,
                   )
    return 0


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
