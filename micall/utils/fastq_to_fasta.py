#! /usr/bin/env python

import sys
import argparse
import logging
from pathlib import Path
from typing import Sequence
from Bio import SeqIO
from micall.utils.user_error import UserError


# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
logger = logging.getLogger(__name__)


class ModuleError(UserError):
    """Base exception for errors in this module."""


class InputNotFound(ModuleError):
    def __init__(self, source_fastq: Path):
        fmt = "Input FASTQ file does not exist %r."
        super().__init__(fmt, str(source_fastq))


class OutputDirectoryError(ModuleError):
    def __init__(self, dir_path: Path):
        fmt = "Could not create output directory %r."
        super().__init__(fmt, str(dir_path))


def main_typed(source_fastq: Path, target_fasta: Path) -> None:
    """
    Convert a FASTQ file to a FASTA file using Biopython.

    Reads sequence records from source_fastq in FASTQ format and writes them
    to target_fasta in FASTA format.
    Raises custom exceptions for various failure modes.
    """

    if not source_fastq.exists():
        raise InputNotFound(source_fastq)

    # Ensure that the output directory exists.
    try:
        target_fasta.parent.mkdir(parents=True, exist_ok=True)
        logger.debug("Output directory verified: %r.", str(target_fasta.parent))
    except Exception as e:
        raise OutputDirectoryError(target_fasta.parent) from e

    with source_fastq.open() as reader:
        records = SeqIO.parse(reader, "fastq")
        logger.debug("Loaded input records.")
        SeqIO.write(records, target_fasta, "fasta")
        logger.debug("Conversion complete.")


def parse_arguments(argv: Sequence[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Convert sequence records from FASTQ format to FASTA format."
    )

    parser.add_argument(
        "fastq_in",
        type=Path,
        help="Input FASTQ file to read sequences from."
    )
    parser.add_argument(
        "fasta_out",
        type=Path,
        help="Output FASTA file to write sequences to."
    )

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
        main_typed(args.fastq_in, args.fasta_out)
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


def cli() -> None:
    sys.exit(main(sys.argv[1:]))


if __name__ == "__main__": cli()  # noqa
