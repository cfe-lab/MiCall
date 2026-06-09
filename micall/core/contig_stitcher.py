import sys
from typing import Sequence, Dict
import logging
from pathlib import Path

import numpy as np
import micall.utils.referencefull_contig_stitcher as referencefull
import micall.utils.referenceless_contig_stitcher as referenceless

# Re-export for backward compatibility with tests
from micall.utils.referencefull_contig_stitcher import read_remap_counts, read_contigs  # noqa: F401

logger = logging.getLogger(__name__)


def main(argv: Sequence[str]) -> int:
    import argparse

    head_parser = argparse.ArgumentParser(description='Contig stitcher.', add_help=False)
    head_parser.add_argument('mode', choices=['with-references', 'without-references'])
    head_parser.add_argument('stitcher_arguments', nargs=argparse.REMAINDER)

    if len(argv) > 0 and argv[0] in ['--help', '-h']:
        head_parser.print_help()
        return 0

    head_args = head_parser.parse_args(argv)
    if head_args.mode == 'with-references':
        parser = argparse.ArgumentParser(description='Reference-based contig stitcher.')
        parser.add_argument('contigs', type=argparse.FileType('r'), help="Input CSV file with assembled contigs.")
        parser.add_argument('stitched_contigs', type=argparse.FileType('w'),
                            help="Output CSV file with stitched contigs.")
        parser.add_argument('--remap-counts', type=argparse.FileType('r'), default=None,
                            help="Optional input CSV file with remap read counts.")
        parser.add_argument('--plot', type=argparse.FileType('w'),
                            help="Output SVG image visualizing the stitching process.")
    else:
        parser = argparse.ArgumentParser(description='Reference-less contig stitcher.')
        parser.add_argument('contigs', type=argparse.FileType('r'), help="Input FASTA file with assembled contigs.")
        parser.add_argument('stitched_contigs', type=argparse.FileType('w'),
                            help="Output FASTA file with stitched contigs.")
        parser.add_argument('--fastq1', type=Path, default=None,
                            help='Forward reads FASTQ for coverage validation (plain or .gz).')
        parser.add_argument('--fastq2', type=Path, default=None,
                            help='Reverse reads FASTQ for coverage validation (plain or .gz).')
        parser.add_argument('--minimum-read-depth', type=int, default=1,
                            help='Minimum read coverage to accept an overlap. '
                                 '0 disables coverage validation. (default: 1)')
        parser.add_argument('--read-length', type=int, default=150,
                            help='Read length used for the boundary margin. '
                                 'Coverage is checked overlap +/- read_length. (default: 150)')

    parser.prog = ' '.join([Path(__file__).stem, head_args.mode])

    verbosity_group = parser.add_mutually_exclusive_group()
    verbosity_group.add_argument('--verbose', action='store_true', help='Increase output verbosity.')
    verbosity_group.add_argument('--no-verbose', action='store_true', help='Normal output verbosity.', default=True)
    verbosity_group.add_argument('--debug', action='store_true', help='Maximum output verbosity.')
    verbosity_group.add_argument('--debug2', action='store_true', help='Even more debug messages.')
    verbosity_group.add_argument('--quiet', action='store_true', help='Minimize output verbosity.')

    args = parser.parse_args(head_args.stitcher_arguments)

    if args.quiet:
        logger.setLevel(logging.ERROR)
    elif args.verbose:
        logger.setLevel(logging.INFO)
    elif args.debug:
        logger.setLevel(logging.DEBUG)
    elif args.debug2:
        logger.setLevel(logging.DEBUG - 1)
    else:
        logger.setLevel(logging.WARN)

    logging.basicConfig(level=logger.level)

    if head_args.mode == 'with-references':
        plot_path = Path(args.plot.name) if args.plot is not None else None
        remap_counts = args.remap_counts
        referencefull.logger = logger
        referencefull.referencefull_contig_stitcher(
            args.contigs, args.stitched_contigs, plot_path, remap_counts)
    else:
        referenceless.logger = logger
        coverage_data: Dict[str, np.ndarray] = {}
        if args.fastq1 is not None and args.fastq2 is not None:
            contig_records = tuple(referenceless.read_contigs(args.contigs))
            args.contigs.seek(0)
            coverage_data = referenceless.compute_coverage_from_fastqs(
                args.fastq1, args.fastq2, contig_records,
            )
        referenceless.referenceless_contig_stitcher(
            args.contigs, args.stitched_contigs,
            coverage_data=coverage_data,
            minimum_read_depth=args.minimum_read_depth,
            read_length=args.read_length,
        )

    return 0


def cli() -> None:
    try:
        sys.exit(main(sys.argv[1:]))
    except KeyboardInterrupt:
        sys.exit(1)


if __name__ == '__main__': cli()  # noqa
