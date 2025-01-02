from typing import Sequence
import logging
from pathlib import Path

from micall.utils.referencefull_contig_stitcher \
    import referencefull_contig_stitcher
from micall.utils.referenceless_contig_stitcher \
    import referenceless_contig_stitcher


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
        parser.add_argument('--plot', type=argparse.FileType('w'),
                            help="Output SVG image visualizing the stitching process.")
    else:
        parser = argparse.ArgumentParser(description='Reference-less contig stitcher.')
        parser.add_argument('contigs', type=argparse.FileType('r'), help="Input FASTA file with assembled contigs.")
        parser.add_argument('stitched_contigs', type=argparse.FileType('w'),
                            help="Output FASTA file with stitched contigs.")

    parser.prog = ' '.join([Path(__file__).stem, head_args.mode])
    verbosity_group = parser.add_mutually_exclusive_group()
    verbosity_group.add_argument('--verbose', action='store_true', help='Increase output verbosity.')
    verbosity_group.add_argument('--no-verbose', action='store_true', help='Normal output verbosity.', default=True)
    verbosity_group.add_argument('--debug', action='store_true', help='Maximum output verbosity.')
    verbosity_group.add_argument('--quiet', action='store_true', help='Minimize output verbosity.')

    args = parser.parse_args(head_args.stitcher_arguments)

    if args.quiet:
        logger.setLevel(logging.ERROR)
    elif args.verbose:
        logger.setLevel(logging.INFO)
    elif args.debug:
        logger.setLevel(logging.DEBUG)
    else:
        logger.setLevel(logging.WARN)

    logging.basicConfig(level=logger.level)

    plot_path = args.plot.name if args.plot is not None else None

    if head_args.mode == 'with-references':
        referencefull_contig_stitcher(args.contigs, args.stitched_contigs, plot_path)
    else:
        referenceless_contig_stitcher(args.contigs, args.stitched_contigs)

    return 0


if __name__ == '__main__':
    import sys
    sys.exit(main(sys.argv[1:]))
