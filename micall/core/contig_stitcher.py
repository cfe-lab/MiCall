from typing import Sequence
import logging

from micall.utils.referencefull_contig_stitcher \
    import referencefull_contig_stitcher
from micall.utils.referenceless_contig_stitcher \
    import referenceless_contig_stitcher


logger = logging.getLogger(__name__)


def main(argv: Sequence[str]):
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('contigs', type=argparse.FileType('r'), help="Input CSV file with assembled contigs.")
    parser.add_argument('stitched_contigs', type=argparse.FileType('w'),
                        help="Output CSV file with stitched contigs.")
    parser.add_argument('--use-references', required=True, choices=['yes', 'no'],
                        help="Whether to use reference genomes during stitching.")
    parser.add_argument('--plot', type=argparse.FileType('w'),
                        help="Output SVG image visualizing the stitching process.")
    verbosity_group = parser.add_mutually_exclusive_group()
    verbosity_group.add_argument('--verbose', action='store_true', help='Increase output verbosity.')
    verbosity_group.add_argument('--no-verbose', action='store_true', help='Normal output verbosity.', default=True)
    verbosity_group.add_argument('--debug', action='store_true', help='Maximum output verbosity.')
    verbosity_group.add_argument('--quiet', action='store_true', help='Minimize output verbosity.')

    args = parser.parse_args(argv)

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

    if args.use_references:
        referencefull_contig_stitcher(args.contigs, args.stitched_contigs, plot_path)
    else:
        referenceless_contig_stitcher(args.contigs, args.stitched_contigs)


if __name__ == '__main__':
    import sys
    main(sys.argv[1:])
