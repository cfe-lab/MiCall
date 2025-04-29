#! /usr/bin/env python

import argparse
import sys
from typing import Sequence
from pathlib import Path
import logging

from micall.utils.dir_path import DirPath
from micall.utils.user_error import UserError
from micall.utils.analyze_kive_batches.logger import logger

import micall.utils.analyze_kive_batches.download
import micall.utils.analyze_kive_batches.make_stats_1


def dir_path(string: str) -> DirPath:
    path = Path(string)
    if (not path.exists()) or path.is_dir():
        return DirPath(path)
    else:
        raise UserError("Path %r is not a directory.", string)


def cli_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser("Analyze a kive run.")
    mode_parsers = parser.add_subparsers(dest='subcommand',
                                         title='subcommands',
                                         required=True,
                                         )

    download = mode_parsers.add_parser("download")
    download.add_argument("--json-file", type=Path,
                          help="The big json file with all the run infos.")
    download.add_argument("--root", type=dir_path,
                          help="Root directory for all output subdirectories.")

    make_stats_1 = mode_parsers.add_parser("make-stats-1")
    make_stats_1.add_argument("--input", type=dir_path,
                              help="Input directory with all the input/output files of a MiCall run.")
    make_stats_1.add_argument("--output", type=Path,
                              help="Output stats file.")

    verbosity_group = parser.add_mutually_exclusive_group()
    verbosity_group.add_argument('--verbose', action='store_true', help='Increase output verbosity.')
    verbosity_group.add_argument('--no-verbose', action='store_true', help='Normal output verbosity.', default=True)
    verbosity_group.add_argument('--debug', action='store_true', help='Maximum output verbosity.')
    verbosity_group.add_argument('--quiet', action='store_true', help='Minimize output verbosity.')

    return parser


def parse_args(argv: Sequence[str]) -> argparse.Namespace:
    parser = cli_parser()
    return parser.parse_args(argv)


def main_typed(args: argparse.Namespace) -> None:
    if args.subcommand == 'download':
        micall.utils.analyze_kive_batches.download.main_typed(json_file=args.json_file, root=args.root)
    elif args.subcommand == 'make-stats-1':
        micall.utils.analyze_kive_batches.make_stats_1.main_typed(input=args.input, output=args.output)
    else:
        raise UserError("Unrecognized subcommand %r.", args.subcommand)


def main(argv: Sequence[str]) -> int:
    args = parse_args(argv)

    if args.quiet:
        logger.setLevel(logging.ERROR)
    elif args.verbose:
        logger.setLevel(logging.INFO)
    elif args.debug:
        logger.setLevel(logging.DEBUG)
    else:
        logger.setLevel(logging.WARN)

    try:
        main_typed(args)
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


if __name__ == '__main__': entry()  # noqa
