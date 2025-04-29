#! /usr/bin/env python

import argparse
import sys
from typing import Sequence
from pathlib import Path
import logging

from micall.utils.dir_path import DirPath
from micall.utils.user_error import UserError
from .logger import logger

from .download import download
from .make_stats_1 import make_stats_1
from .get_batch import get_batch
from .run_all import run_all
from .combine_batches_runs import combine_batches_runs
from .combine_runs_stats import combine_runs_stats


def dir_path(string: str) -> DirPath:
    path = Path(string)
    if (not path.exists()) or path.is_dir():
        return DirPath(path)
    else:
        raise UserError("Path %r is not a directory.", string)


def cli_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Analyze a kive run.")
    mode_parsers = parser.add_subparsers(dest='subcommand',
                                         title='subcommands',
                                         required=True,
                                         )

    all = mode_parsers.add_parser("all", help="The main entry to this script. Runs all other subentries.")
    all.add_argument("--batches-list", type=Path, required=True,
                     help="Path to a text-file containing the list of batches to be analyzed.")
    all.add_argument("--root", type=dir_path, required=True,
                     help="Root directory for all output subdirectories.")
    all.add_argument("--properties", type=Path, required=True,
                     help="Additional properties associated with particular images.")

    get_batch = mode_parsers.add_parser("get-batch", help="The main entry to this script. Runs get_batch other subentries.")
    get_batch.add_argument("--batch", type=str, required=True,
                           help="The name of the batch to download the runs info for.")
    get_batch.add_argument("--target", type=Path, required=True,
                           help="Target file where to put the runs info to.")

    sub = mode_parsers.add_parser("combine-batches-runs", help="Extract batches run infos and combine them.")
    sub.add_argument("--batches", type=Path, required=True, nargs=argparse.ONE_OR_MORE,
                     help="The downloaded batches files.")
    sub.add_argument("--target", type=Path, required=True,
                     help="Target file where to put the runs info to.")

    sub = mode_parsers.add_parser("combine-runs-stats", help="Combine all stats.json files into one.")
    sub.add_argument("--root", type=dir_path, required=True,
                     help="Root directory for all output subdirectories.")
    sub.add_argument("--runs-json", type=Path, required=True,
                     help="The big json file with all the run infos.")
    sub.add_argument("--target", type=Path, required=True,
                     help="Target file where to put the combine stats to.")

    download = mode_parsers.add_parser("download")
    download.add_argument("--json-file", type=Path, required=True,
                          help="The big json file with all the run infos.")
    download.add_argument("--root", type=dir_path, required=True,
                          help="Root directory for all output subdirectories.")

    make_stats_1 = mode_parsers.add_parser("make-stats-1")
    make_stats_1.add_argument("--input", type=Path, required=True,
                              help="Input json file with the run info.")
    make_stats_1.add_argument("--output", type=Path, required=True,
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


def main_typed(subcommand: str, args: argparse.Namespace) -> None:
    if args.subcommand == 'all':
        run_all(batches_list=args.batches_list, root=args.root, properties=args.properties)
    elif args.subcommand == 'get-batch':
        get_batch(batch=args.batch, target=args.target)
    elif args.subcommand == 'combine-batches-runs':
        combine_batches_runs(batches=args.batches, target=args.target)
    elif args.subcommand == 'combine-runs-stats':
        combine_runs_stats(root=args.root, runs_json=args.runs_json, target=args.target)
    elif args.subcommand == 'download':
        download(json_file=args.json_file, root=args.root)
    elif args.subcommand == 'make-stats-1':
        make_stats_1(input=args.input, output=args.output)
    else:
        raise UserError("Unrecognized subcommand %r.", args.subcommand)


def main(argv: Sequence[str]) -> int:
    args = parse_args(argv)
    subcommand: str = args.subcommand

    if args.quiet:
        logger.setLevel(logging.ERROR)
    elif args.verbose:
        logger.setLevel(logging.INFO)
    elif args.debug:
        logger.setLevel(logging.DEBUG)
    else:
        logger.setLevel(logging.WARN)

    try:
        main_typed(subcommand, args)
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
