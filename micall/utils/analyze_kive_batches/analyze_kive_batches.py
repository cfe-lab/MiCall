#! /usr/bin/env python

import argparse
import sys
from typing import Sequence
from pathlib import Path
import logging

from micall.utils.dir_path import DirPath
from micall.utils.user_error import UserError
from .logger import logger

from .make_stats import make_stats
from .get_batch import get_batch
from .run_all import run_all
from .combine_batches_runs import combine_batches_runs
from .combine_runs_stats import combine_runs_stats
from .combine_runs_overlaps import combine_runs_overlaps
from .extract_run_ids import extract_run_ids
from .aggregate_runs_stats import aggregate_runs_stats
from .aggregate_runs_overlaps import aggregate_runs_overlaps
from .stitch_contigs import stitch_contigs
from .calculate_exact_coverage_file import calculate_exact_coverage_file
from .make_properties import make_properties
from .join_tables import join_tables
from .diff_samples_of_two_apps import diff_samples_of_two_apps


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

    sub = mode_parsers.add_parser("diff_samples_of_two_apps", help="Computes differences between samples of two apps.")
    sub.add_argument("--input", type=Path, required=True,
                     help="Input CSV file.")
    sub.add_argument("--app1", type=str, required=True,
                     help="Name of the first app.")
    sub.add_argument("--app2", type=str, required=True,
                     help="Name of the second app.")
    sub.add_argument("--output", type=Path, required=True,
                     help="Target CSV file where to put the results to.")

    sub = mode_parsers.add_parser("get-batch", help="Downloads a batch info.")
    sub.add_argument("--batch", type=str, required=True,
                     help="The name of the batch to download the runs info for.")
    sub.add_argument("--target", type=Path, required=True,
                     help="Target file where to put the runs info to.")

    sub = mode_parsers.add_parser("combine-batches-runs", help="Extract batches run infos and combine them.")
    sub.add_argument("--batches", type=Path, required=True, nargs=argparse.ONE_OR_MORE,
                     help="The downloaded batches files.")
    sub.add_argument("--target", type=Path, required=True,
                     help="Target file where to put the runs info to.")

    sub = mode_parsers.add_parser("combine-runs-stats", help="Combine all stats.json files into one.")
    sub.add_argument("--root", type=dir_path, required=True,
                     help="Root directory for all output subdirectories.")
    sub.add_argument("--runs-txt", type=Path, required=True,
                     help="The txt file with all the run ids in it.")
    sub.add_argument("--target", type=Path, required=True,
                     help="Target file where to put the combine stats to.")

    sub = mode_parsers.add_parser("combine-runs-overlaps", help="Combine all stats.json:overlaps data into one file.")
    sub.add_argument("--root", type=dir_path, required=True,
                     help="Root directory for all output subdirectories.")
    sub.add_argument("--runs-txt", type=Path, required=True,
                     help="The txt file with all the run ids in it.")
    sub.add_argument("--target", type=Path, required=True,
                     help="Target file where to put the combine overlaps to.")

    sub = mode_parsers.add_parser("make-stats")
    sub.add_argument("--input", type=Path, required=True,
                     help="Input JSON file with the run info.")
    sub.add_argument("--output", type=Path, required=True,
                     help="Output stats file.")

    sub = mode_parsers.add_parser("stitch-contigs")
    sub.add_argument("--info-file", type=Path, required=True,
                     help="Input JSON file with the run info.")
    sub.add_argument("--output", type=Path, required=True,
                     help="Output file.")

    sub = mode_parsers.add_parser("calculate-exact-coverage")
    sub.add_argument("--info-file", type=Path, required=True,
                     help="Input JSON file with the run info.")
    sub.add_argument("--output", type=Path, required=True,
                     help="Output CSV file with exact coverage data.")

    sub = mode_parsers.add_parser("aggregate-runs-stats")
    sub.add_argument("--input", type=Path, required=True,
                     help="Input CSV file with the run stats.")
    sub.add_argument("--output", type=Path, required=True,
                     help="Output stats file.")

    sub = mode_parsers.add_parser("aggregate-runs-overlaps")
    sub.add_argument("--input", type=Path, required=True,
                     help="Input CSV file with the run overlaps.")
    sub.add_argument("--output", type=Path, required=True,
                     help="Output stats file.")

    sub = mode_parsers.add_parser("make-properties")
    sub.add_argument("--input", type=Path, required=True,
                     help="Input TOML file with the apps properties.")
    sub.add_argument("--output", type=Path, required=True,
                     help="Output CSV file.")

    sub = mode_parsers.add_parser("extract-run-ids")
    sub.add_argument("--input", type=Path, required=True,
                     help="Input JSON file with the run info.")
    sub.add_argument("--output", type=Path, required=True,
                     help="Output ids file.")

    sub = mode_parsers.add_parser("join-tables")
    sub.add_argument("--inputs", type=Path, nargs=argparse.ONE_OR_MORE,
                     help="Input CSV files to union.")
    sub.add_argument("--column", type=str, required=True,
                     help="The column that will serve as the index.")
    sub.add_argument("--output", type=Path, required=True,
                     help="Output CSV file.")

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
    elif args.subcommand == 'diff_samples_of_two_apps':
        diff_samples_of_two_apps(input=args.input, app1=args.app1, app2=args.app2, output=args.output)
    elif args.subcommand == 'get-batch':
        get_batch(batch=args.batch, target=args.target)
    elif args.subcommand == 'combine-batches-runs':
        combine_batches_runs(batches=args.batches, target=args.target)
    elif args.subcommand == 'combine-runs-stats':
        combine_runs_stats(root=args.root, runs_txt=args.runs_txt, target=args.target)
    elif args.subcommand == 'combine-runs-overlaps':
        combine_runs_overlaps(root=args.root, runs_txt=args.runs_txt, target=args.target)
    elif args.subcommand == 'make-stats':
        make_stats(input=args.input, output=args.output)
    elif args.subcommand == 'stitch-contigs':
        stitch_contigs(info_file=args.info_file, output=args.output)
    elif args.subcommand == 'calculate-exact-coverage':
        calculate_exact_coverage_file(info_file=args.info_file, output=args.output)
    elif args.subcommand == 'extract-run-ids':
        extract_run_ids(input=args.input, output=args.output)
    elif args.subcommand == 'aggregate-runs-stats':
        aggregate_runs_stats(input=args.input, output=args.output)
    elif args.subcommand == 'aggregate-runs-overlaps':
        aggregate_runs_overlaps(input=args.input, output=args.output)
    elif args.subcommand == 'make-properties':
        make_properties(input=args.input, output=args.output)
    elif args.subcommand == 'join-tables':
        join_tables(inputs=args.inputs, column=args.column, output=args.output)
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
