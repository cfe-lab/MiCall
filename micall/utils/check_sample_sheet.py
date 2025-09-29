#!/usr/bin/env python3
"""
Script to check sample name consistency between sample sheet and FASTQ files.
"""

import argparse
import logging
import sys
from pathlib import Path
from typing import Iterable

from micall.utils.sample_sheet_parser import read_sample_sheet_and_overrides

logger = logging.getLogger(__name__)


def check_sample_name_consistency(
    sample_sheet_path: Path, fastq_file_names: Iterable[str], run_path: Path
):
    """
    Check FASTQ file recognition ratio against sample sheet.

    Prints warning when there are fewer or equal recognized FASTQ files
    compared to unrecognized ones. Recognized FASTQ are those that have
    corresponding entries in the sample sheet.

    :param sample_sheet_path: Path to the SampleSheet.csv file
    :param fastq_file_names: List of FASTQ file names
    :param run_path: Path to the run folder for logging
    """
    logger.info("Checking sample name consistency for run: %s", run_path)
    logger.debug("Sample sheet path: %s", sample_sheet_path)
    logger.debug("FASTQ files to check: %s", list(fastq_file_names))

    # Extract sample names from the sample sheet
    logger.debug("Reading sample sheet: %s", sample_sheet_path)
    run_info = read_sample_sheet_and_overrides(sample_sheet_path)
    data_split = run_info.get("DataSplit")
    if data_split is None:
        raise RuntimeError(f"Missing 'DataSplit' section in {sample_sheet_path}")

    assert hasattr(data_split, "__iter__")
    if not hasattr(data_split, "__iter__"):
        raise RuntimeError(f"Invalid 'DataSplit' section in {sample_sheet_path}")

    # Get filenames from sample sheet DataSplit (these are the expected FASTQ file names)
    sheet_filenames = set()
    for row in data_split:
        filename = row.get("filename", "")
        if filename:
            # Store the trimmed version that matches what find_groups() uses
            trimmed_filename = "_".join(filename.split("_")[:2])
            sheet_filenames.add(trimmed_filename)

    logger.debug("Sample sheet contains %d unique sample prefixes: %s", len(sheet_filenames), sorted(sheet_filenames))

    # Extract trimmed names from FASTQ files (same logic as find_groups())
    fastq_trimmed_names = set()
    for file_name in fastq_file_names:
        # This matches the logic in find_groups: '_'.join(file_name.split('_')[:2])
        trimmed_name = "_".join(file_name.split("_")[:2])
        fastq_trimmed_names.add(trimmed_name)

    logger.debug("FASTQ files contain %d unique sample prefixes: %s", len(fastq_trimmed_names), sorted(fastq_trimmed_names))

    # Identify recognized vs unrecognized FASTQ files
    recognized_fastq = fastq_trimmed_names & sheet_filenames  # intersection
    unrecognized_fastq = fastq_trimmed_names - sheet_filenames  # FASTQ only

    logger.info("Found %d recognized and %d unrecognized FASTQ files", len(recognized_fastq), len(unrecognized_fastq))

    if len(recognized_fastq) < len(unrecognized_fastq):
        unrecognized_list = ",".join(sorted(unrecognized_fastq))
        warning_message = f"""\
Large number of unrecognized FASTQ files in run folder {run_path}.
There are {len(recognized_fastq)} recognized FASTQ files.
And {len(unrecognized_fastq)} unrecognized: {unrecognized_list}."""
        logger.warning("%s", warning_message)
        print(warning_message, file=sys.stderr)
    else:
        logger.debug(
            "FASTQ recognition ratio acceptable: %d recognized, %d unrecognized in %s",
            len(recognized_fastq),
            len(unrecognized_fastq),
            run_path,
        )
        print(
            f"FASTQ recognition ratio acceptable: {len(recognized_fastq)} recognized, {len(unrecognized_fastq)} unrecognized in {run_path}"
        )


def main():
    parser = argparse.ArgumentParser(
        description="Check sample name consistency between sample sheet and FASTQ files."
    )
    parser.add_argument(
        "sample_sheet", type=Path, help="Path to the SampleSheet.csv file"
    )
    parser.add_argument(
        "--run-path",
        type=Path,
        help="Path to the run folder (default: inferred from sample_sheet)",
    )
    parser.add_argument(
        "--fastq-files",
        nargs="+",
        help="FASTQ file names (default: inferred from run_path)",
    )

    verbosity_group = parser.add_mutually_exclusive_group()
    verbosity_group.add_argument('--verbose', action='store_true', help='Increase output verbosity.')
    verbosity_group.add_argument('--no-verbose', action='store_true', help='Normal output verbosity.', default=True)
    verbosity_group.add_argument('--debug', action='store_true', help='Maximum output verbosity.')
    verbosity_group.add_argument('--quiet', action='store_true', help='Minimize output verbosity.')

    args = parser.parse_args()

    # Set logging level based on verbosity arguments
    if args.quiet:
        logger.setLevel(logging.ERROR)
    elif args.verbose:
        logger.setLevel(logging.INFO)
    elif args.debug:
        logger.setLevel(logging.DEBUG)
    else:
        logger.setLevel(logging.WARN)

    logging.basicConfig(level=logger.level, format='%(levelname)s: %(message)s')

    # Infer run_path from sample_sheet if not provided
    run_path = args.run_path
    if run_path is None:
        run_path = args.sample_sheet.parent
        logger.debug("Inferred run_path from sample_sheet: %s", run_path)

    # Infer fastq_files from run_path if not provided
    fastq_file_names = args.fastq_files
    if fastq_file_names is None:
        # Try standard MiSeq structure first
        base_calls_path = run_path / "Data" / "Intensities" / "BaseCalls"
        logger.debug("Looking for FASTQ files in standard MiSeq location: %s", base_calls_path)
        if base_calls_path.exists():
            fastq_files = list(base_calls_path.glob("*_R1_*.fastq.gz"))
            fastq_file_names = [f.name for f in fastq_files]
            logger.debug("Found %d FASTQ files in standard location", len(fastq_files))
        else:
            # Fall back to looking for .fastq files directly in run_path (for tests)
            logger.debug("Standard location not found, checking run_path directly: %s", run_path)
            fastq_files = list(run_path.glob("*_R1_*.fastq"))
            if fastq_files:
                fastq_file_names = [f.name for f in fastq_files]
                logger.debug("Found %d FASTQ files in run_path", len(fastq_files))
            else:
                logger.error("No FASTQ files found in %s or %s", base_calls_path, run_path)
                print(
                    f"Error: No FASTQ files found in {base_calls_path} or {run_path}",
                    file=sys.stderr,
                )
                sys.exit(1)
    else:
        logger.debug("Using explicitly provided FASTQ files: %s", fastq_file_names)

    logger.info("Starting consistency check with %d FASTQ files", len(fastq_file_names))

    try:
        check_sample_name_consistency(args.sample_sheet, fastq_file_names, run_path)
        logger.info("Consistency check completed successfully")
    except Exception as e:
        logger.error("Error during consistency check: %s", e)
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
