#! /usr/bin/env python3

"""
Monitor MiSeq run directories and create needsprocessing markers when runs complete.

This module monitors run directories for FASTQ file stability and creates 'needsprocessing'
marker files when a run is complete and ready for processing.

The monitoring runs in an independent daemon thread with a fixed 60-second polling interval.
"""

import argparse
import logging
from dataclasses import dataclass
import os
from pathlib import Path
import sys
from time import sleep
from typing import List, Sequence

from micall.monitor import disk_operations

logger = logging.getLogger(__name__)

# Fixed polling interval in seconds
COMPLETION_POLL_INTERVAL = 60

# Delay before restarting after crash
CRASH_RECOVERY_DELAY = 10


@dataclass
class RunInfo:
    """Information about a run directory that needs monitoring."""

    run_dir: Path
    file_count: int
    glob_pattern: str  # Relative pattern like 'Data/Intensities/BaseCalls/*.gz'


def find_unstable_runs(runs_dir: Path) -> Sequence[RunInfo]:
    """
    Find run directories without needsprocessing or processed markers.

    Searches for run directories that have FASTQ files but haven't been marked
    for processing yet. Looks in two possible locations for FASTQ files:
    1. Traditional location: Data/Intensities/BaseCalls/*.gz
    2. Newer location: Alignment*/*/Fastq/*.gz

    Args:
        runs_dir: Directory containing MiSeq run folders

    Returns:
        List of RunInfo objects for runs that need monitoring
    """
    unstable_runs: List[RunInfo] = []

    try:
        if not runs_dir.exists():
            logger.warning("Runs directory does not exist: %s", runs_dir)
            return unstable_runs

        logger.debug("Searching for runs in %s", runs_dir)

        for run_path in sorted(runs_dir.iterdir()):
            try:
                if not run_path.is_dir():
                    continue

                # Skip if already marked
                needs_processing = run_path / "needsprocessing"
                already_processed = run_path / "processed"

                if needs_processing.exists() or already_processed.exists():
                    logger.debug("Skipping %s (already marked)", run_path.name)
                    continue

                # Check for files in multiple possible locations:
                # 1. Traditional location: Data/Intensities/BaseCalls/*.gz
                basecalls_pattern = "Data/Intensities/BaseCalls/*.gz"
                basecalls_files = list(run_path.glob(basecalls_pattern))
                num_basecalls = len(basecalls_files)

                # 2. Newer location: Alignment*/*/Fastq/*.gz
                alignment_pattern = "Alignment*/*/Fastq/*.gz"
                alignment_files = list(run_path.glob(alignment_pattern))
                num_alignment = len(alignment_files)

                # Use whichever location has files (prefer BaseCalls if both exist)
                if num_basecalls > 0:
                    logger.debug(
                        "Found %s with %d files (BaseCalls)",
                        run_path.name,
                        num_basecalls,
                    )
                    unstable_runs.append(
                        RunInfo(
                            run_dir=run_path,
                            file_count=num_basecalls,
                            glob_pattern=basecalls_pattern,
                        )
                    )
                elif num_alignment > 0:
                    logger.debug(
                        "Found %s with %d files (Alignment/Fastq)",
                        run_path.name,
                        num_alignment,
                    )
                    unstable_runs.append(
                        RunInfo(
                            run_dir=run_path,
                            file_count=num_alignment,
                            glob_pattern=alignment_pattern,
                        )
                    )
                else:
                    logger.debug("Skipping %s (no FASTQ files found)", run_path.name)

            except (OSError, IOError) as e:
                # Skip individual runs that fail to read
                logger.info("Error scanning run %s: %s", run_path.name, e)
                continue

    except (OSError, IOError) as e:
        # If we can't even list the directory, log and return empty
        logger.info("Error scanning runs directory %s: %s", runs_dir, e)
        return []

    return unstable_runs


def monitor_run_completion(runs_dir: Path) -> None:
    """
    Monitor runs for completion and mark them when stable.

    This function runs in an infinite loop, checking for new runs and monitoring
    existing runs for file stability. When a run's file count stops changing,
    it's considered complete and a 'needsprocessing' marker is created.

    This is designed to run in a daemon thread with a fixed 60-second polling interval.
    It includes crash protection - if any unexpected error occurs, it will log the
    error and restart after a short delay.

    Args:
        runs_dir: Directory containing MiSeq run folders (e.g., /data/RAW_DATA/MiSeq/runs)
    """

    # Outer loop for crash recovery
    while True:
        try:
            _monitor_run_completion_inner(runs_dir)
        except KeyboardInterrupt:
            logger.info("Run completion monitor shutting down")
            raise
        except Exception:
            logger.warning(
                "Run completion monitor crashed unexpectedly, restarting in %d seconds",
                CRASH_RECOVERY_DELAY,
                exc_info=True,
            )
            sleep(CRASH_RECOVERY_DELAY)


def check_run_completions(runs_dir: Path, monitoring: List[RunInfo]) -> None:
    # Find new unstable runs
    unstable_runs = find_unstable_runs(runs_dir)

    # Add new runs to monitoring list
    for run_info in unstable_runs:
        if not any(r.run_dir == run_info.run_dir for r in monitoring):
            monitoring.append(run_info)
            logger.debug("Now monitoring run: %s", run_info.run_dir.name)

    if monitoring:
        logger.debug("Currently monitoring %d run(s)", len(monitoring))
    else:
        logger.debug("No runs currently being monitored")

    # Sleep before checking stability
    sleep(COMPLETION_POLL_INTERVAL)

    # Check each monitored run for stability
    completed_indices = []
    for i, run_info in enumerate(monitoring):
        try:
            # Re-glob to get current file count
            current_files = list(run_info.run_dir.glob(run_info.glob_pattern))
            current_count = len(current_files)

            if current_count == run_info.file_count:
                # File count hasn't changed - run is stable
                marker_path = run_info.run_dir / "needsprocessing"
                disk_operations.touch(marker_path)
                logger.info("Marked run as ready: %s", run_info.run_dir.name)
                completed_indices.append(i)
            else:
                # File count changed - still being written
                logger.debug(
                    "Run %s still growing: %d -> %d files",
                    run_info.run_dir.name,
                    run_info.file_count,
                    current_count,
                )
                # Update the file count for next check
                monitoring[i] = RunInfo(
                    run_dir=run_info.run_dir,
                    file_count=current_count,
                    glob_pattern=run_info.glob_pattern,
                )
        except (OSError, IOError) as e:
            # If we can't check a run, log and continue to next
            logger.info("Error checking run %s: %s", run_info.run_dir.name, e)
            continue

    # Remove completed runs from monitoring list (reverse order to preserve indices)
    for i in reversed(completed_indices):
        del monitoring[i]


def _monitor_run_completion_inner(runs_dir: Path) -> None:
    """
    Inner monitoring loop - actual monitoring logic.

    Separated from monitor_run_completion to allow crash recovery.
    """
    # Track runs we're currently monitoring
    monitoring: List[RunInfo] = []

    while True:
        check_run_completions(runs_dir, monitoring)


def main(argv: Sequence[str]) -> int:
    parser = argparse.ArgumentParser(
        description="Monitor MiSeq run directories for completion."
    )
    parser.add_argument("--raw-data", type=Path, help="Directory containing MiSeq run folders")
    parser.add_argument("--once", action="store_true", help="Run once and exit")
    args = parser.parse_args(argv)
    raw_data = args.raw_data
    if raw_data is None:
        raw_data = os.getenv("MISEQ_RAW_DATA")
        if raw_data is None:
            logger.error("Error: --runs_dir argument or MISEQ_RAW_DATA environment variable is required")
            return 1

    run_dir = Path(raw_data) / "MiSeq" / "runs"
    if args.once:
        check_run_completions(run_dir, [])
    else:
        monitor_run_completion(run_dir)
    return 0


def entry() -> None:
    sys.exit(main(sys.argv[1:]))


if __name__ == "__main__": entry() # noqa
