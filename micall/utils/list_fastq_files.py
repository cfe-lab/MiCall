"""Utility functions for finding FASTQ files in sequencing run folders."""

import re
from pathlib import Path
from typing import List, Union, Optional


def _get_base_calls_path(run_path: Union[str, Path]) -> Path:
    """Get the standard BaseCalls directory path for a MiSeq run.

    Args:
        run_path: Path to the sequencing run folder

    Returns:
        Path to the BaseCalls folder (Data/Intensities/BaseCalls)
    """
    if isinstance(run_path, str):
        run_path = Path(run_path)
    return run_path / "Data" / "Intensities" / "BaseCalls"


def find_fastq_source_folder(
    run_path: Union[str, Path], pattern: str = "*_R1_*"
) -> Optional[Path]:
    """Find the folder containing FASTQ files in a sequencing run.

    First tries the standard MiSeq structure (Data/Intensities/BaseCalls),
    then tries Alignment_\\d+/.*/Fastq/ directories (preferring highest alignment
    number, then lexicographically largest subdirectory name),
    then falls back to the run_path directly.

    Args:
        run_path: Path to the sequencing run folder
        pattern: Glob pattern for FASTQ files (default: "*_R1_*")

    Returns:
        Path to the folder containing FASTQ files, or None if no files found
    """
    if isinstance(run_path, str):
        run_path = Path(run_path)

    # Try standard MiSeq structure first
    base_calls_path = _get_base_calls_path(run_path)
    if base_calls_path.exists():
        fastq_files = list(base_calls_path.glob(pattern))
        if fastq_files:
            return base_calls_path

    # Try Alignment_\d+/.*/Fastq/ directories
    # Collect all matching alignment directories and sort by number (highest first)
    alignment_dirs = []
    for alignment_dir in run_path.glob("Alignment_*"):
        if not alignment_dir.is_dir():
            continue
        # Check if the directory name matches Alignment_\d+ pattern
        match = re.match(r"Alignment_(\d+)$", alignment_dir.name)
        if match:
            alignment_num = int(match.group(1))
            alignment_dirs.append((alignment_num, alignment_dir))

    # Sort by alignment number (highest first)
    alignment_dirs.sort(key=lambda x: x[0], reverse=True)

    for alignment_num, alignment_dir in alignment_dirs:
        # Search for Fastq directories within this alignment directory
        # Collect all matching Fastq directories and sort lexicographically (largest first)
        fastq_dirs = []
        for fastq_dir in alignment_dir.glob("*/Fastq"):
            if fastq_dir.is_dir():
                # Get the parent directory name (the X in Alignment_\d+/X/Fastq)
                parent_name = fastq_dir.parent.name
                fastq_dirs.append((parent_name, fastq_dir))

        # Sort by parent name lexicographically (largest first)
        fastq_dirs.sort(key=lambda x: x[0], reverse=True)

        for parent_name, fastq_dir in fastq_dirs:
            fastq_files = list(fastq_dir.glob(pattern))
            if fastq_files:
                return fastq_dir

    # Fall back to run_path directly
    fastq_files = list(run_path.glob(pattern))
    if fastq_files:
        return run_path

    return None


def list_fastq_files(
    run_path: Union[str, Path],
    pattern: str = "*_R1_*.fastq*",
    fallback_to_run_path: bool = True,
) -> List[Path]:
    """List FASTQ files in a sequencing run folder.

    First tries the standard MiSeq structure (Data/Intensities/BaseCalls),
    then tries Alignment_\\d+/.*/Fastq/ directories (preferring highest number),
    then falls back to searching the run_path directly if requested.

    Args:
        run_path: Path to the sequencing run folder
        pattern: Glob pattern for FASTQ files (default: "*_R1_*.fastq*")
        fallback_to_run_path: If True, search run_path directly if BaseCalls
                             doesn't exist or has no matching files

    Returns:
        List of Path objects for matching FASTQ files
    """
    if isinstance(run_path, str):
        run_path = Path(run_path)

    # Use find_fastq_source_folder to locate the source directory
    source_folder = find_fastq_source_folder(run_path, pattern)
    if source_folder is not None:
        return list(source_folder.glob(pattern))

    # If no source folder found and fallback is disabled, return empty list
    if not fallback_to_run_path:
        return []

    # Final fallback: search run_path directly
    return list(run_path.glob(pattern))


def list_fastq_file_names(
    run_path: Union[str, Path],
    pattern: str = "*_R1_*.fastq*",
    fallback_to_run_path: bool = True,
) -> List[str]:
    """List FASTQ file names (without directory path) in a sequencing run folder.

    First tries the standard MiSeq structure (Data/Intensities/BaseCalls),
    then tries Alignment_\\d+/.*/Fastq/ directories (preferring highest number),
    then falls back to searching the run_path directly if requested.

    Args:
        run_path: Path to the sequencing run folder
        pattern: Glob pattern for FASTQ files (default: "*_R1_*.fastq*")
        fallback_to_run_path: If True, search run_path directly if BaseCalls
                             doesn't exist or has no matching files

    Returns:
        List of file names (strings) for matching FASTQ files
    """
    fastq_files = list_fastq_files(run_path, pattern, fallback_to_run_path)
    return [f.name for f in fastq_files]
