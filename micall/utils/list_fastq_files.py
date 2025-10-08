"""Utility functions for finding FASTQ files in sequencing run folders."""

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
    run_path: Union[str, Path],
    pattern: str = "*_R1_*"
) -> Optional[Path]:
    """Find the folder containing FASTQ files in a sequencing run.

    First tries the standard MiSeq structure (Data/Intensities/BaseCalls),
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

    # Fall back to run_path directly
    fastq_files = list(run_path.glob(pattern))
    if fastq_files:
        return run_path

    return None


def list_fastq_files(
    run_path: Union[str, Path],
    pattern: str = "*_R1_*.fastq*",
    fallback_to_run_path: bool = True
) -> List[Path]:
    """List FASTQ files in a sequencing run folder.

    First tries the standard MiSeq structure (Data/Intensities/BaseCalls),
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

    # Try standard MiSeq structure first
    base_calls_path = find_fastq_source_folder(run_path)
    if base_calls_path is not None:
        fastq_files = list(base_calls_path.glob(pattern))
        if fastq_files:
            return fastq_files

    # Fall back to run_path directly if requested
    if fallback_to_run_path:
        fastq_files = list(run_path.glob(pattern))
        if fastq_files:
            return fastq_files

    return []


def list_fastq_file_names(
    run_path: Union[str, Path],
    pattern: str = "*_R1_*.fastq*",
    fallback_to_run_path: bool = True
) -> List[str]:
    """List FASTQ file names (without directory path) in a sequencing run folder.

    First tries the standard MiSeq structure (Data/Intensities/BaseCalls),
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
