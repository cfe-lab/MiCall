

from pathlib import Path
from typing import Sequence, TextIO

from miseqinteropreader.interop_reader import InterOpReader
from miseqinteropreader.models import (
    TileMetricCodes,
    QualityRecord,
    TileMetricRecord,
    ReadLengths4,
)
from miseqinteropreader.error_metrics_parser import write_phix_csv


def summarize_quality_records(
    records: Sequence[QualityRecord],
    summary: dict,
    read_lengths: ReadLengths4
) -> None:
    """Summarize quality records into a dictionary.

    Args:
        records: Sequence of QualityRecord objects
        summary: Dictionary to populate with q30_fwd and q30_rev values
        read_lengths: ReadLengths4 specifying read structure
    """

    good_count = total_count = 0
    good_reverse = total_reverse = 0
    last_forward_cycle = read_lengths.forward_read
    first_reverse_cycle = read_lengths.forward_read + read_lengths.index1 + read_lengths.index2 + 1

    for record in records:
        cycle = record.cycle
        cycle_clusters = sum(record.quality_bins)
        cycle_good = sum(record.quality_bins[29:])

        if cycle <= last_forward_cycle:
            total_count += cycle_clusters
            good_count += cycle_good
        elif cycle >= first_reverse_cycle:
            total_reverse += cycle_clusters
            good_reverse += cycle_good

    if total_count > 0:
        summary['q30_fwd'] = good_count / float(total_count)
    if total_reverse > 0:
        summary['q30_rev'] = good_reverse / float(total_reverse)


def summarize_quality(
    quality_metrics_path: Path,
    summary: dict,
    read_lengths: ReadLengths4
) -> None:
    """Summarize quality metrics from InterOp files.

    Modifies summary dict in place with q30_fwd and q30_rev values.

    Args:
        quality_metrics_path: Path to quality metrics file
        summary: Dictionary to populate
        read_lengths: ReadLengths4 specifying read structure
    """
    # Extract the run folder path (parent of InterOp folder)
    interop_path = quality_metrics_path.parent
    run_path = interop_path.parent

    reader = InterOpReader(run_path)
    quality_records = reader.read_quality_records()
    summarize_quality_records(quality_records, summary, read_lengths)


def summarize_tile_records(records: Sequence[TileMetricRecord], summary: dict) -> None:
    density_sum = 0.0
    density_count = 0
    total_clusters = 0.0
    passing_clusters = 0.0

    for record in records:
        if record.metric_code == TileMetricCodes.CLUSTER_DENSITY:
            density_sum += record.metric_value
            density_count += 1
        elif record.metric_code == TileMetricCodes.CLUSTER_COUNT:
            total_clusters += record.metric_value
        elif record.metric_code == TileMetricCodes.CLUSTER_COUNT_PASSING_FILTERS:
            passing_clusters += record.metric_value

    if density_count > 0:
        summary['cluster_density'] = density_sum / density_count
    if total_clusters > 0:
        summary['pass_rate'] = passing_clusters / total_clusters


def summarize_tiles(tile_metrics_path: Path, summary: dict) -> None:
    """Summarize tile metrics from InterOp files.

    Modifies summary dict in place with cluster_density and pass_rate values.
    """
    # Extract the run folder path (parent of InterOp folder)
    interop_path = tile_metrics_path.parent
    run_path = interop_path.parent

    reader = InterOpReader(run_path)
    tile_records = reader.read_tile_records()
    summarize_tile_records(tile_records, summary)


def summarize_error_metrics(
    error_metrics_path: Path,
    quality_csv_file: TextIO,
    summary: dict,
    read_lengths: ReadLengths4
) -> None:
    """Summarize error metrics from InterOp files and write quality CSV.

    Modifies summary dict in place with error_rate_fwd and error_rate_rev values.
    Also writes the quality data to the provided file handle.

    Args:
        error_metrics_path: Path to error metrics file (ErrorMetricsOut.bin)
        quality_csv_file: Open file handle to write quality CSV data
        summary: Dictionary to populate with error rate values
        read_lengths: ReadLengths4 specifying read structure
    """
    # Extract the run folder path (parent of InterOp folder)
    interop_path = error_metrics_path.parent
    run_path = interop_path.parent

    reader = InterOpReader(run_path)
    records = reader.read_error_records()
    phix_summary = write_phix_csv(quality_csv_file, records, read_lengths)

    # Store error rates in the summary dictionary
    if hasattr(phix_summary, 'error_rate_forward'):
        summary['error_rate_fwd'] = phix_summary.error_rate_forward
    if hasattr(phix_summary, 'error_rate_reverse'):
        summary['error_rate_rev'] = phix_summary.error_rate_reverse

