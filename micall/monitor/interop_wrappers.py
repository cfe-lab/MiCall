

from miseqinteropreader.interop_reader import InterOpReader
from miseqinteropreader.models import TileMetricCodes


def summarize_quality(quality_metrics_path, summary, read_lengths):
    """Summarize quality metrics from InterOp files.

    Modifies summary dict in place with q30_fwd and q30_rev values.
    """
    # Extract the run folder path (parent of InterOp folder)
    import os
    interop_path = os.path.dirname(quality_metrics_path)
    run_path = os.path.dirname(interop_path)

    reader = InterOpReader(run_path)
    quality_records = reader.read_quality_records()

    good_count = total_count = 0
    good_reverse = total_reverse = 0
    last_forward_cycle = read_lengths[0]
    first_reverse_cycle = sum(read_lengths[:-1]) + 1

    for record in quality_records:
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


def summarize_tiles(tile_metrics_path, summary):
    """Summarize tile metrics from InterOp files.

    Modifies summary dict in place with cluster_density and pass_rate values.
    """
    # Extract the run folder path (parent of InterOp folder)
    import os
    interop_path = os.path.dirname(tile_metrics_path)
    run_path = os.path.dirname(interop_path)

    reader = InterOpReader(run_path)
    tile_records = reader.read_tile_records()

    density_sum = 0.0
    density_count = 0
    total_clusters = 0.0
    passing_clusters = 0.0

    for record in tile_records:
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
