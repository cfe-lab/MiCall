from io import BytesIO
from struct import pack
from unittest import TestCase

from miseqinteropreader.read_records import read_tiles
from miseqinteropreader.models import TileMetricCodes


# Alias for compatibility with old tests
MetricCodes = TileMetricCodes


def summarize_tile_records(records, summary):
    """Wrapper to match old API - modifies summary dict in place

    Note: This is a standalone implementation that doesn't require InterOpReader
    since the MiCall tests pass plain dicts, not TileMetricRecord objects.
    """
    density_sum = 0.0
    density_count = 0
    total_clusters = 0.0
    passing_clusters = 0.0

    for record in records:
        if record['metric_code'] == TileMetricCodes.CLUSTER_DENSITY:
            density_sum += record['metric_value']
            density_count += 1
        elif record['metric_code'] == TileMetricCodes.CLUSTER_COUNT:
            total_clusters += record['metric_value']
        elif record['metric_code'] == TileMetricCodes.CLUSTER_COUNT_PASSING_FILTERS:
            passing_clusters += record['metric_value']

    if density_count > 0:
        summary['cluster_density'] = density_sum / density_count
    if total_clusters > 0:
        summary['pass_rate'] = passing_clusters / total_clusters


class TileMetricsParserTest(TestCase):
    def setUp(self):
        self.addTypeEqualityFunc(str, self.assertMultiLineEqual)
        self.sample_data = [4,      # file version
                            10,     # record size
                            1,      # lane
                            2,      # tile
                            100,    # metric code
                            4.0]    # metric value
        format_string = '<BBHHHf'
        self.sample_stream = BytesIO(pack(format_string, *self.sample_data))

    def test_load(self):
        expected_records = [dict(lane=1,
                                 tile=2,
                                 metric_code=100,
                                 metric_value=4.0)]

        records = [r.model_dump() for r in read_tiles(self.sample_stream)]

        self.assertEqual(expected_records, records)

    def test_summarize(self):
        records = [dict(metric_code=MetricCodes.CLUSTER_DENSITY,
                        metric_value=3.0)]
        expected_summary = dict(cluster_density=3.0)

        summary = {}
        summarize_tile_records(records, summary)

        self.assertEqual(expected_summary, summary)

    def test_summarize_empty(self):
        records = []
        expected_summary = {}

        summary = {}
        summarize_tile_records(records, summary)

        self.assertEqual(expected_summary, summary)

    def test_summarize_multiple(self):
        records = [dict(metric_code=MetricCodes.CLUSTER_DENSITY,
                        metric_value=3.0),
                   dict(metric_code=MetricCodes.CLUSTER_DENSITY,
                        metric_value=4.0)]
        expected_summary = dict(cluster_density=3.5)

        summary = {}
        summarize_tile_records(records, summary)

        self.assertEqual(expected_summary, summary)

    def test_summarize_skips_other_codes(self):
        records = [dict(metric_code=MetricCodes.CLUSTER_DENSITY,
                        metric_value=3.0),
                   dict(metric_code=MetricCodes.CLUSTER_DENSITY_PASSING_FILTERS,
                        metric_value=4.0)]
        expected_summary = dict(cluster_density=3.0)

        summary = {}
        summarize_tile_records(records, summary)

        self.assertEqual(expected_summary, summary)

    def test_summarize_passing_clusters(self):
        records = [dict(metric_code=MetricCodes.CLUSTER_COUNT,
                        metric_value=100.0),
                   dict(metric_code=MetricCodes.CLUSTER_COUNT_PASSING_FILTERS,
                        metric_value=50.0)]
        expected_summary = dict(pass_rate=0.5)

        summary = {}
        summarize_tile_records(records, summary)

        self.assertEqual(expected_summary, summary)
