from io import BytesIO
from struct import pack
from unittest import TestCase

from miseqinteropreader.read_records import read_tiles
from miseqinteropreader.models import TileMetricCodes, TileMetricRecord
from micall.utils.interop_wrappers import summarize_tile_records


# Alias for compatibility with old tests
MetricCodes = TileMetricCodes


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
        records = [TileMetricRecord(lane=1, tile=1,
                                     metric_code=MetricCodes.CLUSTER_DENSITY,
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
        records = [TileMetricRecord(lane=1, tile=1,
                                     metric_code=MetricCodes.CLUSTER_DENSITY,
                                     metric_value=3.0),
                   TileMetricRecord(lane=1, tile=2,
                                    metric_code=MetricCodes.CLUSTER_DENSITY,
                                    metric_value=4.0)]
        expected_summary = dict(cluster_density=3.5)

        summary = {}
        summarize_tile_records(records, summary)

        self.assertEqual(expected_summary, summary)

    def test_summarize_skips_other_codes(self):
        records = [TileMetricRecord(lane=1, tile=1,
                                     metric_code=MetricCodes.CLUSTER_DENSITY,
                                     metric_value=3.0),
                   TileMetricRecord(lane=1, tile=1,
                                    metric_code=MetricCodes.CLUSTER_DENSITY_PASSING_FILTERS,
                                    metric_value=4.0)]
        expected_summary = dict(cluster_density=3.0)

        summary = {}
        summarize_tile_records(records, summary)

        self.assertEqual(expected_summary, summary)

    def test_summarize_passing_clusters(self):
        records = [TileMetricRecord(lane=1, tile=1,
                                     metric_code=MetricCodes.CLUSTER_COUNT,
                                     metric_value=100.0),
                   TileMetricRecord(lane=1, tile=1,
                                    metric_code=MetricCodes.CLUSTER_COUNT_PASSING_FILTERS,
                                    metric_value=50.0)]
        expected_summary = dict(pass_rate=0.5)

        summary = {}
        summarize_tile_records(records, summary)

        self.assertEqual(expected_summary, summary)
