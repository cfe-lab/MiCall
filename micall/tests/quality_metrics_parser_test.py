from cStringIO import StringIO
from struct import pack
from unittest import TestCase
from micall.monitor.quality_metrics_parser import read_quality,\
    summarize_quality


class QualityMetricsParserTest(TestCase):
    def setUp(self):
        self.addTypeEqualityFunc(str, self.assertMultiLineEqual)
        self.sample_data = [4,      # file version
                            206,    # record size
                            1,      # lane
                            2,      # tile
                            3]      # cycle
        self.sample_data.extend(range(101, 151))
        format_string = '<BBHHH' + 'L'*50
        self.sample_stream = StringIO(pack(format_string, *self.sample_data))

    def test_load(self):
        expected_records = [dict(lane=1,
                                 tile=2,
                                 cycle=3,
                                 quality_bins=tuple(range(101, 151)))]

        records = list(read_quality(self.sample_stream))

        self.assertEqual(expected_records, records)

    def test_old_version(self):
        self.sample_data[0] = 3
        format_string = '<BBHHH' + 'L'*50
        self.sample_stream = StringIO(pack(format_string, *self.sample_data))

        try:
            list(read_quality(self.sample_stream))

            self.fail('Should have thrown.')
        except IOError as ex:
            self.assertEqual('Old quality metrics file version: 3', str(ex))

    def test_new_version(self):
        self.sample_data[:2] = [5, 207]
        self.sample_data.append(42)
        self.sample_data.extend(self.sample_data[2:])
        format_string = '<BB' + 2*('HHH' + 50*'L' + 'B')
        self.sample_stream = StringIO(pack(format_string, *self.sample_data))
        expected_records = [dict(lane=1,
                                 tile=2,
                                 cycle=3,
                                 quality_bins=tuple(range(101, 151)))] * 2

        records = list(read_quality(self.sample_stream))

        self.maxDiff = 1000
        self.assertEqual(expected_records, records)

    def test_partial_record(self):
        self.sample_data = self.sample_data[:3]
        format_string = '<BBH'
        self.sample_stream = StringIO(pack(format_string, *self.sample_data))

        try:
            list(read_quality(self.sample_stream))

            self.fail('Should have thrown.')
        except IOError as ex:
            self.assertEqual('Partial record of length 2 found in quality metrics file.',
                             str(ex))

    def test_summarize(self):
        records = [dict(cycle=1, quality_bins=[0] * 50),
                   dict(cycle=101, quality_bins=[0] * 50)]
        records[0]['quality_bins'][28] = 1
        records[0]['quality_bins'][29] = 1
        records[1]['quality_bins'][29] = 1
        expected_portion = 2/3.0

        portion = summarize_quality(records)

        self.assertAlmostEqual(expected_portion, portion)

    def test_summarize_reverse(self):
        read_lengths = [95, 5, 95]
        records = [dict(cycle=1, quality_bins=[0] * 50),
                   dict(cycle=101, quality_bins=[0] * 50)]
        records[0]['quality_bins'][28] = 1
        records[0]['quality_bins'][29] = 1
        records[1]['quality_bins'][29] = 1
        expected_portions = (0.5, 1.0)

        portions = summarize_quality(records, read_lengths)

        self.assertEqual(expected_portions, portions)

    def test_summarize_ignores_index(self):
        read_lengths = [95, 5, 95]
        records = [dict(cycle=96, quality_bins=[0] * 50),
                   dict(cycle=100, quality_bins=[0] * 50)]
        records[0]['quality_bins'][29] = 1
        records[1]['quality_bins'][29] = 1
        expected_portions = (0, 0)

        portions = summarize_quality(records, read_lengths)

        self.assertEqual(expected_portions, portions)

    def test_summarize_blank(self):
        records = []
        expected_portion = 0

        portion = summarize_quality(records)

        self.assertEqual(expected_portion, portion)
