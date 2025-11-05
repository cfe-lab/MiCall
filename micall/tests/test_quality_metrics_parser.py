from io import BytesIO
from struct import pack
from unittest import TestCase
from miseqinteropreader.read_records import read_quality
from miseqinteropreader.models import QualityRecord


def quality_record_to_dict(record: QualityRecord) -> dict:
    """Convert QualityRecord to dict with quality_bins tuple"""
    return {
        'lane': record.lane,
        'tile': record.tile,
        'cycle': record.cycle,
        'quality_bins': tuple(record.quality_bins)
    }


def summarize_quality_records(records, summary, read_lengths=None):
    """Wrapper to match old API - modifies summary dict in place

    Note: This is a standalone implementation that doesn't require InterOpReader
    since the MiCall tests pass plain dicts, not QualityRecord objects.
    """
    good_count = total_count = 0
    good_reverse = total_reverse = 0
    if read_lengths is None:
        last_forward_cycle = first_reverse_cycle = None
    else:
        last_forward_cycle = read_lengths[0]
        first_reverse_cycle = sum(read_lengths[:-1]) + 1
    for record in records:
        cycle = record['cycle']
        cycle_clusters = sum(record['quality_bins'])
        cycle_good = sum(record['quality_bins'][29:])

        if last_forward_cycle is None or cycle <= last_forward_cycle:
            total_count += cycle_clusters
            good_count += cycle_good
        elif cycle >= first_reverse_cycle:
            total_reverse += cycle_clusters
            good_reverse += cycle_good

    if total_count > 0:
        summary['q30_fwd'] = good_count/float(total_count)
    if total_reverse > 0:
        summary['q30_rev'] = good_reverse/float(total_reverse)


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
        self.sample_stream = BytesIO(pack(format_string, *self.sample_data))

    def test_load(self):
        expected_records = [dict(lane=1,
                                 tile=2,
                                 cycle=3,
                                 quality_bins=tuple(range(101, 151)))]

        records = [quality_record_to_dict(r) for r in read_quality(self.sample_stream)]

        self.assertEqual(expected_records, records)

    def test_new_version(self):
        self.sample_data[:2] = [5, 207]
        self.sample_data.append(42)
        self.sample_data.extend(self.sample_data[2:])
        format_string = '<BB' + 2*('HHH' + 50*'L' + 'B')
        self.sample_stream = BytesIO(pack(format_string, *self.sample_data))
        expected_records = [dict(lane=1,
                                 tile=2,
                                 cycle=3,
                                 quality_bins=tuple(range(101, 151)))] * 2

        records = [quality_record_to_dict(r) for r in read_quality(self.sample_stream)]

        self.maxDiff = 1000
        self.assertEqual(expected_records, records)

    def test_summarize(self):
        records = [dict(cycle=1, quality_bins=[0] * 50),
                   dict(cycle=101, quality_bins=[0] * 50)]
        records[0]['quality_bins'][28] = 1
        records[0]['quality_bins'][29] = 1
        records[1]['quality_bins'][29] = 1
        expected_summary = dict(q30_fwd=2/3.0)

        summary = {}
        summarize_quality_records(records, summary)

        self.assertEqual(expected_summary, summary)

    def test_summarize_reverse(self):
        read_lengths = [95, 5, 95]
        records = [dict(cycle=1, quality_bins=[0] * 50),
                   dict(cycle=101, quality_bins=[0] * 50)]
        records[0]['quality_bins'][28] = 1
        records[0]['quality_bins'][29] = 1
        records[1]['quality_bins'][29] = 1
        expected_summary = dict(q30_fwd=0.5, q30_rev=1.0)

        summary = {}
        summarize_quality_records(records, summary, read_lengths)

        self.assertEqual(expected_summary, summary)

    def test_summarize_ignores_index_reads(self):
        read_lengths = [95, 5, 95]
        records = [dict(cycle=96, quality_bins=[0] * 50),
                   dict(cycle=100, quality_bins=[0] * 50)]
        records[0]['quality_bins'][29] = 1
        records[1]['quality_bins'][29] = 1
        expected_summary = {}

        summary = {}
        summarize_quality_records(records, summary, read_lengths)

        self.assertEqual(expected_summary, summary)

    def test_summarize_blank(self):
        records = []
        expected_summary = {}

        summary = {}
        summarize_quality_records(records, summary)

        self.assertEqual(expected_summary, summary)
