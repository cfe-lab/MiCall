from io import BytesIO
from struct import pack
from unittest import TestCase
from miseqinteropreader.read_records import read_quality
from miseqinteropreader.models import QualityRecord
from miseqinteropreader import ReadLengths4
from micall.utils.interop_wrappers import summarize_quality_records


def quality_record_to_dict(record: QualityRecord) -> dict:
    """Convert QualityRecord to dict with quality_bins tuple"""
    return {
        'lane': record.lane,
        'tile': record.tile,
        'cycle': record.cycle,
        'quality_bins': tuple(record.quality_bins)
    }


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
        records = [QualityRecord(lane=1, tile=1, cycle=1, quality_bins=[0] * 50),
                   QualityRecord(lane=1, tile=1, cycle=101, quality_bins=[0] * 50)]
        records[0].quality_bins[28] = 1
        records[0].quality_bins[29] = 1
        records[1].quality_bins[29] = 1
        expected_summary = dict(q30_fwd=2/3.0)

        summary = {}
        # Use ReadLengths4 for single combined read
        read_lengths = ReadLengths4(forward_read=200, index1=0, index2=0, reverse_read=0)
        summarize_quality_records(records, summary, read_lengths)

        self.assertEqual(expected_summary, summary)

    def test_summarize_reverse(self):
        # Use ReadLengths4 format
        read_lengths = ReadLengths4(forward_read=95, index1=3, index2=2, reverse_read=95)
        records = [QualityRecord(lane=1, tile=1, cycle=1, quality_bins=[0] * 50),
                   QualityRecord(lane=1, tile=1, cycle=101, quality_bins=[0] * 50)]
        records[0].quality_bins[28] = 1
        records[0].quality_bins[29] = 1
        records[1].quality_bins[29] = 1
        expected_summary = dict(q30_fwd=0.5, q30_rev=1.0)

        summary = {}
        summarize_quality_records(records, summary, read_lengths)

        self.assertEqual(expected_summary, summary)

    def test_summarize_ignores_index_reads(self):
        # Use ReadLengths4 format
        read_lengths = ReadLengths4(forward_read=95, index1=3, index2=2, reverse_read=95)
        records = [QualityRecord(lane=1, tile=1, cycle=96, quality_bins=[0] * 50),
                   QualityRecord(lane=1, tile=1, cycle=100, quality_bins=[0] * 50)]
        records[0].quality_bins[29] = 1
        records[1].quality_bins[29] = 1
        expected_summary = {}

        summary = {}
        summarize_quality_records(records, summary, read_lengths)

        self.assertEqual(expected_summary, summary)

    def test_summarize_blank(self):
        records = []
        expected_summary = {}

        summary = {}
        # Use ReadLengths4 for single combined read
        read_lengths = ReadLengths4(forward_read=200, index1=0, index2=0, reverse_read=0)
        summarize_quality_records(records, summary, read_lengths)

        self.assertEqual(expected_summary, summary)

    def test_summarize_with_readlengths4(self):
        # Test using ReadLengths4 format
        read_lengths = ReadLengths4(
            forward_read=95,
            index1=3,
            index2=2,
            reverse_read=95
        )
        records = [QualityRecord(lane=1, tile=1, cycle=1, quality_bins=[0] * 50),
                   QualityRecord(lane=1, tile=1, cycle=101, quality_bins=[0] * 50)]
        records[0].quality_bins[28] = 1
        records[0].quality_bins[29] = 1
        records[1].quality_bins[29] = 1
        expected_summary = dict(q30_fwd=0.5, q30_rev=1.0)

        summary = {}
        summarize_quality_records(records, summary, read_lengths)

        self.assertEqual(expected_summary, summary)
