from io import BytesIO, StringIO
from struct import pack
from unittest import TestCase

from miseqinteropreader.error_metrics_parser import write_phix_csv
from miseqinteropreader.models import ErrorRecord, ReadLengths4
from miseqinteropreader.read_records import read_errors, read_records


class RecordsParserTest(TestCase):
    def pack_data(self):
        format_string = '<bbcccc'
        self.sample_stream = BytesIO(pack(format_string, *self.sample_data))
        self.sample_stream.name = 'test_file'

    def setUp(self):
        self.addTypeEqualityFunc(str, self.assertMultiLineEqual)
        self.sample_data = [1,      # file version
                            4,      # record size
                            b'A',    # field 1
                            b'B',    # field 2
                            b'C',    # field 3
                            b'D']    # field 4
        self.pack_data()

    def test_load(self):
        expected_records = [b'ABCD']

        records = [data for data, _ in read_records(self.sample_stream, min_version=1)]

        self.assertEqual(expected_records, records)

    def test_load_multiple_records(self):
        self.sample_data[1] = 2  # record size
        self.pack_data()
        expected_records = [b'AB', b'CD']

        records = [data for data, _ in read_records(self.sample_stream, min_version=1)]

        self.assertEqual(expected_records, records)

    def test_old_version(self):

        records = read_records(self.sample_stream, min_version=3)

        self.assertRaisesRegex(
            ValueError,
            'File version 1 is less than minimum version 3 in test_file.',
            records.__next__)

    def test_partial_record(self):
        self.sample_data[1] = 3
        self.pack_data()
        records = read_records(self.sample_stream, min_version=1)
        data, _ = next(records)

        self.assertEqual(b'ABC', data)
        self.assertRaisesRegex(
            RuntimeError,
            'Partial record of length 1 found in test_file.',
            records.__next__)


class ErrorMetricsParserTest(TestCase):
    def setUp(self):
        self.addTypeEqualityFunc(str, self.assertMultiLineEqual)
        self.sample_data = [3,      # file version
                            30,     # record size
                            1,      # lane
                            2,      # tile
                            3,      # cycle
                            0.5,    # error rate
                            4,      # num reads with 0 errors
                            5,      # num reads with 1 error
                            6,      # num reads with 2 errors
                            7,      # num reads with 3 errors
                            8]      # num reads with 4 errors
        format_string = '<bbHHHfLLLLL'
        self.sample_stream = BytesIO(pack(format_string, *self.sample_data))

    def test_load(self):
        expected_records = [dict(lane=1,
                                 tile=2,
                                 cycle=3,
                                 error_rate=0.5,
                                 num_0_errors=4,
                                 num_1_errors=5,
                                 num_2_errors=6,
                                 num_3_errors=7,
                                 num_4_errors=8)]

        records = [r.model_dump() for r in read_errors(self.sample_stream)]

        self.assertEqual(expected_records, records)

    def test_new_version(self):
        self.sample_data[:2] = [4, 31]
        self.sample_data.append(42)
        self.sample_data.extend(self.sample_data[2:])
        format_string = '<bbHHHfLLLLLbHHHfLLLLLb'
        self.sample_stream = BytesIO(pack(format_string, *self.sample_data))
        expected_records = [dict(lane=1,
                                 tile=2,
                                 cycle=3,
                                 error_rate=0.5,
                                 num_0_errors=4,
                                 num_1_errors=5,
                                 num_2_errors=6,
                                 num_3_errors=7,
                                 num_4_errors=8)] * 2

        records = [r.model_dump() for r in read_errors(self.sample_stream)]

        self.maxDiff = 1000
        self.assertEqual(expected_records, records)

    def test_write(self):
        out_file = StringIO()
        records = [
            ErrorRecord(lane=1, tile=2, cycle=1, error_rate=0.1, num_0_errors=0, num_1_errors=0, num_2_errors=0, num_3_errors=0, num_4_errors=0),
            ErrorRecord(lane=1, tile=2, cycle=2, error_rate=0.2, num_0_errors=0, num_1_errors=0, num_2_errors=0, num_3_errors=0, num_4_errors=0)
        ]
        expected_csv = """\
tile,cycle,errorrate
2,1,0.1
2,2,0.2
"""

        write_phix_csv(out_file, records)

        self.assertEqual(expected_csv, out_file.getvalue())

    def test_write_rounded(self):
        out_file = StringIO()
        records = [
            ErrorRecord(lane=1, tile=2, cycle=1, error_rate=0.12345, num_0_errors=0, num_1_errors=0, num_2_errors=0, num_3_errors=0, num_4_errors=0),
            ErrorRecord(lane=1, tile=2, cycle=2, error_rate=0.23001, num_0_errors=0, num_1_errors=0, num_2_errors=0, num_3_errors=0, num_4_errors=0)
        ]
        expected_csv = """\
tile,cycle,errorrate
2,1,0.1235
2,2,0.23
"""

        write_phix_csv(out_file, records)

        self.assertEqual(expected_csv, out_file.getvalue())

    def test_write_sorted(self):
        out_file = StringIO()
        records = [
            ErrorRecord(lane=1, tile=2, cycle=2, error_rate=0.4, num_0_errors=0, num_1_errors=0, num_2_errors=0, num_3_errors=0, num_4_errors=0),
            ErrorRecord(lane=1, tile=2, cycle=1, error_rate=0.5, num_0_errors=0, num_1_errors=0, num_2_errors=0, num_3_errors=0, num_4_errors=0)
        ]
        expected_csv = """\
tile,cycle,errorrate
2,1,0.5
2,2,0.4
"""

        write_phix_csv(out_file, records)

        self.assertEqual(expected_csv, out_file.getvalue())

    def test_write_reverse(self):
        out_file = StringIO()
        records = [
            ErrorRecord(lane=1, tile=2, cycle=1, error_rate=0.1, num_0_errors=0, num_1_errors=0, num_2_errors=0, num_3_errors=0, num_4_errors=0),
            ErrorRecord(lane=1, tile=2, cycle=2, error_rate=0.2, num_0_errors=0, num_1_errors=0, num_2_errors=0, num_3_errors=0, num_4_errors=0),
            ErrorRecord(lane=1, tile=2, cycle=3, error_rate=0.3, num_0_errors=0, num_1_errors=0, num_2_errors=0, num_3_errors=0, num_4_errors=0),
            ErrorRecord(lane=1, tile=2, cycle=4, error_rate=0.4, num_0_errors=0, num_1_errors=0, num_2_errors=0, num_3_errors=0, num_4_errors=0)
        ]
        read_lengths = ReadLengths4(forward_read=2, index1=0, index2=0, reverse_read=2)
        expected_csv = """\
tile,cycle,errorrate
2,1,0.1
2,2,0.2
2,-1,0.3
2,-2,0.4
"""

        write_phix_csv(out_file, records, read_lengths)

        self.assertEqual(expected_csv, out_file.getvalue())

    def test_write_skip_indexes(self):
        out_file = StringIO()
        records = [
            ErrorRecord(lane=1, tile=2, cycle=1, error_rate=0.1, num_0_errors=0, num_1_errors=0, num_2_errors=0, num_3_errors=0, num_4_errors=0),
            ErrorRecord(lane=1, tile=2, cycle=2, error_rate=0.2, num_0_errors=0, num_1_errors=0, num_2_errors=0, num_3_errors=0, num_4_errors=0),
            ErrorRecord(lane=1, tile=2, cycle=3, error_rate=0.3, num_0_errors=0, num_1_errors=0, num_2_errors=0, num_3_errors=0, num_4_errors=0),
            ErrorRecord(lane=1, tile=2, cycle=4, error_rate=0.4, num_0_errors=0, num_1_errors=0, num_2_errors=0, num_3_errors=0, num_4_errors=0),
            ErrorRecord(lane=1, tile=2, cycle=5, error_rate=0.5, num_0_errors=0, num_1_errors=0, num_2_errors=0, num_3_errors=0, num_4_errors=0),
            ErrorRecord(lane=1, tile=2, cycle=6, error_rate=0.6, num_0_errors=0, num_1_errors=0, num_2_errors=0, num_3_errors=0, num_4_errors=0)
        ]
        read_lengths = ReadLengths4(forward_read=2, index1=1, index2=1, reverse_read=2)
        expected_csv = """\
tile,cycle,errorrate
2,1,0.1
2,2,0.2
2,-1,0.5
2,-2,0.6
"""

        write_phix_csv(out_file, records, read_lengths)

        self.assertEqual(expected_csv, out_file.getvalue())

    def test_write_missing(self):
        out_file = StringIO()
        records = [
            ErrorRecord(lane=1, tile=2, cycle=1, error_rate=0.1, num_0_errors=0, num_1_errors=0, num_2_errors=0, num_3_errors=0, num_4_errors=0),
            ErrorRecord(lane=1, tile=2, cycle=4, error_rate=0.4, num_0_errors=0, num_1_errors=0, num_2_errors=0, num_3_errors=0, num_4_errors=0)
        ]
        expected_csv = """\
tile,cycle,errorrate
2,1,0.1
2,2,
2,3,
2,4,0.4
"""

        write_phix_csv(out_file, records)

        self.assertEqual(expected_csv, out_file.getvalue())

    def test_write_missing_end(self):
        out_file = StringIO()
        records = [
            ErrorRecord(lane=1, tile=2, cycle=1, error_rate=0.1, num_0_errors=0, num_1_errors=0, num_2_errors=0, num_3_errors=0, num_4_errors=0),
            ErrorRecord(lane=1, tile=2, cycle=2, error_rate=0.2, num_0_errors=0, num_1_errors=0, num_2_errors=0, num_3_errors=0, num_4_errors=0),
            ErrorRecord(lane=1, tile=2, cycle=4, error_rate=0.4, num_0_errors=0, num_1_errors=0, num_2_errors=0, num_3_errors=0, num_4_errors=0),
            ErrorRecord(lane=1, tile=2, cycle=5, error_rate=0.5, num_0_errors=0, num_1_errors=0, num_2_errors=0, num_3_errors=0, num_4_errors=0)
        ]
        read_lengths = ReadLengths4(forward_read=3, index1=0, index2=0, reverse_read=3)
        expected_csv = """\
tile,cycle,errorrate
2,1,0.1
2,2,0.2
2,3,
2,-1,0.4
2,-2,0.5
2,-3,
"""

        write_phix_csv(out_file, records, read_lengths)

        self.assertEqual(expected_csv, out_file.getvalue())

    def test_summary(self):
        out_file = StringIO()
        records = [
            ErrorRecord(lane=1, tile=2, cycle=1, error_rate=0.25, num_0_errors=0, num_1_errors=0, num_2_errors=0, num_3_errors=0, num_4_errors=0),
            ErrorRecord(lane=1, tile=2, cycle=2, error_rate=0.75, num_0_errors=0, num_1_errors=0, num_2_errors=0, num_3_errors=0, num_4_errors=0)
        ]

        summary = write_phix_csv(out_file, records)

        self.assertEqual(0.5, summary.error_rate_forward)

    def test_summary_reverse(self):
        out_file = StringIO()
        records = [
            ErrorRecord(lane=1, tile=2, cycle=1, error_rate=0.75, num_0_errors=0, num_1_errors=0, num_2_errors=0, num_3_errors=0, num_4_errors=0),
            ErrorRecord(lane=1, tile=2, cycle=4, error_rate=0.375, num_0_errors=0, num_1_errors=0, num_2_errors=0, num_3_errors=0, num_4_errors=0),
            ErrorRecord(lane=1, tile=2, cycle=5, error_rate=0.125, num_0_errors=0, num_1_errors=0, num_2_errors=0, num_3_errors=0, num_4_errors=0)
        ]
        read_lengths = ReadLengths4(forward_read=3, index1=0, index2=0, reverse_read=3)

        summary = write_phix_csv(out_file, records, read_lengths)

        self.assertEqual(0.75, summary.error_rate_forward)
        self.assertEqual(0.25, summary.error_rate_reverse)
