import csv
import os
from io import BytesIO
from io import StringIO
import unittest
from pathlib import Path

from micall.core.trim_fastqs import censor, trim
from micall.utils.translation import reverse_and_complement


class CensorTest(unittest.TestCase):
    def setUp(self):
        self.addTypeEqualityFunc(str, self.assertMultiLineEqual)
        self.original_bytes = b"""\
@M01841:45:000000000-A5FEG:1:1101:5296:13227 1:N:0:9
ACGT
+
AAAA
"""
        self.original_unicode = self.original_bytes.decode()
        self.original_file = BytesIO(self.original_bytes)
        self.bad_cycles = []
        self.censored_file = StringIO()
        self.summary_file = StringIO()
        self.summary_writer = csv.DictWriter(self.summary_file,
                                             ['avg_quality', 'base_count'],
                                             lineterminator=os.linesep)
        self.summary_writer.writeheader()

    def testNoBadCycles(self):
        expected_text = self.original_unicode

        censor(self.original_file,
               self.bad_cycles,
               self.censored_file,
               use_gzip=False)

        self.assertEqual(expected_text, self.censored_file.getvalue())

    def testBadCycle(self):
        self.bad_cycles = [{'tile': '1101', 'cycle': '3'}]
        expected_text = """\
@M01841:45:000000000-A5FEG:1:1101:5296:13227 1:N:0:9
ACNT
+
AA#A
"""

        censor(self.original_file,
               self.bad_cycles,
               self.censored_file,
               use_gzip=False)

        self.assertEqual(expected_text, self.censored_file.getvalue())

    def testBadTail(self):
        self.bad_cycles = [{'tile': '1101', 'cycle': '3'},
                           {'tile': '1101', 'cycle': '4'}]
        expected_text = """\
@M01841:45:000000000-A5FEG:1:1101:5296:13227 1:N:0:9
AC
+
AA
"""

        censor(self.original_file,
               self.bad_cycles,
               self.censored_file,
               use_gzip=False)

        self.assertEqual(expected_text, self.censored_file.getvalue())

    def testDifferentTile(self):
        self.bad_cycles = [{'tile': '1102', 'cycle': '3'}]
        expected_text = self.original_unicode

        censor(self.original_file,
               self.bad_cycles,
               self.censored_file,
               use_gzip=False)

        self.assertEqual(expected_text, self.censored_file.getvalue())

    def testDifferentDirection(self):
        self.original_bytes = b"""\
@M01841:45:000000000-A5FEG:1:1101:5296:13227 2:N:0:9
ACGT
+
AAAA
"""
        self.original_file = BytesIO(self.original_bytes)
        self.bad_cycles = [{'tile': '1101', 'cycle': '3'}]
        expected_text = self.original_bytes.decode()

        censor(self.original_file,
               self.bad_cycles,
               self.censored_file,
               use_gzip=False)

        self.assertEqual(expected_text, self.censored_file.getvalue())

    def testReverseDirection(self):
        self.original_bytes = b"""\
@M01841:45:000000000-A5FEG:1:1101:5296:13227 2:N:0:9
ACGT
+
AAAA
"""
        self.original_file = BytesIO(self.original_bytes)
        self.bad_cycles = [{'tile': '1101', 'cycle': '-3'}]
        expected_text = """\
@M01841:45:000000000-A5FEG:1:1101:5296:13227 2:N:0:9
ACNT
+
AA#A
"""

        censor(self.original_file,
               self.bad_cycles,
               self.censored_file,
               use_gzip=False)

        self.assertEqual(expected_text, self.censored_file.getvalue())

    def testTwoReads(self):
        self.original_bytes = b"""\
@M01841:45:000000000-A5FEG:1:1101:5296:13227 1:N:0:9
ACGT
+
AAAA
@M01841:45:000000000-A5FEG:1:1102:1234:12345 1:N:0:9
TGCA
+
BBBB
"""
        self.original_file = BytesIO(self.original_bytes)
        self.bad_cycles = [{'tile': '1101', 'cycle': '2'},
                           {'tile': '1102', 'cycle': '3'}]
        expected_text = """\
@M01841:45:000000000-A5FEG:1:1101:5296:13227 1:N:0:9
ANGT
+
A#AA
@M01841:45:000000000-A5FEG:1:1102:1234:12345 1:N:0:9
TGNA
+
BB#B
"""

        censor(self.original_file,
               self.bad_cycles,
               self.censored_file,
               use_gzip=False)

        self.assertEqual(expected_text, self.censored_file.getvalue())

    def testSummary(self):
        self.bad_cycles = [{'tile': '1101', 'cycle': '3'}]
        expected_summary = """\
avg_quality,base_count
32.0,4
"""

        censor(self.original_file,
               self.bad_cycles,
               self.censored_file,
               use_gzip=False,
               summary_writer=self.summary_writer)

        self.assertEqual(expected_summary, self.summary_file.getvalue())

    def testSummaryAverage(self):
        self.original_bytes = b"""\
@M01841:45:000000000-A5FEG:1:1101:5296:13227 1:N:0:9
ACGT
+
AACC
"""
        self.original_file = BytesIO(self.original_bytes)
        self.bad_cycles = [{'tile': '1101', 'cycle': '3'}]
        expected_summary = """\
avg_quality,base_count
33.0,4
"""

        censor(self.original_file,
               self.bad_cycles,
               self.censored_file,
               use_gzip=False,
               summary_writer=self.summary_writer)

        self.assertEqual(expected_summary, self.summary_file.getvalue())

    def testSummaryEmpty(self):
        self.original_bytes = b""
        self.original_file = BytesIO(self.original_bytes)
        expected_summary = """\
avg_quality,base_count
,0
"""

        censor(self.original_file,
               self.bad_cycles,
               self.censored_file,
               use_gzip=False,
               summary_writer=self.summary_writer)

        self.assertEqual(expected_summary, self.summary_file.getvalue())


def test_trim(tmpdir):
    core_path = Path(__file__).parent.parent / 'core'
    adapters_read1_path = core_path / 'adapters_read1.fasta'
    adapters_read2_path = core_path / 'adapters_read1.fasta'
    adapter_read1 = adapters_read1_path.read_text().splitlines()[1]
    adapter_read2 = adapters_read2_path.read_text().splitlines()[1]
    read1_content = 'ACTCTACTAACTGTCGTGTTT'
    read2_content = reverse_and_complement(read1_content)
    read1_extra = 'GGGAAATTT'
    read2_extra = 'CCCAAATTT'
    content_quality = 'A' * len(read1_content)
    adapter_quality = 'A' * len(adapter_read1)
    extra_quality = 'A' * len(read1_extra)

    tmp_path = Path(tmpdir)
    fastq1_path = tmp_path / 'read1.fastq'
    fastq2_path = tmp_path / 'read2.fastq'
    trimmed1_path = tmp_path / 'trimmed1.fastq'
    trimmed2_path = tmp_path / 'trimmed2.fastq'
    fastq1_path.write_text(f'''\
@pair1::::tile1:: 1:::
{read1_content}
+
{content_quality}
@pair2::::tile1:: 1:::
{read1_content}{read1_extra}
+
{content_quality}{extra_quality}
@pair3::::tile1:: 1:::
{read1_content}{adapter_read1}{read1_extra}
+
{content_quality}{adapter_quality}{extra_quality}
''')
    fastq2_path.write_text(f'''\
@pair1::::tile1:: 2:::
{read2_content}
+
{content_quality}
@pair2::::tile1:: 2:::
{read2_content}{read2_extra}
+
{content_quality}{extra_quality}
@pair3::::tile1:: 2:::
{read2_content}{adapter_read2}{read2_extra}
+
{content_quality}{adapter_quality}{extra_quality}
''')
    expected_trimmed1 = f'''\
@pair1::::tile1:: 1:::
{read1_content}
+
{content_quality}
@pair2::::tile1:: 1:::
{read1_content}{read1_extra}
+
{content_quality}{extra_quality}
@pair3::::tile1:: 1:::
{read1_content}
+
{content_quality}
'''
    expected_trimmed2 = f'''\
@pair1::::tile1:: 2:::
{read2_content}
+
{content_quality}
@pair2::::tile1:: 2:::
{read2_content}{read2_extra}
+
{content_quality}{extra_quality}
@pair3::::tile1:: 2:::
{read2_content}{adapter_read2}{read2_extra}
+
{content_quality}{adapter_quality}{extra_quality}
'''

    trim([fastq1_path, fastq2_path],
         'no_bad_cycles.csv',
         [str(trimmed1_path), str(trimmed2_path)],
         use_gzip=False)

    trimmed1 = trimmed1_path.read_text()
    trimmed2 = trimmed2_path.read_text()
    assert trimmed1 == expected_trimmed1
    assert trimmed2 == expected_trimmed2
