import csv
import os
from io import BytesIO
from io import StringIO
import unittest
from pathlib import Path

import pytest

from micall.core.trim_fastqs import censor, trim, cut_all
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
        """ Bad cycle doesn't match this read. """
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
               use_gzip=False,
               cycle_sign=-1)

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
               use_gzip=False,
               cycle_sign=-1)

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
    read1_content = 'TATCTACTAACTGTCGGTCTAC'
    read2_content = reverse_and_complement(read1_content)
    expected1 = build_fastq(read1_content)
    expected2 = build_fastq(read2_content)

    tmp_path = Path(tmpdir)
    fastq1_path = tmp_path / 'read1.fastq'
    fastq2_path = tmp_path / 'read2.fastq'
    trimmed1_path = tmp_path / 'trimmed1.fastq'
    trimmed2_path = tmp_path / 'trimmed2.fastq'
    fastq1_path.write_text(expected1)
    fastq2_path.write_text(expected2)

    trim([fastq1_path, fastq2_path],
         'no_bad_cycles.csv',
         [str(trimmed1_path), str(trimmed2_path)],
         use_gzip=False)

    trimmed1 = trimmed1_path.read_text()
    trimmed2 = trimmed2_path.read_text()
    assert trimmed1 == expected1
    assert trimmed2 == expected2


@pytest.mark.parametrize(
    "scenario,read1,read2,expected1,expected2",
    [
     ('no adapter',
      'TGGAAGGGCTAATTCACTCCCAACG',
      # REF
      'CGTTGGGAGTGAATTAGCCCTTCCA',
      # rev(REF)
      'TGGAAGGGCTAATTCACTCCCAACG',
      # unchanged
      'CGTTGGGAGTGAATTAGCCCTTCCA',
      # unchanged
      ),
     ('full adapters',
      'TGGAAGGGCTAATTCACTCCCAACGCTGTCTCTTATACACATCTCCGAGCCCACGAGAC',
      # REF                    ][ rev(ADAPTER2)
      'CGTTGGGAGTGAATTAGCCCTTCCACTGTCTCTTATACACATCTGACGCTGCCGACGA',
      # rev(REF)               ][ rev(ADAPTER1)
      'TGGAAGGGCTAATTCACTCCCAACG',
      # REF
      'CGTTGGGAGTGAATTAGCCCTTCCA',
      # rev(REF)
      ),
     ('full adapters plus garbage',
      'TGGAAGGGCTAATTCACTCCCAACGCTGTCTCTTATACACATCTCCGAGCCCACGAGACCAGTACGCA',
      # REF                    ][ rev(ADAPTER2)                  ][ garbage
      'CGTTGGGAGTGAATTAGCCCTTCCACTGTCTCTTATACACATCTGACGCTGCCGACGAAAGTAGCAAC',
      # rev(REF)               ][ rev(ADAPTER1)                 ][ garbage
      'TGGAAGGGCTAATTCACTCCCAACG',
      # REF
      'CGTTGGGAGTGAATTAGCCCTTCCA',
      # rev(REF)
      ),
     ('partial adapters',
      'TGGAAGGGCTAATTCACTCCCAACGCTGTCTCTTATACACATCTCCGAG',
      # REF                    ][ partial rev(ADAPTER2)
      'CGTTGGGAGTGAATTAGCCCTTCCACTGTCTCTTATACACATCTGACGC',
      # rev(REF)               ][ partial rev(ADAPTER1)
      'TGGAAGGGCTAATTCACTCCCAACG',
      # REF
      'CGTTGGGAGTGAATTAGCCCTTCCA',
      # rev(REF)
      ),
     ('partial adapters plus garbage',
      'TGGAAGGGCTAATTCACTCCCAACGCTGTCTCTTATACACATCTCCGAGCCAGTACGCA',
      # REF                    ][ partial rev(ADAPTER2)][ garbage
      'CGTTGGGAGTGAATTAGCCCTTCCACTGTCTCTTATACACATCTGACGCAAGTAGCAAC',
      # rev(REF)               ][ partial rev(ADAPTER1)][ garbage
      'TGGAAGGGCTAATTCACTCCCAACGCTGTCTCTTATACACATCTCCGAGCCAGTACGCA',
      # unchanged, because partial adapters only trimmed off the end
      'CGTTGGGAGTGAATTAGCCCTTCCACTGTCTCTTATACACATCTGACGCAAGTAGCAAC',
      # unchanged
      ),
     ('no primers',
      'TGGAAGGGCTAATTCACTCCCAACG',
      # REF
      'CGTTGGGAGTGAATTAGCCCTTCCA',
      # rev(REF)
      'TGGAAGGGCTAATTCACTCCCAACG',
      # unchanged
      'CGTTGGGAGTGAATTAGCCCTTCCA',
      # unchanged
      ),
     ('full right primers',
      'TGGAAGGGCTAATTCACTCCCAACGGAGGCACGTCAACATCTTAAAGATG',
      # REF                    ][ rev(RIGHT)
      'CATCTTTAAGATGTTGACGTGCCTCCGTTGGGAGTGAATTAGCCCTTCCA',
      # RIGHT                  ][ rev(REF)
      'TGGAAGGGCTAATTCACTCCCAACG',
      # REF
      'CGTTGGGAGTGAATTAGCCCTTCCA',
      # rev(REF)
      ),
     ('partial right primers',
      'TGGAAGGGCTAATTCACTCCCAACGGAGGCACGTCAACATCTTAA',
      # REF                    ][ partial rev(RIGHT)
      'TTAAGATGTTGACGTGCCTCCGTTGGGAGTGAATTAGCCCTTCCA',
      # partial RIGHT     ][ rev(REF)
      'TGGAAGGGCTAATTCACTCCCAACG',
      # REF
      'CGTTGGGAGTGAATTAGCCCTTCCA',
      # rev(REF)
      ),
     ('full left primers',
      'ACCAACCAACTTTCGATCTCTTGTTGGAAGGGCTAATTCACTCCCAACG',
      # LEFT                  ][ REF
      'CGTTGGGAGTGAATTAGCCCTTCCAACAAGAGATCGAAAGTTGGTTGGT',
      # rev(REF)               ][ rev(LEFT)
      'TGGAAGGGCTAATTCACTCCCAACG',
      # REF
      'CGTTGGGAGTGAATTAGCCCTTCCA',
      # rev(REF)
      ),
     ('partial left primers',
      'ACCAACTTTCGATCTCTTGTTGGAAGGGCTAATTCACTCCCAACG',
      # partial LEFT      ][ REF
      'CGTTGGGAGTGAATTAGCCCTTCCAACAAGAGATCGAAAGTTGG',
      # rev(REF)               ][ rev(partial LEFT)
      'TGGAAGGGCTAATTCACTCCCAACG',
      # REF
      'CGTTGGGAGTGAATTAGCCCTTCCA',
      # rev(REF)
      ),
     ('partial left primers plus garbage',
      'CATAAGGATACCAACTTTCGATCTCTTGTTGGAAGGGCTAATTCACTCCCAACG',
      # garbage][ partial LEFT     ][ REF
      'CGTTGGGAGTGAATTAGCCCTTCCAACAAGAGATCGAAAGTTGGATCCTTATG',
      # rev(REF)               ][ rev(part LEFT)  ][garbage
      'CATAAGGATACCAACTTTCGATCTCTTGTTGGAAGGGCTAATTCACTCCCAACG',
      # unchanged
      'CGTTGGGAGTGAATTAGCCCTTCCAACAAGAGATCGAAAGTTGGATCCTTATG',
      # unchanged
      ),
     ('full left primers plus garbage',
      'CATAAGGATACCAACCAACTTTCGATCTCTTGTTGGAAGGGCTAATTCACTCCCAACG',
      # garbage][ LEFT                 ][ REF
      'CGTTGGGAGTGAATTAGCCCTTCCAACAAGAGATCGAAAGTTGGTTGGTATCCTTATG',
      # rev(REF)               ][ rev(LEFT)            ][garbage
      'CATAAGGATACCAACCAACTTTCGATCTCTTGTTGGAAGGGCTAATTCACTCCCAACG',
      # unchanged
      'CGTTGGGAGTGAATTAGCCCTTCCAACAAGAGATCGAAAGTTGGTTGGTATCCTTATG',
      # unchanged
      ),
     ('full right primers plus garbage',
      'TGGAAGGGCTAATTCACTCCCAACGGAGGCACGTCAACATCTTAAAGATGATGCACTT',
      # REF                    ][ rev(RIGHT)            ][garbage
      'TACCGGACTCATCTTTAAGATGTTGACGTGCCTCCGTTGGGAGTGAATTAGCCCTTCCA',
      # garbage][ RIGHT                 ][ rev(REF)
      'TGGAAGGGCTAATTCACTCCCAACGGAGGCACGTCAACATCTTAAAGATGATGCACTT',
      # unchanged
      'TACCGGACTCATCTTTAAGATGTTGACGTGCCTCCGTTGGGAGTGAATTAGCCCTTCCA',
      # unchanged
      ),
     ('left and right primers',
      'ACCAACCAACTTTCGATCTCTTGTTGGAAGGGCTAATTCACTCCCAACGGAGGCACGTCAACATCTTAAAGATG',
      # LEFT                  ][ REF                   ][ rev(RIGHT)            ]
      'CATCTTTAAGATGTTGACGTGCCTCCGTTGGGAGTGAATTAGCCCTTCCAACAAGAGATCGAAAGTTGGTTGGT',
      # RIGHT                  ][ rev(REF)              ][ rev(LEFT)
      'TGGAAGGGCTAATTCACTCCCAACG',
      # REF
      'CGTTGGGAGTGAATTAGCCCTTCCA',
      # rev(REF)
      ),
     ('right primers plus adapters',
      'TGGAAGGGCTAATTCACTCCCAACGGAGGCACGTCAACATCTTAAAGATGCTGTCTCTTATACACATCTCCGAGCCCACGAGAC',
      # REF                    ][ rev(RIGHT)            ][ rev(ADAPTER2)
      'CATCTTTAAGATGTTGACGTGCCTCCGTTGGGAGTGAATTAGCCCTTCCACTGTCTCTTATACACATCTGACGCTGCCGACGA',
      # RIGHT                  ][ rev(REF)              ][ rev(ADAPTER1)
      'TGGAAGGGCTAATTCACTCCCAACG',
      # REF
      'CGTTGGGAGTGAATTAGCCCTTCCA',
      # rev(REF)
      ),
     ('primer dimer',
      'TGGAAATACCCACAAGTTAATGGTTTAACAGGCACAGGTGTCTGTCTCTTATACACATCTCCGAGCCCACGAGACACTACCTGGAA',
      # nCoV-2019_18_LEFT          ]            [ rev(ADAPTER2)                  ][ garbage
      #                   [ rev(..._76_RIGHT)  ]
      'ACACCTGTGCCTGTTAAACCATTAACTTGTGGGTATTTCCACTGTCTCTTATACACATCTGACGCTGCCGACGAAGGTTCTCAGGA',
      # nCoV-2019_76_RIGHT  ]                   [ rev(ADAPTER1)                 ][ garbage
      #            [ rev(nCoV-2019_18_LEFT)    ]
      '',
      # Trimmed to nothing
      '',
      # Trimmed to nothing
      ),
     ('primer dimer with partial right match',
      'TGGCTATTGATTATAAACACTACACACCCTGCACAAGAAAAGAACTTCACCTGTCTCTTATACACATCTCCGAGCCCACGAGACACTACCTGGAA',
      #  nCoV-2019_21_LEFT         ]                     [ rev(ADAPTER2)                  ][ garbage ]
      #                          [ rev(..._81_RIGHT)    ] last 5 of 81_RIGHT match start of HCV_Pr3
      'GTGAAGTTCTTTTCTTGTGCAGGGTGTGTAGTGTTTATAATCAATAGCCACTGTCTCTTATACACATCTGACGCTGCCGACGAAGGTTCTCAGGA',
      # nCoV-2019_81_RIGHT    ]                          [ rev(ADAPTER1)                 ][ garbage
      #                     [ rev(nCoV-2019_21_LEFT)    ]
      '',
      # Trimmed to nothing
      '',
      # Trimmed to nothing
      )
     ])
def test_cut_adapters(tmpdir: str,
                      scenario: str,
                      read1: str,
                      read2: str,
                      expected1: str,
                      expected2: str):
    """ Cut adapter sequence from a read pair.

    The reference section is pulled from the start of HXB2:
    TGGAAGGGCTAATTCACTCCCAACG
    Reverse complement of that is:
    CGTTGGGAGTGAATTAGCCCTTCCA
    Nextera Read 1 adapter:
    TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG
    Reverse complement:
    CTGTCTCTTATACACATCTGACGCTGCCGACGA
    Nextera Read 2 adapter:
    GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG
    Reverse complement:
    CTGTCTCTTATACACATCTCCGAGCCCACGAGAC
    Left primer:
    ACCAACCAACTTTCGATCTCTTGT
    Reverse complement:
    ACAAGAGATCGAAAGTTGGTTGGT
    Right primer (matches reverse complement of reference):
    CATCTTTAAGATGTTGACGTGCCTC
    Reverse complement (matches forward reference):
    GAGGCACGTCAACATCTTAAAGATG
    """
    tmp_path = Path(tmpdir)
    fastq1_path = tmp_path / 'read1.fastq'
    fastq2_path = tmp_path / 'read2.fastq'
    trimmed1_path = tmp_path / 'trimmed1.fastq'
    trimmed2_path = tmp_path / 'trimmed2.fastq'
    fastq1_path.write_text(build_fastq(read1))
    fastq2_path.write_text(build_fastq(read2))
    expected_trimmed1 = build_fastq(expected1)
    expected_trimmed2 = build_fastq(expected2)

    cut_all(fastq1_path, fastq2_path, trimmed1_path, trimmed2_path)

    assert trimmed1_path.read_text() == expected_trimmed1
    assert trimmed2_path.read_text() == expected_trimmed2


def build_fastq(read_sequence):
    # Write two reads in the file to test primer dimer caching.
    expected_quality1 = 'A' * len(read_sequence)
    expected_trimmed1 = f'''\
@pair1
{read_sequence}
+
{expected_quality1}
@pair2
{read_sequence}
+
{expected_quality1}
'''
    return expected_trimmed1
