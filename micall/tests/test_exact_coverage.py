"""
Tests for the exact_coverage tool.
"""

import unittest
from io import StringIO
from micall.utils.exact_coverage import (
    reverse_complement,
    build_kmer_index,
    find_exact_matches,
    read_fastq_pairs,
)


class TestReverseComplement(unittest.TestCase):
    def test_simple_sequence(self):
        seq = "ACGT"
        expected = "ACGT"
        result = reverse_complement(seq)
        self.assertEqual(expected, result)

    def test_reverse_complement(self):
        seq = "AAAA"
        expected = "TTTT"
        result = reverse_complement(seq)
        self.assertEqual(expected, result)

    def test_complex_sequence(self):
        seq = "ATCGATCG"
        expected = "CGATCGAT"
        result = reverse_complement(seq)
        self.assertEqual(expected, result)


class TestBuildKmerIndex(unittest.TestCase):
    def test_simple_contig(self):
        contigs = {"contig1": "ACGTACGT"}
        kmer_size = 3
        index = build_kmer_index(contigs, kmer_size)

        self.assertIn("ACG", index)
        self.assertIn("CGT", index)
        self.assertIn("GTA", index)
        self.assertIn("TAC", index)
        self.assertEqual(len(index["ACG"]), 2)  # ACG appears twice

    def test_skip_n_bases(self):
        contigs = {"contig1": "ACGNACGT"}
        kmer_size = 3
        index = build_kmer_index(contigs, kmer_size)

        # K-mers with N should not be indexed
        self.assertNotIn("CGN", index)
        self.assertNotIn("GNA", index)
        self.assertNotIn("NAC", index)


class TestFindExactMatches(unittest.TestCase):
    def test_exact_match(self):
        contigs = {"contig1": "ACGTACGTACGT"}
        kmer_index = build_kmer_index(contigs, 4)
        read_seq = "ACGT"

        matches = find_exact_matches(read_seq, kmer_index, contigs, 4, 3)

        # Should find 3 matches at positions 0, 4, and 8
        self.assertEqual(len(matches), 3)
        positions = [start for _, start, _ in matches]
        self.assertEqual(sorted(positions), [0, 4, 8])

    def test_no_match(self):
        contigs = {"contig1": "AAAAAAAAAA"}
        kmer_index = build_kmer_index(contigs, 4)
        read_seq = "CCCC"

        matches = find_exact_matches(read_seq, kmer_index, contigs, 4, 3)

        self.assertEqual(len(matches), 0)

    def test_partial_match_not_returned(self):
        contigs = {"contig1": "ACGTGGGG"}
        kmer_index = build_kmer_index(contigs, 4)
        read_seq = "ACGTTTTT"  # Only first 4 bases match

        matches = find_exact_matches(read_seq, kmer_index, contigs, 4, 3)

        # Should not match because the entire read doesn't match exactly
        self.assertEqual(len(matches), 0)


class TestReadFastqPairs(unittest.TestCase):
    def test_read_one_pair(self):
        fastq1 = StringIO("""\
@read1
ACGT
+
IIII
""")
        fastq2 = StringIO("""\
@read1
GGGG
+
IIII
""")

        pairs = list(read_fastq_pairs(fastq1, fastq2))

        self.assertEqual(len(pairs), 1)
        self.assertEqual(pairs[0], ("ACGT", "GGGG"))

    def test_read_multiple_pairs(self):
        fastq1 = StringIO("""\
@read1
ACGT
+
IIII
@read2
TTTT
+
IIII
""")
        fastq2 = StringIO("""\
@read1
GGGG
+
IIII
@read2
CCCC
+
IIII
""")

        pairs = list(read_fastq_pairs(fastq1, fastq2))

        self.assertEqual(len(pairs), 2)
        self.assertEqual(pairs[0], ("ACGT", "GGGG"))
        self.assertEqual(pairs[1], ("TTTT", "CCCC"))
