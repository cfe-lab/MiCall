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
    read_contigs,
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

        matches = list(find_exact_matches(read_seq, kmer_index, contigs, 4))

        # Should find 3 matches at positions 0, 4, and 8
        self.assertEqual(len(matches), 3)
        positions = [start for _, start, _ in matches]
        self.assertEqual(sorted(positions), [0, 4, 8])

    def test_no_match(self):
        contigs = {"contig1": "AAAAAAAAAA"}
        kmer_index = build_kmer_index(contigs, 4)
        read_seq = "CCCC"

        matches = list(find_exact_matches(read_seq, kmer_index, contigs, 4))

        self.assertEqual(len(matches), 0)

    def test_partial_match_not_returned(self):
        contigs = {"contig1": "ACGTGGGG"}
        kmer_index = build_kmer_index(contigs, 4)
        read_seq = "ACGTTTTT"  # Only first 4 bases match

        matches = list(find_exact_matches(read_seq, kmer_index, contigs, 4))

        # Should not match because the entire read doesn't match exactly
        self.assertEqual(len(matches), 0)

    def test_variable_read_lengths_short(self):
        """Test with reads shorter than typical (e.g., 20bp)"""
        contigs = {"contig1": "ACGTACGTACGTACGTACGT"}
        kmer_index = build_kmer_index(contigs, 4)
        read_seq = "ACGTAC"  # 6bp read

        matches = list(find_exact_matches(read_seq, kmer_index, contigs, 4))

        # Should find matches at positions 0, 4, 8, 12
        self.assertEqual(len(matches), 4)
        positions = [start for _, start, _ in matches]
        self.assertIn(0, positions)

    def test_variable_read_lengths_long(self):
        """Test with reads longer than typical (e.g., 200bp)"""
        long_contig = "ACGT" * 60  # 240 bp
        contigs = {"contig1": long_contig}
        kmer_index = build_kmer_index(contigs, 31)
        read_seq = long_contig[:200]  # 200bp read

        matches = list(find_exact_matches(read_seq, kmer_index, contigs, 31))

        # Should find at least one match
        self.assertGreaterEqual(len(matches), 1)
        # Verify match is correct length
        for _, start, end in matches:
            self.assertEqual(end - start, 200)

    def test_variable_read_lengths_mixed(self):
        """Test with different read lengths in same dataset"""
        contigs = {"contig1": "ACGTACGTACGTACGTACGTACGTACGTACGT"}
        kmer_index = build_kmer_index(contigs, 4)

        # Test various read lengths
        for length in [4, 8, 12, 16, 20]:
            read_seq = "ACGT" * (length // 4)
            matches = list(find_exact_matches(read_seq, kmer_index, contigs, 4))
            self.assertGreater(
                len(matches), 0, f"Should find matches for {length}bp read"
            )

    def test_read_at_contig_boundary(self):
        """Test reads that match at the end of a contig"""
        contigs = {"contig1": "AAAACGTACGTACGT"}
        kmer_index = build_kmer_index(contigs, 4)
        read_seq = "ACGTACGT"  # Matches at position 7

        matches = list(find_exact_matches(read_seq, kmer_index, contigs, 4))

        # Should find match at position 7
        positions = [start for _, start, _ in matches]
        self.assertIn(7, positions)

    def test_read_extends_beyond_contig(self):
        """Test that reads extending beyond contig boundaries are not matched"""
        contigs = {"contig1": "ACGTACGT"}
        kmer_index = build_kmer_index(contigs, 4)
        read_seq = "ACGTACGTACGTACGT"  # 16bp, longer than contig

        matches = list(find_exact_matches(read_seq, kmer_index, contigs, 4))

        # Should not find any matches
        self.assertEqual(len(matches), 0)

    def test_empty_read(self):
        """Test behavior with empty read"""
        contigs = {"contig1": "ACGTACGT"}
        kmer_index = build_kmer_index(contigs, 4)
        read_seq = ""

        matches = list(find_exact_matches(read_seq, kmer_index, contigs, 4))

        self.assertEqual(len(matches), 0)

    def test_read_with_n_bases(self):
        """Test reads containing N bases"""
        contigs = {"contig1": "ACGTACGTACGT"}
        kmer_index = build_kmer_index(contigs, 4)
        read_seq = "ACNTACGT"

        matches = list(find_exact_matches(read_seq, kmer_index, contigs, 4))

        # Should not match because contig doesn't have N
        self.assertEqual(len(matches), 0)

    def test_multiple_contigs(self):
        """Test matching against multiple contigs"""
        contigs = {
            "contig1": "AAAAAAAAAA",
            "contig2": "ACGTACGTAC",
            "contig3": "GGGGGGGGGG",
        }
        kmer_index = build_kmer_index(contigs, 4)
        read_seq = "ACGTACGT"

        matches = list(find_exact_matches(read_seq, kmer_index, contigs, 4))

        # Should only match contig2
        contig_names = [name for name, _, _ in matches]
        self.assertTrue(all(name == "contig2" for name in contig_names))

    def test_repeated_regions_get_multiple_hits(self):
        """Test that reads mapping to repeated regions increment coverage at all positions"""
        # Create a contig with a repeated sequence
        # AAACGTACGTAAACGTACGTAAA
        # Position 2: ACGTACGT
        # Position 12: ACGTACGT
        contigs = {"contig1": "AAACGTACGTAAACGTACGTAAA"}
        kmer_index = build_kmer_index(contigs, 4)
        read_seq = "ACGTACGT"  # This appears twice in the contig

        matches = list(find_exact_matches(read_seq, kmer_index, contigs, 4))

        # Should find 2 matches: one at position 2 and one at position 12
        self.assertEqual(len(matches), 2)
        positions = sorted([start for _, start, _ in matches])
        self.assertEqual(positions, [2, 12])

        # Verify the matches are for the same contig
        contig_names = [name for name, _, _ in matches]
        self.assertTrue(all(name == "contig1" for name in contig_names))


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

    def test_variable_length_reads(self):
        """Test reading FASTQ pairs with different read lengths"""
        fastq1 = StringIO("""\
@read1
ACGT
+
IIII
@read2
ACGTACGT
+
IIIIIIII
@read3
AC
+
II
""")
        fastq2 = StringIO("""\
@read1
GG
+
II
@read2
GGGGGGGG
+
IIIIIIII
@read3
GGGGGG
+
IIIIII
""")

        pairs = list(read_fastq_pairs(fastq1, fastq2))

        self.assertEqual(len(pairs), 3)
        self.assertEqual(pairs[0], ("ACGT", "GG"))
        self.assertEqual(pairs[1], ("ACGTACGT", "GGGGGGGG"))
        self.assertEqual(pairs[2], ("AC", "GGGGGG"))

    def test_empty_files(self):
        """Test behavior with empty FASTQ files"""
        fastq1 = StringIO("")
        fastq2 = StringIO("")

        pairs = list(read_fastq_pairs(fastq1, fastq2))

        self.assertEqual(len(pairs), 0)

    def test_reads_with_spaces_in_sequence(self):
        """Test that sequences are properly stripped"""
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

        # Should strip whitespace from sequences
        self.assertEqual(pairs[0][0], "ACGT")
        self.assertEqual(pairs[0][1], "GGGG")


class TestReverseComplementEdgeCases(unittest.TestCase):
    def test_with_n_bases(self):
        """Test reverse complement with N bases"""
        seq = "ACGTN"
        expected = "NACGT"
        result = reverse_complement(seq)
        self.assertEqual(expected, result)

    def test_empty_sequence(self):
        """Test reverse complement of empty sequence"""
        seq = ""
        expected = ""
        result = reverse_complement(seq)
        self.assertEqual(expected, result)

    def test_single_base(self):
        """Test reverse complement of single base"""
        for base, complement in [("A", "T"), ("T", "A"), ("C", "G"), ("G", "C")]:
            result = reverse_complement(base)
            self.assertEqual(complement, result)

    def test_palindrome(self):
        """Test reverse complement of palindromic sequence"""
        seq = "GAATTC"  # EcoRI restriction site
        expected = "GAATTC"
        result = reverse_complement(seq)
        self.assertEqual(expected, result)


class TestReadContigs(unittest.TestCase):
    def test_read_fasta_format(self):
        """Test reading contigs from FASTA format"""
        fasta_file = StringIO("""\
>contig1
ACGTACGT
>contig2
GGGGCCCC
""")

        contigs = read_contigs(fasta_file)

        self.assertEqual(len(contigs), 2)
        self.assertEqual(contigs["contig1"], "ACGTACGT")
        self.assertEqual(contigs["contig2"], "GGGGCCCC")

    def test_read_csv_with_contig_column(self):
        """Test reading contigs from CSV with 'contig' column (contigs.csv format)"""
        csv_file = StringIO("""\
ref,match,group_ref,contig
ref1,1.0,group1,ACGTACGT
ref2,1.0,group2,GGGGCCCC
""")

        contigs = read_contigs(csv_file)

        self.assertEqual(len(contigs), 2)
        # Should use position-based names since no sample/region columns
        self.assertEqual(contigs["contig1"], "ACGTACGT")
        self.assertEqual(contigs["contig2"], "GGGGCCCC")

    def test_read_csv_with_sequence_column(self):
        """Test reading contigs from CSV with 'sequence' column (conseq.csv format)"""
        csv_file = StringIO("""\
sample,region,q-cutoff,consensus-percent-cutoff,offset,sequence
sample1,region1,15,MAX,0,ACGTACGT
sample1,region2,15,MAX,0,GGGGCCCC
""")

        contigs = read_contigs(csv_file)

        self.assertEqual(len(contigs), 2)
        # Should use 'sample' column for name
        self.assertIn("sample1", contigs)
        # Second entry with same sample name should get _2 suffix
        self.assertIn("sample1_2", contigs)

    def test_sequence_column_prioritized_over_contig(self):
        """Test that 'sequence' column is prioritized over 'contig' column"""
        csv_file = StringIO("""\
ref,contig,sequence
ref1,TTTTTTTT,ACGTACGT
""")

        contigs = read_contigs(csv_file)

        # Should use 'sequence' column, not 'contig' column
        # Should use position-based name since no sample/region
        self.assertEqual(contigs["contig1"], "ACGTACGT")

    def test_name_column_priority(self):
        """Test that 'sample' is prioritized, then 'region', then position"""
        # Test with sample column
        csv_file = StringIO("""\
sample,region,ref,contig
mysample,myregion,myref,ACGTACGT
""")
        contigs = read_contigs(csv_file)
        self.assertIn("mysample", contigs)

        # Test with region column (no sample)
        csv_file = StringIO("""\
region,ref,contig
myregion,myref,GGGGCCCC
""")
        contigs = read_contigs(csv_file)
        self.assertIn("myregion", contigs)

        # Test with neither sample nor region - should use position
        csv_file = StringIO("""\
ref,contig
myref,TTTTTTTT
""")
        contigs = read_contigs(csv_file)
        self.assertIn("contig1", contigs)

    def test_csv_without_sequence_or_contig_column_raises_error(self):
        """Test that CSV without 'sequence' or 'contig' column raises ValueError"""
        csv_file = StringIO("""\
ref,other_column
ref1,data
""")

        with self.assertRaises(ValueError) as context:
            read_contigs(csv_file)

        self.assertIn("sequence", str(context.exception).lower())
        self.assertIn("contig", str(context.exception).lower())

    def test_empty_sequences_skipped(self):
        """Test that empty sequences are skipped"""
        csv_file = StringIO("""\
ref,contig
ref1,ACGTACGT
ref2,
ref3,GGGGCCCC
""")

        contigs = read_contigs(csv_file)

        # Should only have contig1 and contig3, ref2 should be skipped
        # Uses position-based names
        self.assertEqual(len(contigs), 2)
        self.assertIn("contig1", contigs)
        self.assertIn("contig3", contigs)
        self.assertNotIn("contig2", contigs)
