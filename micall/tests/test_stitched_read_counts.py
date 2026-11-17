"""Sequence-based remap evidence survives referenceless stitching changes."""
from io import StringIO
import unittest

from micall.utils.referencefull_contig_stitcher import read_contigs
from micall.utils.stitched_read_counts import project_stitched_counts


def contigs(*rows):
    lines = ["ref,match,group_ref,contig"]
    lines += [f"{ref},1,{ref},{sequence}" for ref, sequence in rows]
    return StringIO("\n".join(lines) + "\n")


def remap(*rows):
    lines = ["type,count"]
    lines += [f"remap {i}-{ref},{count}" for i, ref, count in rows]
    return StringIO("\n".join(lines) + "\n")


class StitchedReadCountsTests(unittest.TestCase):
    def test_remapper_partial_and_unknown_names_keep_read_counts(self):
        from micall.core.remap import get_contig_name

        def named_rows():
            return StringIO(
                "ref,match,group_ref,contig\n"
                "HIV1-B-FR-K03455-seed,0.1,HIV1-B-FR-K03455-seed,AAAAAC\n"
                "unknown,0.1,unknown,CCCCCA\n"
            )

        names = [
            get_contig_name(1, "HIV1-B-FR-K03455-seed", False, False),
            get_contig_name(2, "unknown", False, False),
        ]
        self.assertTrue(all(name.endswith("-partial") for name in names))
        counts = StringIO(
            "type,count\n"
            + "".join(f"remap {name},{depth}\n"
                      for name, depth in zip(names, (100, 200)))
        )
        projected = project_stitched_counts(named_rows(), named_rows(), counts)
        self.assertEqual(projected,
                         {1: 100, 2: 200})

    def test_referencefull_cli_reader_uses_actual_remapper_names(self):
        from micall.core.remap import get_contig_name

        def rows():
            return StringIO(
                "ref,match,group_ref,contig\n"
                "HIV1-B-FR-K03455-seed,0.1,HIV1-B-FR-K03455-seed,AAAAAC\n"
                "unknown,0.1,unknown,CCCCCA\n"
            )

        names = [
            get_contig_name(1, "HIV1-B-FR-K03455-seed", False, False),
            get_contig_name(2, "unknown", False, False),
        ]
        counts = {names[0]: 31, names[1]: 59}
        self.assertEqual(
            [c.reads_count for c in read_contigs(rows(), counts)],
            [31, 59])
        self.assertEqual(
            [c.reads_count for c in read_contigs(
                rows(), {}, row_counts={1: 31, 2: 59})],
            [31, 59])

    def test_referencefull_cli_counts_and_native_row_counts_are_exclusive(self):
        from micall.utils.referencefull_contig_stitcher import (
            referencefull_contig_stitcher,
        )

        self.assertEqual(referencefull_contig_stitcher(
            contigs(("main", "A" * 160)), None, None,
            remap((1, "main", 42))), 1)
        with self.assertRaisesRegex(ValueError, "either"):
            referencefull_contig_stitcher(
                contigs(("main", "A" * 160)), None, None,
                remap((1, "main", 42)), row_counts={1: 42})
        self.assertEqual(referencefull_contig_stitcher(
            contigs(("main", "A" * 160)), None, None,
            row_counts={1: 42}), 1)

    def test_excluded_seed_name_is_retained_and_not_guessed(self):
        from micall.core.remap import get_contig_name

        def rows():
            return StringIO(
                "ref,match,group_ref,contig\n"
                "HIV1-B-FR-K03455-seed,0.9,HIV1-B-FR-K03455-seed,AAAAAC\n"
            )

        seed = "HIV1-B-FR-K03455-seed"
        excluded = get_contig_name(1, seed, True, False, {seed})
        ordinary = get_contig_name(1, seed, True, False)
        counts = StringIO(f"type,count\nremap {excluded},53\n")
        self.assertEqual(project_stitched_counts(rows(), rows(), counts),
                         {1: 53})
        # Two incompatible names for the same source row are ambiguous.
        conflicting = StringIO(
            f"type,count\nremap {excluded},53\nremap {ordinary},92\n"
        )
        self.assertEqual(project_stitched_counts(
            rows(), rows(), conflicting), {})

    def test_reordered_unique_sequences_keep_exact_support(self):
        original = contigs(("main", "A" * 324), ("small", "C" * 251),
                           ("other", "G" * 306))
        stitched = contigs(("other", "G" * 306), ("main", "A" * 324),
                           ("small", "C" * 251))
        counts = remap((1, "main", 412), (2, "small", 39),
                       (3, "other", 119))
        projected = project_stitched_counts(original, stitched, counts)
        self.assertEqual(projected,
                         {1: 119, 2: 412, 3: 39})

    def test_duplicate_originals_may_share_reads_do_not_sum_counts(self):
        original = contigs(("main", "A" * 324), ("small", "C" * 251),
                           ("other", "G" * 306), ("main", "A" * 324))
        stitched = contigs(("main", "A" * 324), ("other", "G" * 306),
                           ("small", "C" * 251))
        counts = remap((1, "main", 412), (2, "small", 39),
                       (3, "other", 119), (4, "main", 392))
        # 412+392 is not necessarily the union of read IDs when reads
        # can multi-map. An unknown count is safer than the false 804.
        projected = project_stitched_counts(original, stitched, counts)
        self.assertEqual(projected,
                         {2: 119, 3: 39})

    def test_duplicate_stitched_sequences_cannot_copy_one_count_twice(self):
        original = contigs(("main", "AAACTG" * 30))
        stitched = contigs(("left", "AAACTG" * 30),
                           ("right", "CAGTTT" * 30))
        projected = project_stitched_counts(
            original, stitched, remap((1, "main", 157)))
        self.assertEqual(projected, {})

    def test_stitched_join_does_not_inherit_an_unrelated_row_count(self):
        original = contigs(("left", "A" * 160), ("right", "C" * 160))
        stitched = contigs(("joined", "A" * 160 + "C" * 160),
                           ("right", "C" * 160))
        projected = project_stitched_counts(
            original, stitched, remap((1, "left", 90), (2, "right", 70)))
        self.assertEqual(projected, {2: 70})

    def test_missing_original_remap_count_is_not_silently_zero(self):
        original = contigs(("left", "A" * 160), ("right", "C" * 160))
        stitched = contigs(("right", "C" * 160), ("left", "A" * 160))
        projected = project_stitched_counts(
            original, stitched, remap((1, "left", 90)))
        self.assertEqual(projected, {2: 90})

    def test_duplicate_with_missing_count_cannot_claim_complete_evidence(self):
        original = contigs(("a", "A" * 160), ("b", "A" * 160))
        stitched = contigs(("a", "A" * 160))
        projected = project_stitched_counts(
            original, stitched, remap((1, "a", 90)))
        self.assertEqual(projected, {})

    def test_reverse_complement_of_original_contig_keeps_its_count(self):
        original = contigs(("forward", "AAACTG" * 30))
        stitched = contigs(("reverse", "CAGTTT" * 30))
        projected = project_stitched_counts(
            original, stitched, remap((1, "forward", 157)))
        self.assertEqual(projected, {1: 157})

    def test_reverse_complement_duplicates_are_ambiguous(self):
        original = contigs(("forward", "AAACTG" * 30),
                           ("reverse", "CAGTTT" * 30))
        stitched = contigs(("merged", "AAACTG" * 30))
        projected = project_stitched_counts(
            original, stitched,
            remap((1, "forward", 20), (2, "reverse", 30)))
        self.assertEqual(projected, {})

    def test_remap_final_counts_override_prior_remap_counts(self):
        original = contigs(("a", "A" * 160))
        stitched = contigs(("a", "A" * 160))
        source = StringIO("type,count\nremap 1-a,12\nremap-final 1-a,15\n")
        projected = project_stitched_counts(original, stitched, source)
        self.assertEqual(projected, {1: 15})


if __name__ == "__main__":
    unittest.main()
