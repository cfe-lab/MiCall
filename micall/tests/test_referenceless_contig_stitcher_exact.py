"""
This module includes exact tests for the referenceless contig stitcher.
These tests produce detailed logs of every action taken by the stitcher and assert an
exact match against previously recorded outputs. These tests are intentionally
brittle to catch any behavioural changes from refactoring or other
non-semantic modifications.

If these brittle tests fail due to unexpected log changes, you may revert
the logs to their prior state with:

    git checkout HEAD micall/tests/data
"""

import pytest
from pathlib import Path
from typing import Iterator, Iterable

from micall.utils.contig_stitcher_context import ReferencelessStitcherContext

# Load autouse fixtures
from micall.tests.referenceless_tests_utils import disable_acceptable_prob_check, log_check, random_fasta_file, run_full_pipeline, load_projects, disable_kmer_filter, disable_min_overlap_size

# to avoid linter warnings
assert disable_acceptable_prob_check is not None
assert log_check is not None
assert random_fasta_file is not None
assert load_projects is not None
assert disable_kmer_filter is not None
assert disable_min_overlap_size is not None


def params(good: Iterable[int], bad: Iterable[object], reason_fmt: str) -> Iterator[object]:
    bad = frozenset(bad)

    for testcase in sorted(good):
        if testcase in bad:
            reason = reason_fmt.format(testcase=testcase)
            yield pytest.param(testcase, marks=pytest.mark.xfail(reason=reason, strict=True))
        else:
            yield testcase


# TODO: ensure that every random seed can be stitched.
@pytest.mark.parametrize("random_seed", params(range(50), [8, 13, 15, 29, 32, 36, 42], "Probably gaps that are too small."))
def test_full_pipeline_small_values(log_check, tmp_path: Path, random_fasta_file, random_seed: int, monkeypatch, disable_acceptable_prob_check):
    monkeypatch.setattr("micall.utils.referenceless_contig_stitcher.MAX_ALTERNATIVES", 1)
    assert not ReferencelessStitcherContext.get().is_debug2
    ReferencelessStitcherContext.get().is_debug2 = True
    converted_fasta_file, ref_seqs = random_fasta_file(6, random_seed)
    run_full_pipeline(log_check, tmp_path, converted_fasta_file, ref_seqs)


# TODO: ensure that every random seed can be stitched.
@pytest.mark.parametrize("random_seed", params(range(999), [8, 23, 57, 59, 64, 67, 88, 95, 116, 122, 146, 180, 233, 282, 324, 331, 335, 342, 473, 492, 494, 520, 531, 552, 561, 569, 581, 631, 639, 666, 732, 860, 861, 874, 876, 956, 966, 988], "Probably gaps that are too small."))
def test_full_pipeline_tiny_values(log_check, tmp_path: Path, random_fasta_file, random_seed: int, monkeypatch, disable_acceptable_prob_check):
    monkeypatch.setattr("micall.utils.referenceless_contig_stitcher.MAX_ALTERNATIVES", 1)
    assert not ReferencelessStitcherContext.get().is_debug2
    ReferencelessStitcherContext.get().is_debug2 = True
    converted_fasta_file, ref_seqs = random_fasta_file(3, random_seed)
    run_full_pipeline(log_check, tmp_path, converted_fasta_file, ref_seqs)


@pytest.mark.parametrize("random_seed", params(range(10), [], "Probably gaps that are too small."))
def test_full_pipeline(log_check, tmp_path: Path, random_fasta_file, random_seed: int, monkeypatch):
    monkeypatch.setattr("micall.utils.referenceless_contig_stitcher.MIN_MATCHES", 70)
    assert not ReferencelessStitcherContext.get().is_debug2
    converted_fasta_file, ref_seqs = random_fasta_file(99, random_seed)
    run_full_pipeline(log_check, tmp_path, converted_fasta_file, ref_seqs)
