from dataclasses import dataclass
import logging
import os
import pytest
from typing import Tuple, List

from aligntools import CigarActions, CigarHit, Cigar

import micall.utils.registry as registry
import micall.utils.referencefull_contig_stitcher as stitcher
from micall.utils.referencefull_contig_stitcher import (
    split_contigs_with_gaps,
    stitch_contigs,
    GenotypedContig,
    merge_intervals,
    find_covered_contig,
    stitch_consensus,
    align_all_to_reference,
    lstrip,
    rstrip,
)
from micall.core.plot_contigs import plot_stitcher_coverage
from micall.tests.utils import mock_align_consensus, MockAlignment
from micall.tests.test_fasta_to_csv import (
    check_hcv_db,
    DEFAULT_DATABASE,
)  # activates the fixture
from micall.tests.test_remap import load_projects  # activates the "projects" fixture


logging.getLogger("micall.utils.referencefull_contig_stitcher").setLevel(logging.DEBUG)
logging.getLogger("micall.core.plot_contigs").setLevel(logging.DEBUG)


# make linters not complain about unused imports.
assert check_hcv_db is not None
assert DEFAULT_DATABASE is not None
assert load_projects is not None


@pytest.fixture()
def exact_aligner(monkeypatch):
    monkeypatch.setattr(
        "micall.utils.referencefull_contig_stitcher.align_consensus",
        mock_align_consensus
    )


@pytest.fixture
def visualizer(request, tmp_path):
    stitcher.ReferencefullStitcherContext.set(stitcher.ReferencefullStitcherContext())
    registry.set(registry.Registry())
    test_name = request.node.name
    plot_name = test_name + ".svg"
    pwd = os.path.dirname(__file__)
    plots_dir = os.path.join(pwd, "data", "stitcher_plots")
    os.makedirs(plots_dir, exist_ok=True)
    path_to_expected = os.path.join(plots_dir, plot_name)
    path_to_produced = os.path.join(tmp_path, plot_name)

    def check():
        logs = stitcher.ReferencefullStitcherContext.get().events
        figure = plot_stitcher_coverage(logs, path_to_produced)

        with open(path_to_produced, "r") as produced_file:
            produced_data = produced_file.read()
        with open(path_to_expected, "r") as expected_file:
            expected_data = expected_file.read()

        assert produced_data == expected_data, (
            "The contents of the stitched contigs plot"
            " does not match the expected contents."
        )

        return figure

    return check


def test_identical_stitching_of_one_contig(exact_aligner, visualizer):
    # Scenario: When stitching one contig, it remains the same.

    contigs = [
        GenotypedContig(
            name="a",
            seq="ACTGACTG" * 100,
            ref_name="testref",
            group_ref="testref",
            ref_seq="T" * 20 + "ACTGACTG" * 110 + "T" * 20,
            match_fraction=1.0,
            reads_count=None,
        ),
    ]

    results = list(stitch_contigs(contigs))
    assert len(results) == 1
    assert results[0].seq == contigs[0].seq

    assert len(visualizer().elements) > len(contigs)


def test_separate_stitching_of_non_overlapping_contigs_1(exact_aligner, visualizer):
    # Scenario: When stitching multiple non-overlapping contigs, the order doesn't matter.

    ref_seq = "A" * 100 + "C" * 100

    contigs = [
        GenotypedContig(
            name="a",
            seq="A" * 70,
            ref_name="testref",
            group_ref="testref",
            ref_seq=ref_seq,
            match_fraction=0.5,
            reads_count=None,
        ),
        GenotypedContig(
            name="b",
            seq="C" * 70,
            ref_name="testref",
            group_ref="testref",
            ref_seq=ref_seq,
            match_fraction=0.5,
            reads_count=None,
        ),
    ]

    results = list(stitch_contigs(contigs))

    # No claims about the output order, so wrap into set()
    assert {contig.seq for contig in contigs} == {contig.seq for contig in results}

    assert len(visualizer().elements) > len(contigs)


def test_separate_stitching_of_non_overlapping_contigs_2(exact_aligner, visualizer):
    # Scenario: When stitching multiple non-overlapping contigs,
    # the order doesn't matter.

    ref_seq = "A" * 100 + "C" * 100

    contigs = [
        GenotypedContig(
            name="b",
            seq="C" * 70,
            ref_name="testref",
            group_ref="testref",
            ref_seq=ref_seq,
            match_fraction=0.5,
            reads_count=None,
        ),
        GenotypedContig(
            name="a",
            seq="A" * 70,
            ref_name="testref",
            group_ref="testref",
            ref_seq=ref_seq,
            match_fraction=0.5,
            reads_count=None,
        ),
    ]

    results = list(stitch_contigs(contigs))

    # No claims about the output order, so wrap into set()
    assert {contig.seq for contig in contigs} == {contig.seq for contig in results}

    assert len(visualizer().elements) > len(contigs)


def test_correct_stitching_of_two_partially_overlapping_contigs(
    exact_aligner, visualizer
):
    # Scenario: Two partially overlapping contigs are stitched correctly into a single sequence.

    ref_seq = "A" * 100 + "C" * 100

    contigs = [
        GenotypedContig(
            name="a",
            seq="A" * 50 + "C" * 20,
            ref_name="testref",
            group_ref="testref",
            ref_seq=ref_seq,
            match_fraction=0.5,
            reads_count=None,
        ),
        GenotypedContig(
            name="b",
            seq="A" * 20 + "C" * 50,
            ref_name="testref",
            group_ref="testref",
            ref_seq=ref_seq,
            match_fraction=0.5,
            reads_count=None,
        ),
    ]

    results = list(stitch_contigs(contigs))
    assert len(results) == 1

    result = results[0]

    assert 100 == len(result.seq)
    assert result.seq == "A" * 50 + "C" * 50

    assert len(visualizer().elements) > len(contigs)


def test_correct_stitching_of_two_partially_overlapping_contigs_with_padding(
    exact_aligner, visualizer
):
    # Scenario: Two partially overlapping contigs are stitched correctly into a single sequence.

    ref_seq = "A" * 100 + "C" * 100

    contigs = [
        GenotypedContig(
            name="a",
            seq="M" * 10 + "A" * 50 + "C" * 20 + "Z" * 10,
            ref_name="testref",
            group_ref="testref",
            ref_seq=ref_seq,
            match_fraction=0.5,
            reads_count=None,
        ),
        GenotypedContig(
            name="b",
            seq="J" * 10 + "A" * 20 + "C" * 50 + "N" * 10,
            ref_name="testref",
            group_ref="testref",
            ref_seq=ref_seq,
            match_fraction=0.5,
            reads_count=None,
        ),
    ]

    results = list(stitch_contigs(contigs))
    assert len(results) == 1

    result = results[0]

    assert 120 == len(result.seq)
    assert result.seq == "M" * 10 + "A" * 50 + "C" * 50 + "N" * 10

    assert len(visualizer().elements) > len(contigs)


def test_correct_stitching_of_two_partially_overlapping_contigs_real_hiv(
    projects, visualizer
):
    # Scenario: Two partially overlapping contigs are stitched
    # correctly into a single sequence. Not using exact aligner this time.

    ref_name = "HIV1-B-ZA-KP109515-seed"
    ref = projects.getReference(ref_name)

    contigs = [
        GenotypedContig(
            name="a",
            seq=ref[1700:2000],
            ref_name=ref_name,
            group_ref=ref_name,
            ref_seq=ref,
            match_fraction=0.5,
            reads_count=None,
        ),
        GenotypedContig(
            name="b",
            seq=ref[1900:2200],
            ref_name=ref_name,
            group_ref=ref_name,
            ref_seq=ref,
            match_fraction=0.5,
            reads_count=None,
        ),
    ]

    results = list(stitch_contigs(contigs))
    assert len(results) == 1

    result = results[0]

    assert 500 == len(result.seq)
    assert result.seq == ref[1700:2200]

    assert len(visualizer().elements) > len(contigs)


def test_correct_stitching_of_two_partially_overlapping_contigs_real_hcv(
    projects, visualizer
):
    # Scenario: Two partially overlapping contigs are stitched
    # correctly into a single sequence. Not using exact aligner this time.

    ref_name = "HCV-1a"
    ref = projects.getReference(ref_name)
    group_ref = ref_name

    contigs = [
        GenotypedContig(
            name="a",
            seq=ref[1700:2000],
            ref_name=ref_name,
            group_ref=group_ref,
            ref_seq=ref,
            match_fraction=0.5,
            reads_count=None,
        ),
        GenotypedContig(
            name="b",
            seq=ref[1900:2200],
            ref_name=ref_name,
            group_ref=group_ref,
            ref_seq=ref,
            match_fraction=0.5,
            reads_count=None,
        ),
    ]

    results = list(stitch_contigs(contigs))
    assert len(results) == 1

    result = results[0]

    assert 500 == len(result.seq)
    assert result.seq == ref[1700:2200]

    assert len(visualizer().elements) > len(contigs)


def test_correct_stitching_of_two_partially_overlapping_contigs_with_insignificant_gaps(
    projects, visualizer
):
    # Scenario: Two partially overlapping contigs are stitched
    # correctly into a single sequence, with insignificant gaps.

    hxb2_name = "HIV1-B-FR-K03455-seed"
    ref = projects.getReference(hxb2_name)
    gap_ref = "".join(c if i % 30 > 2 else "" for i, c in enumerate(ref))

    contigs = [
        GenotypedContig(
            name="a",
            seq=gap_ref[1700:2000],
            ref_name=hxb2_name,
            group_ref=hxb2_name,
            ref_seq=ref,
            match_fraction=0.5,
            reads_count=None,
        ),
        GenotypedContig(
            name="b",
            seq=gap_ref[1900:2200],
            ref_name=hxb2_name,
            group_ref=hxb2_name,
            ref_seq=ref,
            match_fraction=0.5,
            reads_count=None,
        ),
    ]

    results = list(stitch_contigs(contigs))
    assert len(results) == 1

    result = results[0]

    assert 500 == len(result.seq)
    assert result.seq == gap_ref[1700:2200]

    assert len(visualizer().elements) > len(contigs)


def test_correct_processing_of_two_overlapping_and_one_separate_contig(
    exact_aligner, visualizer
):
    # Scenario: Two overlapping contigs are stitched together, the non-overlapping is kept separate.
    # One contig on the right, and two on the left.

    ref_seq = "Z" * 5 + "A" * 100 + "C" * 100 + "T" * 100 + "Y" * 5

    contigs = [
        GenotypedContig(
            name="a",
            seq="M" * 5 + "A" * 50 + "C" * 20 + "J" * 5,
            ref_name="testref",
            group_ref="testref",
            ref_seq=ref_seq,
            match_fraction=0.5,
            reads_count=None,
        ),
        GenotypedContig(
            name="b",
            seq="Q" * 5 + "A" * 20 + "C" * 50 + "I" * 5,
            ref_name="testref",
            group_ref="testref",
            ref_seq=ref_seq,
            match_fraction=0.5,
            reads_count=None,
        ),
        GenotypedContig(
            name="c",
            seq="N" * 5 + "C" * 20 + "T" * 50 + "H" * 5,
            ref_name="testref",
            group_ref="testref",
            ref_seq=ref_seq,
            match_fraction=0.5,
            reads_count=None,
        ),
    ]

    results = list(stitch_consensus(contigs))
    assert len(results) == 1
    assert results[0].seq == contigs[0].seq.rstrip("J") + "C" * 30 + contigs[
        2
    ].seq.lstrip("N")
    assert len(visualizer().elements) > len(contigs)


def test_correct_processing_of_two_overlapping_and_one_separate_contig_2(
    exact_aligner, visualizer
):
    # Scenario: Two overlapping contigs are stitched together, the non-overlapping is kept separate.
    # One contig on the left, and two on the right.

    ref_seq = "Z" * 5 + "A" * 100 + "C" * 100 + "T" * 100 + "Y" * 5

    contigs = [
        GenotypedContig(
            name="a",
            seq="N" * 5 + "A" * 50 + "C" * 20 + "H" * 5,
            ref_name="testref",
            group_ref="testref",
            ref_seq=ref_seq,
            match_fraction=0.5,
            reads_count=None,
        ),
        GenotypedContig(
            name="b",
            seq="M" * 5 + "C" * 50 + "T" * 20 + "J" * 5,
            ref_name="testref",
            group_ref="testref",
            ref_seq=ref_seq,
            match_fraction=0.5,
            reads_count=None,
        ),
        GenotypedContig(
            name="c",
            seq="Q" * 5 + "C" * 20 + "T" * 50 + "I" * 5,
            ref_name="testref",
            group_ref="testref",
            ref_seq=ref_seq,
            match_fraction=0.5,
            reads_count=None,
        ),
    ]

    results = list(stitch_consensus(contigs))
    assert len(results) == 1
    assert results[0].seq == contigs[0].seq.rstrip("H") + "C" * 30 + contigs[
        2
    ].seq.lstrip("Q")
    assert len(visualizer().elements) > len(contigs)


def test_stitching_of_all_overlapping_contigs_into_one_sequence(
    exact_aligner, visualizer
):
    # Scenario: All contigs have some overlapping parts, resulting in one continuous sequence after stitching.

    ref_seq = "A" * 100 + "C" * 100 + "T" * 100

    contigs = [
        GenotypedContig(
            name="a",
            seq="A" * 50 + "C" * 20,
            ref_name="testref",
            group_ref="testref",
            ref_seq=ref_seq,
            match_fraction=0.5,
            reads_count=None,
        ),
        GenotypedContig(
            name="b",
            seq="A" * 20 + "C" * 100 + "T" * 20,
            ref_name="testref",
            group_ref="testref",
            ref_seq=ref_seq,
            match_fraction=0.5,
            reads_count=None,
        ),
        GenotypedContig(
            name="c",
            seq="C" * 20 + "T" * 50,
            ref_name="testref",
            group_ref="testref",
            ref_seq=ref_seq,
            match_fraction=0.5,
            reads_count=None,
        ),
    ]

    results = list(stitch_contigs(contigs))
    assert len(results) == 1

    result = results[0]

    assert 200 == len(result.seq)
    assert result.seq == "A" * 50 + "C" * 100 + "T" * 50

    assert len(visualizer().elements) > len(contigs)


def test_stitching_with_empty_contigs(exact_aligner, visualizer):
    # Scenario: The function is able to handle and ignore empty contigs.

    ref_seq = "A" * 100

    contigs = [
        GenotypedContig(
            name="a",
            seq=ref_seq,
            ref_name="testref",
            group_ref="testref",
            ref_seq=ref_seq,
            match_fraction=0.5,
            reads_count=None,
        ),
        GenotypedContig(
            name="b",
            seq="",
            ref_name="testref",
            group_ref="testref",
            ref_seq=ref_seq,
            match_fraction=0.5,
            reads_count=None,
        ),
    ]

    results = list(stitch_contigs(contigs))
    assert {contig.seq for contig in contigs} == {contig.seq for contig in results}

    assert len(visualizer().elements) > len(contigs)


def test_stitching_of_identical_contigs(exact_aligner, visualizer):
    # Scenario: The function correctly handles and avoids duplication when identical contigs are stitched together.

    contigs = [
        GenotypedContig(
            name=name,
            seq="ACTGACTG" * 100,
            ref_name="testref",
            group_ref="testref",
            ref_seq="ACTGACTG" * 100,
            match_fraction=1.0,
            reads_count=reads_count,
        )
        for name, reads_count in [("a", 10), ("b", 20), ("c", 100)]
    ]

    results = list(stitch_contigs(contigs))
    assert len(results) == 1
    assert results[0].seq == contigs[2].seq

    assert len(visualizer().elements) > len(contigs)


def test_stitching_of_completely_identical_contigs(exact_aligner, visualizer):
    # Scenario: The function correctly handles and avoids duplication when completely identical contigs
    # are stitched together.

    contigs = [
        GenotypedContig(
            name="x",
            seq="ACTGACTG" * 100,
            ref_name="testref",
            group_ref="testref",
            ref_seq="ACTGACTG" * 100,
            match_fraction=1.0,
            reads_count=copy * 10,  # Different read counts: 10, 20, 30
        )
        for copy in [1, 2, 3]
    ]

    results = list(stitch_contigs(contigs))
    assert len(results) == 1
    assert results[0].seq == contigs[2].seq

    assert len(visualizer().elements) > len(contigs)


def test_stitching_of_zero_contigs(exact_aligner, visualizer):
    # Scenario: The function does not crash if no contigs given.

    contigs = []
    results = list(stitch_contigs(contigs))
    assert results == contigs

    assert len(visualizer().elements) > 0


def test_correct_stitching_of_two_partially_overlapping_different_organism_contigs(
    exact_aligner, visualizer
):
    # Scenario: Two partially overlapping contigs, but which come from different organism,
    # are not stitched into a single sequence.

    ref_seq = "A" * 100 + "C" * 100

    contigs = [
        GenotypedContig(
            name="a",
            seq="A" * 50 + "C" * 20,
            ref_name="testref-1",
            group_ref="testref-1",
            ref_seq=ref_seq,
            match_fraction=0.5,
            reads_count=None,
        ),
        GenotypedContig(
            name="b",
            seq="A" * 20 + "C" * 50,
            ref_name="testref-2",
            group_ref="testref-2",
            ref_seq=ref_seq,
            match_fraction=0.5,
            reads_count=None,
        ),
    ]

    results = list(stitch_contigs(contigs))
    assert len(results) == 2

    assert {contig.seq for contig in contigs} == {contig.seq for contig in results}

    assert len(visualizer().elements) > len(contigs)


def test_correct_processing_complex_nogaps(exact_aligner, visualizer):
    # Scenario: There are two reference organisms.
    # Each with 4 contigs.
    # For each, three overlapping contigs are stitched together, the non-overlapping is kept separate.

    ref_seq = "A" * 100 + "C" * 100 + "T" * 100 + "G" * 100

    contigs = [
        [
            GenotypedContig(
                name="a" + ref_name,
                seq="A" * 50 + "C" * 20,
                ref_name=ref_name,
                group_ref=ref_name,
                ref_seq=ref_seq,
                match_fraction=0.5,
                reads_count=None,
            ),
            GenotypedContig(
                name="b" + ref_name,
                seq="A" * 20 + "C" * 50,
                ref_name=ref_name,
                group_ref=ref_name,
                ref_seq=ref_seq,
                match_fraction=0.5,
                reads_count=None,
            ),
            GenotypedContig(
                name="c" + ref_name,
                seq="C" * 70 + "T" * 20,
                ref_name=ref_name,
                group_ref=ref_name,
                ref_seq=ref_seq,
                match_fraction=0.5,
                reads_count=None,
            ),
            GenotypedContig(
                name="d" + ref_name,
                seq="T" * 20 + "G" * 50,
                ref_name=ref_name,
                group_ref=ref_name,
                ref_seq=ref_seq,
                match_fraction=0.5,
                reads_count=None,
            ),
        ]
        for ref_name in ["testref-1", "testref-2"]
    ]

    contigs = sum(contigs, start=[])

    results = list(stitch_contigs(contigs))
    assert len(results) == 4

    assert 170 == len(results[0].seq)
    assert results[0].seq == "A" * 50 + "C" * 100 + "T" * 20
    assert results[0].group_ref == "testref-1"

    assert 170 == len(results[1].seq)
    assert results[1].seq == "A" * 50 + "C" * 100 + "T" * 20
    assert results[1].group_ref == "testref-2"

    assert results[2].seq == contigs[3].seq
    assert results[3].seq == contigs[7].seq

    assert len(visualizer().elements) > len(contigs)


def test_stitching_when_one_contig_completely_covered_by_another(
    exact_aligner, visualizer
):
    # Scenario: If one contig is completely covered by another contig,
    # the completely covered contig must be dropped.

    ref_seq = "A" * 100 + "C" * 100

    contigs = [
        GenotypedContig(
            name="a",
            seq="M" * 10 + "A" * 20 + "C" * 20 + "O" * 10,
            ref_name="testref",
            group_ref="testref",
            ref_seq=ref_seq,
            match_fraction=0.5,
            reads_count=None,
        ),
        GenotypedContig(
            name="b",
            seq="P" * 10 + "A" * 50 + "C" * 50 + "Z" * 10,
            ref_name="testref",
            group_ref="testref",
            ref_seq=ref_seq,
            match_fraction=0.5,
            reads_count=None,
        ),
    ]

    results = list(stitch_contigs(contigs))
    assert len(results) == 1

    # Test to ensure that the final result contains the contig 'b' and
    # does not contain the completely covered contig 'a'.
    assert results[0].seq == contigs[1].seq

    assert len(visualizer().elements) > len(contigs)


def test_stitching_when_multiple_contigs_completely_covered_by_other_contigs(
    exact_aligner, visualizer
):
    # Scenario: If two contigs are completely covered by another two contigs.

    ref_seq = "A" * 100 + "B" * 100 + "C" * 100 + "D" * 100

    contigs = [
        GenotypedContig(
            name="a",
            seq="M" * 10 + "A" * 20 + "B" * 100 + "C" * 20 + "O" * 10,
            ref_name="testref",
            group_ref="testref",
            ref_seq=ref_seq,
            match_fraction=0.5,
            reads_count=None,
        ),
        GenotypedContig(
            name="b",
            seq="K" * 10 + "B" * 20 + "C" * 100 + "D" * 20 + "J" * 10,
            ref_name="testref",
            group_ref="testref",
            ref_seq=ref_seq,
            match_fraction=0.5,
            reads_count=None,
        ),
        GenotypedContig(
            name="c",
            seq="I" * 10 + "B" * 60 + "C" * 80 + "P" * 10,
            ref_name="testref",
            group_ref="testref",
            ref_seq=ref_seq,
            match_fraction=0.5,
            reads_count=None,
        ),
        GenotypedContig(
            name="d",
            seq="Z" * 10 + "B" * 80 + "C" * 60 + "F" * 10,
            ref_name="testref",
            group_ref="testref",
            ref_seq=ref_seq,
            match_fraction=0.5,
            reads_count=None,
        ),
    ]

    results = list(stitch_contigs(contigs))
    assert len(results) == 1

    assert len(visualizer().elements) > len(contigs)


def test_stitching_contig_with_big_noncovered_gap(exact_aligner, visualizer):
    # Scenario: One contig has a big gap, which is however not covered by anything else.

    ref_seq = "A" * 100 + "C" * 100 + "T" * 100

    contigs = [
        GenotypedContig(
            name="a",
            seq="A" * 50 + "T" * 50,  # mind the C gap
            ref_name="testref",
            group_ref="testref",
            ref_seq=ref_seq,
            match_fraction=0.5,
            reads_count=None,
        ),
    ]

    results = list(stitch_contigs(contigs))

    assert {contig.seq for contig in contigs} == {contig.seq for contig in results}

    assert len(visualizer().elements) > len(contigs)


def test_stitching_contig_with_big_noncovered_gap_2(exact_aligner, visualizer):
    # Scenario: One contig has a big gap, which is however not covered by anything else.

    ref_seq = "A" * 100 + "C" * 100 + "T" * 100 + "G" * 100

    contigs = [
        GenotypedContig(
            name="a",
            seq="A" * 50 + "T" * 50,  # mind the C gap
            ref_name="testref",
            group_ref="testref",
            ref_seq=ref_seq,
            match_fraction=0.5,
            reads_count=None,
        ),
        GenotypedContig(
            name="B",
            seq="G" * 50,
            ref_name="testref",
            group_ref="testref",
            ref_seq=ref_seq,
            match_fraction=0.5,
            reads_count=None,
        ),
    ]

    results = list(stitch_contigs(contigs))

    assert {contig.seq for contig in contigs} == {contig.seq for contig in results}

    assert len(visualizer().elements) > len(contigs)


def test_stitching_contig_with_big_covered_gap(exact_aligner, visualizer):
    # Scenario: If one contig has a big gap covered by another contig.

    ref_seq = "G" * 100 + "A" * 100 + "C" * 100 + "T" * 100 + "G" * 100

    contigs = [
        GenotypedContig(
            name="a",
            seq="G" * 50 + "A" * 50 + "T" * 100,  # mind the gap
            ref_name="testref",
            group_ref="testref",
            ref_seq=ref_seq,
            match_fraction=0.5,
            reads_count=None,
        ),
        GenotypedContig(
            name="b",
            seq="A" * 100 + "C" * 100 + "T" * 50,
            ref_name="testref",
            group_ref="testref",
            ref_seq=ref_seq,
            match_fraction=0.5,
            reads_count=None,
        ),
    ]

    contigs = list(align_all_to_reference(contigs))
    assert len(contigs) == 2
    assert len(list(contigs[0].alignment.deletions())) == 1
    assert len(list(contigs[1].alignment.deletions())) == 0

    results = list(split_contigs_with_gaps(contigs))
    assert len(results) == 3
    assert all(list(contig.alignment.deletions()) == [] for contig in results)

    assert len(visualizer().elements) > len(contigs)


def test_stitching_contig_with_small_covered_gap(exact_aligner, visualizer):
    # Scenario: If one contig has a small gap covered by another contig.

    ref_seq = "G" * 100 + "A" * 29 + "C" * 100 + "T" * 100

    contigs = [
        GenotypedContig(
            name="a",
            seq="G" * 100 + "A" * 0 + "C" * 100,  # mind the gap
            ref_name="testref",
            group_ref="testref",
            ref_seq=ref_seq,
            match_fraction=0.5,
            reads_count=None,
        ),
        GenotypedContig(
            name="b",
            seq="A" * 29 + "C" * 50,
            ref_name="testref",
            group_ref="testref",
            ref_seq=ref_seq,
            match_fraction=0.5,
            reads_count=None,
        ),
    ]

    contigs = list(align_all_to_reference(contigs))
    assert len(contigs) == 2
    assert len(list(contigs[0].alignment.deletions())) == 1
    assert len(list(contigs[1].alignment.deletions())) == 0
    results = list(split_contigs_with_gaps(contigs))
    assert len(results) == 3

    assert len(visualizer().elements) > len(contigs)

    assert all(x.seq == lstrip(rstrip(x)).seq for x in results)
    assert {contig.seq for contig in contigs} != {contig.seq for contig in results}


def test_stitching_partial_align(exact_aligner, visualizer):
    # Scenario: A single contig has a sequence that partially aligns to the reference sequence.

    contigs = [
        GenotypedContig(
            name="a",
            seq="T" * 10 + "C" * 20 + "A" * 10,
            ref_name="testref",
            group_ref="testref",
            ref_seq="A" * 20 + "C" * 20 + "T" * 20,
            match_fraction=0.3,
            reads_count=None,
        ),
    ]

    results = list(stitch_contigs(contigs))
    assert len(results) == len(contigs)
    for result in results:
        assert any(result.seq in contig.seq for contig in contigs)

    assert len(visualizer().elements) > len(contigs)

    assert all(x.seq != lstrip(rstrip(x)).seq for x in results)

    assert {contig.seq for contig in contigs} != {
        lstrip(rstrip(contig)).seq for contig in results
    }


def test_partial_align_consensus(exact_aligner, visualizer):
    # Scenario: A single contig partially aligns to the reference sequence, and a consensus sequence is being stitched.

    contigs = [
        GenotypedContig(
            name="a",
            seq="T" * 10 + "C" * 20 + "A" * 10,
            ref_name="testref",
            group_ref="testref",
            ref_seq="A" * 20 + "C" * 20 + "T" * 20,
            match_fraction=0.3,
            reads_count=None,
        ),
    ]

    results = list(stitch_consensus(contigs))
    assert len(results) == len(contigs)
    assert {contig.seq for contig in contigs} == {contig.seq for contig in results}

    assert len(visualizer().elements) > len(contigs)


def test_stitching_partial_align_multiple_sequences(exact_aligner, visualizer):
    # Scenario: Multiple contigs have sequences that partially align to the same reference sequence.

    ref_seq = "A" * 20 + "C" * 20 + "T" * 20

    contigs = [
        GenotypedContig(
            name="a",
            seq="Z" * 5 + "C" * 20 + "T" * 5 + "U" * 5,
            ref_name="testref",
            group_ref="testref",
            ref_seq=ref_seq,
            match_fraction=0.3,
            reads_count=None,
        ),
        GenotypedContig(
            name="b",
            seq="M" * 5 + "C" * 5 + "T" * 10 + "G" * 10,
            ref_name="testref",
            group_ref="testref",
            ref_seq=ref_seq,
            match_fraction=0.3,
            reads_count=None,
        ),
    ]

    results = list(stitch_contigs(contigs))
    assert len(results) == 1
    assert results[0].seq == "Z" * 5 + "C" * 20 + "T" * 10 + "G" * 10
    assert len(visualizer().elements) > len(contigs)


def test_partial_align_consensus_multiple_sequences(exact_aligner, visualizer):
    # Scenario: Multiple contigs partially align to the same reference sequence,
    # and a consensus sequence is being stitched from them.

    ref_seq = "A" * 20 + "C" * 20 + "T" * 20

    contigs = [
        GenotypedContig(
            name="a",
            seq="A" * 20,
            ref_name="testref",
            group_ref="testref",
            ref_seq=ref_seq,
            match_fraction=0.3,
            reads_count=None,
        ),
        GenotypedContig(
            name="b",
            seq="T" * 20,
            ref_name="testref",
            group_ref="testref",
            ref_seq=ref_seq,
            match_fraction=0.3,
            reads_count=None,
        ),
    ]

    results = list(stitch_consensus(contigs))
    assert len(results) == 1
    assert results[0].seq == contigs[0].seq + contigs[1].seq

    assert len(visualizer().elements) > len(contigs)


def test_partial_align_consensus_multiple_overlaping_sequences(
    exact_aligner, visualizer
):
    # Scenario: Multiple contigs partially align to the same reference sequence,
    # and a consensus sequence is being stitched from them.

    ref_seq = "A" * 20 + "C" * 20 + "T" * 20

    contigs = [
        GenotypedContig(
            name="a",
            seq="T" * 10 + "A" * 5 + "C" * 20 + "A" * 10,
            ref_name="testref",
            group_ref="testref",
            ref_seq=ref_seq,
            match_fraction=0.3,
            reads_count=None,
        ),
        GenotypedContig(
            name="b",
            seq="C" * 20 + "T" * 5 + "A" * 10 + "G" * 10,
            ref_name="testref",
            group_ref="testref",
            ref_seq=ref_seq,
            match_fraction=0.3,
            reads_count=None,
        ),
    ]

    results = list(stitch_consensus(contigs))
    assert len(results) == 1
    assert (
        results[0].seq == "T" * 10 + "A" * 5 + "C" * 20 + "T" * 5 + "A" * 10 + "G" * 10
    )
    assert results[0].seq == contigs[0].seq[:-10] + contigs[1].seq[20:]

    assert len(visualizer().elements) > len(contigs)


def test_big_insertion_in_a_single_contig(projects, visualizer):
    # Scenario: Single contig produces many alignments.

    hxb2_name = "HIV1-B-FR-K03455-seed"
    ref_seq = projects.getReference(hxb2_name)
    seq = ref_seq[2000:3000] + "C" * 300 + ref_seq[3100:4000]

    contigs = [
        GenotypedContig(
            name="a",
            seq=seq,
            ref_name="testref",
            group_ref="testref",
            ref_seq=ref_seq,
            match_fraction=0.3,
            reads_count=None,
        ),
    ]

    results = list(stitch_consensus(contigs))
    assert len(results) == 1
    assert results[0].seq == contigs[0].seq

    assert len(visualizer().elements) > len(contigs)


def test_big_insertion_in_a_single_contig_2(exact_aligner, visualizer):
    # Scenario: Single contig produces many alignments.

    ref_seq = "A" * 10 + "B" * 20 + "C" * 10

    contigs = [
        GenotypedContig(
            name="a",
            seq="A" * 10 + "D" * 100 + "C" * 10,
            ref_name="testref",
            group_ref="testref",
            ref_seq=ref_seq,
            match_fraction=0.3,
            reads_count=None,
        ),
    ]

    results = list(stitch_consensus(contigs))
    assert len(results) == 1
    assert results[0].seq == contigs[0].seq

    assert len(visualizer().elements) > len(contigs)


def test_gap_around_small_insertion(exact_aligner, visualizer):
    # Scenario: Contig is split around its gap, then stripped.

    ref_seq = "A" * 10 + "B" * 29 + "C" * 10

    contigs = [
        GenotypedContig(
            name="a",
            seq="P" * 5 + "A" * 10 + "D" * 6 + "C" * 10 + "Z" * 5,
            ref_name="testref",
            group_ref="testref",
            ref_seq=ref_seq,
            match_fraction=0.3,
            reads_count=None,
        ),
        GenotypedContig(
            name="b",
            seq="Q" * 5 + "B" * 29 + "J" * 5,
            ref_name="testref",
            group_ref="testref",
            ref_seq=ref_seq,
            match_fraction=0.3,
            reads_count=None,
        ),
    ]

    results = list(stitch_consensus(contigs))
    assert len(results) == 1
    assert results[0].seq == "P" * 5 + "A" * 10 + "B" * 29 + "C" * 10 + "Z" * 5
    assert len(visualizer().elements) > len(contigs)


def test_gap_around_big_insertion(exact_aligner, visualizer):
    # Scenario: Contig is split around its gap, then stripped.

    ref_seq = "A" * 10 + "B" * 29 + "C" * 10

    contigs = [
        GenotypedContig(
            name="a",
            seq="P" * 5 + "A" * 10 + "D" * 100 + "C" * 10 + "Z" * 5,
            ref_name="testref",
            group_ref="testref",
            ref_seq=ref_seq,
            match_fraction=0.3,
            reads_count=None,
        ),
        GenotypedContig(
            name="b",
            seq="Q" * 5 + "B" * 29 + "J" * 5,
            ref_name="testref",
            group_ref="testref",
            ref_seq=ref_seq,
            match_fraction=0.3,
            reads_count=None,
        ),
    ]

    results = list(stitch_consensus(contigs))
    assert len(results) == 1
    assert results[0].seq == "P" * 5 + "A" * 10 + "B" * 29 + "C" * 10 + "Z" * 5
    assert len(visualizer().elements) > len(contigs)


def test_stitch_with_insertion(exact_aligner, visualizer):
    # Scenario: Contig is aligned with multiple hits, and the borders are correctly handled.

    ref_seq = "X" * 5 + "A" * 10 + "B" * 20 + "C" * 10 + "M" * 5

    contigs = [
        GenotypedContig(
            name="a",
            seq="P" * 5 + "A" * 10 + "D" * 6 + "C" * 10 + "Z" * 5,
            ref_name="testref",
            group_ref="testref",
            ref_seq=ref_seq,
            match_fraction=0.3,
            reads_count=None,
        ),
    ]

    results = list(stitch_consensus(contigs))
    assert len(results) == 1
    assert results[0].seq == "PPPPPAAAAAAAAAADDDDDDCCCCCCCCCCZZZZZ"
    assert len(visualizer().elements) > len(contigs)


def test_stitch_cross_alignment(exact_aligner, visualizer):
    # Scenario: Single contig is cross-aligned.

    ref_seq = "X" * 5 + "A" * 10 + "B" * 20 + "C" * 10 + "M" * 5

    contigs = [
        GenotypedContig(
            name="a",
            seq="P" * 5 + "C" * 10 + "D" * 6 + "A" * 10 + "Z" * 5,
            ref_name="testref",
            group_ref="testref",
            ref_seq=ref_seq,
            match_fraction=0.3,
            reads_count=None,
        ),
    ]

    results = list(stitch_consensus(contigs))
    assert len(results) == 1
    assert results[0].seq == "AAAAAAAAAACCCCCCCCCC"
    assert len(visualizer().elements) > len(contigs)


def test_cross_alignment_around_small_insertion(exact_aligner, visualizer):
    # Scenario: Single contig is cross-aligned, then combined with another contig that is between its aligned parts.

    ref_seq = "X" * 5 + "A" * 10 + "B" * 20 + "C" * 10 + "M" * 5

    contigs = [
        GenotypedContig(
            name="a",
            seq="P" * 5 + "C" * 10 + "D" * 6 + "A" * 10 + "Z" * 5,
            ref_name="testref",
            group_ref="testref",
            ref_seq=ref_seq,
            match_fraction=0.3,
            reads_count=None,
        ),
        GenotypedContig(
            name="b",
            seq="Q" * 5 + "B" * 20 + "J" * 5,
            ref_name="testref",
            group_ref="testref",
            ref_seq=ref_seq,
            match_fraction=0.3,
            reads_count=None,
        ),
    ]

    results = list(stitch_consensus(contigs))
    assert len(results) == 1
    assert results[0].seq == "A" * 10 + "B" * 20 + "C" * 10
    assert len(visualizer().elements) > len(contigs)


def test_reverse_complement_match(projects, visualizer):
    # Scenario: Single contig is aligned in the reverse strand.

    from mappy import revcomp

    hxb2_name = "HIV1-B-FR-K03455-seed"
    ref = projects.getReference(hxb2_name)
    ref_part = ref[2000:2200]
    seq = revcomp(ref_part)

    contigs = [
        GenotypedContig(
            name="a",
            seq=seq,
            ref_name="testref",
            group_ref="testref",
            ref_seq=ref,
            match_fraction=0.3,
            reads_count=None,
        ),
    ]

    results = list(stitch_consensus(contigs))
    assert len(results) == 1
    assert results[0].seq == ref_part
    assert len(visualizer().elements) > len(contigs)


def test_reverse_complement_match_with_padding(projects, visualizer):
    # Scenario: Single contig is aligned in the reverse strand.

    from mappy import revcomp

    hxb2_name = "HIV1-B-FR-K03455-seed"
    ref = projects.getReference(hxb2_name)
    ref_part = "T" * 24 + ref[2000:2200] + "G" * 27
    seq = revcomp(ref_part)

    contigs = [
        GenotypedContig(
            name="a",
            seq=seq,
            ref_name="testref",
            group_ref="testref",
            ref_seq=ref,
            match_fraction=0.3,
            reads_count=None,
        ),
    ]

    results = list(stitch_consensus(contigs))
    assert len(results) == 1
    assert results[0].seq == ref_part
    assert len(lstrip(results[0]).seq) == len(ref_part) - 24
    assert len(rstrip(results[0]).seq) == len(ref_part) - 27
    assert rstrip(results[0]).seq == ref_part[:-27]  # 27 Gs on the right
    assert lstrip(results[0]).seq == ref_part[24:]  # 24 Ts on the left
    assert len(visualizer().elements) > len(contigs)


def test_multiple_reverse_complement_matches(projects, visualizer):
    # Scenario: Single contig is aligned in the reverse strand in multiple places.

    from mappy import revcomp

    hxb2_name = "HIV1-B-FR-K03455-seed"
    ref = projects.getReference(hxb2_name)
    ref_part = (
        "T" * 24
        + ref[2000:2600]
        + "A" * 9
        + ref[3000:3600]
        + "T" * 9
        + ref[4000:4600]
        + "G" * 27
    )
    seq = revcomp(ref_part)

    contigs = [
        GenotypedContig(
            name="a",
            seq=seq,
            ref_name="testref",
            group_ref="testref",
            ref_seq=ref,
            match_fraction=0.3,
            reads_count=None,
        ),
    ]

    results = list(stitch_consensus(contigs))
    assert len(results) == 1
    assert len(results[0].seq) == len(ref_part)
    assert results[0].seq == ref_part
    assert len(lstrip(results[0]).seq) == len(ref_part) - 24
    assert len(rstrip(results[0]).seq) == len(ref_part) - 27
    assert lstrip(results[0]).seq == ref_part[24:]
    assert rstrip(results[0]).seq == ref_part[:-27]

    assert len(visualizer().elements) > len(contigs)


def test_multiple_reverse_complement_matches_out_of_order(projects, visualizer):
    # Scenario: Single contig is aligned in the reverse strand in multiple places, producing an out of order alignment.

    from mappy import revcomp

    hxb2_name = "HIV1-B-FR-K03455-seed"
    ref = projects.getReference(hxb2_name)
    ref_part = (
        "T" * 24
        + ref[2000:2600]
        + "A" * 9
        + ref[3000:3600]
        + "T" * 9
        + ref[4000:4600]
        + "G" * 27
    )
    seq = revcomp(ref_part)

    contigs = [
        GenotypedContig(
            name="a",
            seq=seq,
            ref_name="testref",
            group_ref="testref",
            ref_seq=ref,
            match_fraction=0.3,
            reads_count=None,
        ),
    ]

    results = list(stitch_consensus(contigs))
    assert len(results) == 1
    assert len(results[0].seq) == len(ref_part)
    assert results[0].seq == ref_part
    assert len(lstrip(results[0]).seq) == len(ref_part) - 24
    assert len(rstrip(results[0]).seq) == len(ref_part) - 27
    assert lstrip(results[0]).seq == ref_part[24:]
    assert rstrip(results[0]).seq == ref_part[:-27]
    assert len(visualizer().elements) > len(contigs)


def test_forward_and_reverse_match(projects, visualizer):
    # Scenario: Single contig is aligned in both strands.

    from mappy import revcomp

    hxb2_name = "HIV1-B-FR-K03455-seed"
    ref = projects.getReference(hxb2_name)
    seq = ref[1000:1100] + revcomp(ref[2000:2200])

    contigs = [
        GenotypedContig(
            name="a",
            seq=seq,
            ref_name="testref",
            group_ref="testref",
            ref_seq=ref,
            match_fraction=0.3,
            reads_count=None,
        ),
    ]

    results = list(stitch_consensus(contigs))
    assert len(results) == 1
    assert results[0].seq == seq
    assert len(visualizer().elements) > len(contigs)


def test_overlaping_in_reference_space(projects, visualizer, monkeypatch):
    # Scenario: Single contig is aligned in two parts that overlap in reference space.

    def mock_align(
        reference_seq: str, consensus: str
    ) -> Tuple[List[MockAlignment], str]:
        alignments = [
            MockAlignment(
                ctg="N/A",
                ctg_len=0,
                strand=1,
                mapq=60,
                is_primary=True,
                q_st=100,
                q_en=300,
                r_st=200,
                r_en=400,
                cigar=[(200, CigarActions.MATCH)],
                cigar_str="200M",
            ),
            MockAlignment(
                ctg="N/A",
                ctg_len=0,
                strand=1,
                mapq=60,
                is_primary=True,
                q_st=300,
                q_en=500,
                r_st=300,
                r_en=500,
                cigar=[(200, CigarActions.MATCH)],
                cigar_str="200M",
            ),
        ]
        algorithm = "mock"
        return (alignments, algorithm)

    monkeypatch.setattr("micall.utils.referencefull_contig_stitcher.align_consensus", mock_align)

    ref = "A" * 700
    seq = "C" * 600

    contigs = [
        GenotypedContig(
            name="a",
            seq=seq,
            ref_name="testref",
            group_ref="testref",
            ref_seq=ref,
            match_fraction=0.3,
            reads_count=None,
        ),
    ]

    results = list(stitch_consensus(contigs))
    assert len(results) == 1
    assert results[0].seq == "C" * 500

    assert isinstance(results[0], stitcher.AlignedContig)
    assert results[0].alignment == CigarHit(
        Cigar.parse("300M"), r_st=200, r_ei=499, q_st=100, q_ei=399
    )

    assert len(visualizer().elements) > len(contigs)


def test_correct_stitching_of_one_normal_and_one_unknown(exact_aligner, visualizer):
    # Scenario: Two partially overlapping contigs are stitched correctly into a single sequence.

    ref_seq = "A" * 100 + "C" * 100

    contigs = [
        GenotypedContig(
            name="a",
            seq="A" * 50 + "C" * 20,
            ref_name="testref",
            group_ref="testref",
            ref_seq=ref_seq,
            match_fraction=0.5,
            reads_count=None,
        ),
        GenotypedContig(
            name="b",
            seq="A" * 20 + "C" * 50,
            ref_name=None,
            group_ref=None,
            ref_seq=None,
            match_fraction=0.5,
            reads_count=None,
        ),
    ]

    results = list(stitch_contigs(contigs))
    assert len(results) == 2

    assert 70 == len(results[0].seq)
    assert 70 == len(results[1].seq)

    assert {result.seq for result in results} == {contig.seq for contig in contigs}

    assert len(visualizer().elements) > len(contigs)


def test_main_invocation(exact_aligner, tmp_path, hcv_db):
    from micall.core.contig_stitcher import main
    pwd = os.path.dirname(__file__)
    contigs = os.path.join(pwd, "data", "exact_parts_contigs.csv")
    remap_counts = os.path.join(pwd, "data", "exact_parts_contigs_remap_counts.csv")
    stitched_contigs = os.path.join(tmp_path, "stitched.csv")
    main(['with-references', contigs, stitched_contigs, '--remap-counts', remap_counts])

    assert os.path.exists(contigs)
    assert os.path.exists(stitched_contigs)

    # Check the contents of stitched_contigs
    with open(stitched_contigs, "r") as stitched_file:
        stitched_data = stitched_file.read()

    expected_file_path = os.path.join(pwd, "data", "exact_parts_contigs_stitched.csv")
    with open(expected_file_path, "r") as expected_file:
        expected_data = expected_file.read()

    assert stitched_data == expected_data, (
        "The contents of the stitched contigs file do not match the expected contents."
    )


def test_visualizer_simple(exact_aligner, tmp_path, hcv_db):
    from micall.core.contig_stitcher import main
    pwd = os.path.dirname(__file__)
    contigs = os.path.join(pwd, "data", "exact_parts_contigs.csv")
    remap_counts = os.path.join(pwd, "data", "exact_parts_contigs_remap_counts.csv")
    stitched_contigs = os.path.join(tmp_path, "stitched.csv")
    plot = os.path.join(tmp_path, "exact_parts_contigs.plot.svg")
    main(['with-references', contigs, stitched_contigs, "--remap-counts", remap_counts, "--debug", "--plot", plot])

    assert os.path.exists(contigs)
    assert os.path.exists(stitched_contigs)

    # Check the contents of stitched_contigs
    with open(stitched_contigs, "r") as stitched_file:
        stitched_data = stitched_file.read()

    expected_file_path = os.path.join(pwd, "data", "exact_parts_contigs_stitched.csv")
    with open(expected_file_path, "r") as expected_file:
        expected_data = expected_file.read()
        assert stitched_data == expected_data, (
            "The contents of the stitched contigs file do not match the expected contents."
        )

    # Check the contents of stitched_contigs
    expected_plot = os.path.join(pwd, "data", "exact_parts_contigs.plot.svg")
    with open(plot, "r") as stitched_file, open(expected_plot, "r") as expected_file:
        stitched_data = stitched_file.read()
        expected_data = expected_file.read()
        assert stitched_data == expected_data, (
            "The contents of the stitched plot file do not match the expected contents."
        )


def test_main_invocation_without_remap_counts(exact_aligner, tmp_path, hcv_db):
    """Test that main() works without remap_counts argument (all reads_count initialized to None)."""
    from micall.core.contig_stitcher import main
    pwd = os.path.dirname(__file__)
    contigs = os.path.join(pwd, "data", "exact_parts_contigs.csv")
    stitched_contigs = os.path.join(tmp_path, "stitched.csv")

    # Call main without --remap-counts
    main(['with-references', contigs, stitched_contigs])

    assert os.path.exists(contigs)
    assert os.path.exists(stitched_contigs)

    # The file should exist and have content (though results may differ without read counts)
    with open(stitched_contigs, "r") as stitched_file:
        stitched_data = stitched_file.read()

    # Should have header and at least some contigs
    assert "ref,match,group_ref,contig" in stitched_data
    assert len(stitched_data.strip().split('\n')) > 1  # More than just header


def test_visualizer_correct_labeling_of_different_organism_contigs(
    exact_aligner, visualizer
):
    # Scenario: Some discarded and anomaly contigs correctly labelled.

    ref_seq = "A" * 100 + "C" * 100

    contigs = [
        GenotypedContig(
            name="a",
            seq="A" * 50,
            ref_name="testref-1",
            group_ref="testref-1",
            ref_seq=ref_seq,
            match_fraction=0.5,
            reads_count=None,
        ),
        GenotypedContig(
            name="b",
            seq="C" * 50,
            ref_name="testref-2",
            group_ref="testref-2",
            ref_seq=ref_seq,
            match_fraction=0.5,
            reads_count=None,
        ),
        GenotypedContig(
            name="a_anomaly",
            seq="D" * 50,
            ref_name="testref-1",
            group_ref="testref-1",
            ref_seq=ref_seq,
            match_fraction=0.5,
            reads_count=None,
        ),
        GenotypedContig(
            name="b_discarded",
            seq="C" * 20,
            ref_name="testref-2",
            group_ref="testref-2",
            ref_seq=ref_seq,
            match_fraction=0.5,
            reads_count=None,
        ),
        GenotypedContig(
            name="some_anomaly",
            seq="T" * 20,
            ref_name="unknown",
            group_ref=None,
            ref_seq=ref_seq,
            match_fraction=0.5,
            reads_count=None,
        ),
        GenotypedContig(
            name="some_unknown",
            seq="T" * 20,
            ref_name="unknown",
            group_ref=None,
            ref_seq=None,
            match_fraction=0.5,
            reads_count=None,
        ),
    ]

    results = list(stitch_contigs(contigs))
    assert len(results) == 5

    assert len(visualizer().elements) > len(contigs)


#  _   _       _ _     _            _
# | | | |_ __ (_) |_  | |_ ___  ___| |_ ___
# | | | | '_ \| | __| | __/ _ \/ __| __/ __|
# | |_| | | | | | |_  | ||  __/\__ \ |_\__ \
#  \___/|_| |_|_|\__|  \__\___||___/\__|___/
#


@pytest.mark.parametrize(
    "intervals, expected",
    [
        ([], []),
        ([(1, 3)], [(1, 3)]),
        # Non-overlapping intervals
        ([(1, 3), (5, 6)], [(1, 3), (5, 6)]),
        # Directly overlapping intervals
        ([(1, 3), (2, 5)], [(1, 5)]),
        # Adjacent intervals that exactly touch each other
        ([(1, 2), (3, 4)], [(1, 4)]),
        # Nested intervals
        ([(1, 10), (2, 5)], [(1, 10)]),
        # Multiple merged intervals
        ([(1, 3), (2, 4), (6, 8), (10, 11), (11, 12)], [(1, 4), (6, 8), (10, 12)]),
        # Intervals out of initial order
        ([(4, 6), (1, 2)], [(1, 2), (4, 6)]),
        # Overlapping intervals with out of order inputs
        ([(1, 4), (3, 5), (2, 3), (7, 10), (9, 12)], [(1, 5), (7, 12)]),
        # Large set of intervals with various overlaps
        (
            [(1, 4), (2, 6), (5, 8), (7, 8), (10, 15), (11, 12), (13, 14), (17, 18)],
            [(1, 8), (10, 15), (17, 18)],
        ),
        # Intervals where end is less than start should return
        # as is or be handled explicitly depending on implementation
        ([(5, 3), (1, 2)], [(1, 2), (5, 3)]),
        # Intervals that are exactly one after the other in sequence / Intervals that are completely disjoint
        ([(1, 2), (4, 5), (7, 8)], [(1, 2), (4, 5), (7, 8)]),
        # Overlapping intervals that merge into one large interval
        ([(2, 6), (4, 10), (5, 15), (14, 20)], [(2, 20)]),
        # Same interval repeated multiple times
        ([(1, 5), (1, 5), (1, 5)], [(1, 5)]),
        # Single point intervals
        ([(1, 1), (5, 5), (3, 3)], [(1, 1), (3, 3), (5, 5)]),
        ([(1, 1), (5, 5), (3, 3), (1, 1), (1, 1)], [(1, 1), (3, 3), (5, 5)]),
        ([(1, 1), (2, 3)], [(1, 3)]),
        # Intervals that start with negative numbers
        ([(-5, 0), (-2, 3), (1, 7), (9, 12)], [(-5, 7), (9, 12)]),
    ],
)
def test_merge_intervals(intervals, expected):
    assert merge_intervals(intervals) == expected


class MockAlignedContig:
    @dataclass
    class TestMockAlignment:
        r_st: int
        r_ei: int

        @property
        def ref_length(self) -> int:
            """Match aligntools.CigarHit API for tests that rely on ref_length."""

            return abs(self.r_ei - self.r_st)

    def __init__(self, ref_name, group_ref, r_st, r_ei, name="contig", reads_count=None):
        self.ref_name = ref_name
        self.group_ref = group_ref
        self.alignment = MockAlignedContig.TestMockAlignment(r_st, r_ei)
        self.name = name
        self.reads_count = reads_count
        self.id = id(self)
        self.reads_count = reads_count
        self._unique_name = name  # Simple version for testing

    @property
    def unique_name(self):
        return self._unique_name


# Simple function to create mock AlignedContig objects for testing, including ref_name.
def create_mock_aligned_contig(ref_name, r_st, r_ei, name="contig", reads_count=None):
    return MockAlignedContig(ref_name, ref_name, r_st, r_ei, name, reads_count)


@pytest.mark.parametrize(
    "contigs, expected_covered_name",
    [
        # No contigs are completely covered.
        ([("ref1", 0, 100), ("ref1", 101, 200)], None),
        ([("ref1", 0, 50), ("ref1", 51, 100)], None),
        # A single contig is completely covered by one other contig.
        ([("ref1", 0, 100), ("ref1", 0, 200)], "contig1"),
        ([("ref1", 50, 150), ("ref1", 0, 200)], "contig1"),
        # A single contig completely covers another, but with different reference names.
        ([("ref1", 0, 50), ("ref2", 0, 100)], None),
        # Single coverage with exact match - now requires read count comparison
        # Since neither contig has reads_count set (None), total_coverage (0) > current (0) is False
        # So the contig is NOT removed
        ([("ref1", 0, 100), ("ref1", 0, 100)], None),
        # A single contig is completely covered at the beginning by one and at the end by another contig.
        ([("ref1", 0, 50), ("ref1", 50, 100), ("ref1", 25, 75)], "contig3"),
        # Contigs overlap but none are completely covered.
        ([("ref1", 0, 50), ("ref1", 40, 90), ("ref1", 80, 120)], None),
        # Multiple contigs with some covered completely by a single other contig.
        ([("ref1", 0, 200), ("ref1", 10, 30), ("ref1", 170, 190)], "contig2"),
        # Multiple contigs with complex overlaps and one completely covered.
        (
            [("ref1", 30, 60), ("ref1", 0, 50), ("ref1", 20, 70), ("ref1", 60, 90)],
            "contig1",
        ),
        # Edge case where a contig starts where another ends.
        ([("ref1", 0, 50), ("ref1", 50, 100)], None),
        # Contigs are completely covered in a nested fashion.
        ([("ref1", 0, 200), ("ref1", 50, 150), ("ref1", 100, 125)], "contig3"),
        # Contigs are adjacent and cover each other completely.
        ([("ref1", 0, 100), ("ref1", 101, 200), ("ref1", 0, 200)], "contig2"),
        # Single large contig covers several smaller non-adjacent contigs.
        (
            [
                ("ref1", 0, 500),
                ("ref1", 50, 100),
                ("ref1", 200, 250),
                ("ref1", 300, 350),
            ],
            "contig2",
        ),
        # Single large contig covers several smaller adjacent contigs.
        (
            [
                ("ref1", 50, 100),
                ("ref1", 70, 300),
                ("ref1", 101, 199),
                ("ref1", 200, 350),
            ],
            "contig3",
        ),
        # Single small contig is covered by several larger contigs.
        (
            [
                ("ref1", 0, 250),
                ("ref1", 200, 300),
                ("ref1", 600, 800),
                ("ref1", 250, 700),
            ],
            "contig2",
        ),
        # Complex case with multiple contigs and complete coverage by combinations.
        (
            [
                ("ref1", 0, 100),
                ("ref1", 30, 130),
                ("ref1", 60, 160),
                ("ref1", 90, 190),
                ("ref1", 120, 220),
            ],
            "contig2",
        ),
        # Contigs with same start but different end, where one is covered.
        ([("ref1", 0, 100), ("ref1", 0, 50)], "contig2"),
        # Contigs with same end but different start, where one is covered.
        ([("ref1", 50, 100), ("ref1", 0, 100)], "contig1"),
        # Contig covered by two overlapping contigs that don't individually cover the whole range.
        ([("ref1", 0, 75), ("ref1", 25, 100), ("ref1", 0, 100)], "contig1"),
        # Two contigs are covered completely by one large contig.
        ([("ref1", 0, 300), ("ref1", 50, 100), ("ref1", 200, 250)], "contig2"),
        # No contigs at all.
        ([], None),
    ],
)
def test_find_covered(contigs, expected_covered_name):
    stitcher.ReferencefullStitcherContext.set(stitcher.ReferencefullStitcherContext())
    mock_contigs = [
        create_mock_aligned_contig(ref_name, r_st, r_ei, f"contig{i + 1}")
        for i, (ref_name, r_st, r_ei) in enumerate(contigs)
    ]
    covered, covering = find_covered_contig(mock_contigs)
    if expected_covered_name is None:
        assert covered is None
    else:
        assert covered is not None
        assert covered.name == expected_covered_name


def test_find_covered_prioritizes_by_reads_count():
    """Test that find_covered_contig prioritizes by reads_count when available."""
    stitcher.ReferencefullStitcherContext.set(stitcher.ReferencefullStitcherContext())

    # Create three contigs where contig2 and contig3 are both fully covered by contig1:
    # - contig1: ref_length=200, reads_count=100 (large, covers both others)
    # - contig2: ref_length=50,  reads_count=30  (medium ref_length, medium reads)
    # - contig3: ref_length=40,  reads_count=10  (small ref_length, low reads)
    # Both contig2 and contig3 are fully inside contig1's range.

    mock_contigs = [
        create_mock_aligned_contig("ref1", 0, 200, "contig1", reads_count=100),
        create_mock_aligned_contig("ref1", 50, 100, "contig2", reads_count=30),
        create_mock_aligned_contig("ref1", 120, 160, "contig3", reads_count=10),
    ]

    # With reads_count available, should remove contig3 first (lowest reads_count=10)
    # even though contig2 has larger ref_length
    covered, covering = find_covered_contig(mock_contigs)
    assert covered is not None
    assert covered.name == "contig3"


def test_find_covered_prioritizes_by_ref_length_when_no_reads():
    """Test that find_covered_contig falls back to ref_length when reads_count is unavailable."""
    stitcher.ReferencefullStitcherContext.set(stitcher.ReferencefullStitcherContext())

    # Create three contigs where contig2 and contig3 are both fully covered by contig1:
    # - contig1: ref_length=200 (large, covers both others)
    # - contig2: ref_length=50  (medium)
    # - contig3: ref_length=40  (smallest)

    mock_contigs = [
        create_mock_aligned_contig("ref1", 0, 200, "contig1", reads_count=None),
        create_mock_aligned_contig("ref1", 50, 100, "contig2", reads_count=None),
        create_mock_aligned_contig("ref1", 120, 160, "contig3", reads_count=None),
    ]

    # Without reads_count, should remove contig3 first (smallest ref_length=40)
    covered, covering = find_covered_contig(mock_contigs)
    assert covered is not None
    assert covered.name == "contig3"


def test_find_covered_mixed_reads_count_uses_ref_length():
    """Test that find_covered_contig uses ref_length when reads_count is mixed (some None)."""
    stitcher.ReferencefullStitcherContext.set(stitcher.ReferencefullStitcherContext())

    # Create contigs with mixed reads_count (some None, some set):
    # - contig1: ref_length=200, reads_count=100
    # - contig2: ref_length=50,  reads_count=None (medium ref_length)
    # - contig3: ref_length=40,  reads_count=10  (smallest ref_length)

    mock_contigs = [
        create_mock_aligned_contig("ref1", 0, 200, "contig1", reads_count=100),
        create_mock_aligned_contig("ref1", 50, 100, "contig2", reads_count=None),
        create_mock_aligned_contig("ref1", 120, 160, "contig3", reads_count=10),
    ]

    # Mixed reads_count means fall back to ref_length ordering
    # Should remove contig3 first (smallest ref_length=40)
    covered, covering = find_covered_contig(mock_contigs)
    assert covered is not None
    assert covered.name == "contig3"
