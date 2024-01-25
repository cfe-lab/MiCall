
import random
import logging
import os
import pytest

import micall.core.contig_stitcher as stitcher
from micall.core.contig_stitcher import split_contigs_with_gaps, stitch_contigs, GenotypedContig, merge_intervals, find_covered_contig, stitch_consensus, calculate_concordance, align_all_to_reference, main, AlignedContig, disambiguate_concordance
from micall.core.plot_contigs import plot_stitcher_coverage
from micall.tests.utils import MockAligner, fixed_random_seed
from micall.tests.test_denovo import check_hcv_db # activates the fixture


logging.getLogger("micall.core.contig_stitcher").setLevel(logging.DEBUG)
logging.getLogger("micall.core.plot_contigs").setLevel(logging.DEBUG)


@pytest.fixture()
def exact_aligner(monkeypatch):
    monkeypatch.setattr('micall.core.contig_stitcher.Aligner', MockAligner)


@pytest.fixture
def visualizer(request, tmp_path):
    logs = stitcher.context.set(stitcher.StitcherContext())
    test_name = request.node.name
    plot_name = test_name + ".svg"
    pwd = os.path.dirname(__file__)
    plots_dir = os.path.join(pwd, "data", "stitcher_plots")
    os.makedirs(plots_dir, exist_ok=True)
    path_to_expected = os.path.join(plots_dir, plot_name)
    path_to_produced = os.path.join(tmp_path, plot_name)

    def check():
        logs = stitcher.context.get().events
        figure = plot_stitcher_coverage(logs, path_to_produced)

        with open(path_to_produced, 'r') as produced_file:
            produced_data = produced_file.read()
        with open(path_to_expected, 'r') as expected_file:
            expected_data = expected_file.read()

        assert produced_data == expected_data, \
            "The contents of the stitched contigs plot" \
            " does not match the expected contents."

        return figure

    return check


def test_identical_stitching_of_one_contig(exact_aligner, visualizer):
    # Scenario: When stitching one contig, it remains the same.

    contigs = [
        GenotypedContig(name='a',
                        seq='ACTGACTG' * 100,
                        ref_name='testref',
                        group_ref='testref',
                        ref_seq='T' * 20 + 'ACTGACTG' * 110 + 'T' * 20,
                        match_fraction=1.0,
                        ),
        ]

    results = list(stitch_contigs(contigs))
    assert len(results) == 1
    assert results[0].seq == contigs[0].seq

    assert len(visualizer().elements) > len(contigs)


def test_separate_stitching_of_non_overlapping_contigs_1(exact_aligner, visualizer):
    # Scenario: When stitching multiple non-overlapping contigs, the order doesn't matter.

    ref_seq = 'A' * 100 + 'C' * 100

    contigs = [
        GenotypedContig(name='a',
                        seq='A' * 70,
                        ref_name='testref',
                        group_ref='testref',
                        ref_seq=ref_seq,
                        match_fraction=0.5,
                        ),
        GenotypedContig(name='b',
                        seq='C' * 70,
                        ref_name='testref',
                        group_ref='testref',
                        ref_seq=ref_seq,
                        match_fraction=0.5,
                        ),
        ]

    results = list(stitch_contigs(contigs))

    # No claims about the output order, so wrap into set()
    assert { contig.seq for contig in contigs } \
        == { contig.seq for contig in results }

    assert len(visualizer().elements) > len(contigs)


def test_separate_stitching_of_non_overlapping_contigs_2(exact_aligner, visualizer):
    # Scenario: When stitching multiple non-overlapping contigs, the order doesn't matter.

    ref_seq = 'A' * 100 + 'C' * 100

    contigs = [
        GenotypedContig(name='b',
                        seq='C' * 70,
                        ref_name='testref',
                        group_ref='testref',
                        ref_seq=ref_seq,
                        match_fraction=0.5,
                        ),
        GenotypedContig(name='a',
                        seq='A' * 70,
                        ref_name='testref',
                        group_ref='testref',
                        ref_seq=ref_seq,
                        match_fraction=0.5,
                        ),
        ]

    results = list(stitch_contigs(contigs))

    # No claims about the output order, so wrap into set()
    assert { contig.seq for contig in contigs } \
        == { contig.seq for contig in results }

    assert len(visualizer().elements) > len(contigs)


def test_correct_stitching_of_two_partially_overlapping_contigs(exact_aligner, visualizer):
    # Scenario: Two partially overlapping contigs are stitched correctly into a single sequence.

    ref_seq = 'A' * 100 + 'C' * 100

    contigs = [
        GenotypedContig(name='a',
                        seq='A' * 50 + 'C' * 20,
                        ref_name='testref',
                        group_ref='testref',
                        ref_seq=ref_seq,
                        match_fraction=0.5,
                        ),
        GenotypedContig(name='b',
                        seq='A' * 20 + 'C' * 50,
                        ref_name='testref',
                        group_ref='testref',
                        ref_seq=ref_seq,
                        match_fraction=0.5,
                        ),
        ]

    results = list(stitch_contigs(contigs))
    assert len(results) == 1

    result = results[0]

    assert 100 == len(result.seq)
    assert result.seq == 'A' * 50 + 'C' * 50

    assert len(visualizer().elements) > len(contigs)


def test_correct_processing_of_two_overlapping_and_one_separate_contig(exact_aligner, visualizer):
    # Scenario: Two overlapping contigs are stitched together, the non-overlapping is kept separate.

    ref_seq = 'A' * 100 + 'C' * 100 + 'T' * 100

    contigs = [
        GenotypedContig(name='a',
                        seq='A' * 50 + 'C' * 20,
                        ref_name='testref',
                        group_ref='testref',
                        ref_seq=ref_seq,
                        match_fraction=0.5,
                        ),
        GenotypedContig(name='b',
                        seq='A' * 20 + 'C' * 50,
                        ref_name='testref',
                        group_ref='testref',
                        ref_seq=ref_seq,
                        match_fraction=0.5,
                        ),
        GenotypedContig(name='c',
                        seq='C' * 20 + 'T' * 50,
                        ref_name='testref',
                        group_ref='testref',
                        ref_seq=ref_seq,
                        match_fraction=0.5,
                        ),
        ]

    results = list(stitch_contigs(contigs))
    assert len(results) == 2

    assert 100 == len(results[0].seq)
    assert results[0].seq == 'A' * 50 + 'C' * 50

    assert results[1].seq == contigs[2].seq

    assert len(visualizer().elements) > len(contigs)


def test_stitching_of_all_overlapping_contigs_into_one_sequence(exact_aligner, visualizer):
    # Scenario: All contigs have some overlapping parts, resulting in one continuous sequence after stitching.

    ref_seq = 'A' * 100 + 'C' * 100 + 'T' * 100

    contigs = [
        GenotypedContig(name='a',
                        seq='A' * 50 + 'C' * 20,
                        ref_name='testref',
                        group_ref='testref',
                        ref_seq=ref_seq,
                        match_fraction=0.5,
                        ),
        GenotypedContig(name='b',
                        seq='A' * 20 + 'C' * 100 + 'T' * 20,
                        ref_name='testref',
                        group_ref='testref',
                        ref_seq=ref_seq,
                        match_fraction=0.5,
                        ),
        GenotypedContig(name='c',
                        seq='C' * 20 + 'T' * 50,
                        ref_name='testref',
                        group_ref='testref',
                        ref_seq=ref_seq,
                        match_fraction=0.5,
                        ),
        ]

    results = list(stitch_contigs(contigs))
    assert len(results) == 1

    result = results[0]

    assert 200 == len(result.seq)
    assert result.seq == 'A' * 50 + 'C' * 100 + 'T' * 50

    assert len(visualizer().elements) > len(contigs)


def test_stitching_with_empty_contigs(exact_aligner, visualizer):
    # Scenario: The function is able to handle and ignore empty contigs.

    ref_seq = 'A' * 100

    contigs = [
        GenotypedContig(name='a',
                        seq=ref_seq,
                        ref_name='testref',
                        group_ref='testref',
                        ref_seq=ref_seq,
                        match_fraction=0.5,
                        ),
        GenotypedContig(name='b',
                        seq='',
                        ref_name='testref',
                        group_ref='testref',
                        ref_seq=ref_seq,
                        match_fraction=0.5,
                        ),
        ]

    results = list(stitch_contigs(contigs))
    assert { contig.seq for contig in contigs } \
        == { contig.seq for contig in results }

    assert len(visualizer().elements) > len(contigs)


def test_stitching_of_identical_contigs(exact_aligner, visualizer):
    # Scenario: The function correctly handles and avoids duplication when identical contigs are stitched together.

    contigs = [
        GenotypedContig(name=name,
                        seq='ACTGACTG' * 100,
                        ref_name='testref',
                        group_ref='testref',
                        ref_seq='ACTGACTG' * 100,
                        match_fraction=1.0,
                        )
        for name in ["a", "b", "c"]]

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


def test_correct_stitching_of_two_partially_overlapping_different_organism_contigs(exact_aligner, visualizer):
    # Scenario: Two partially overlapping contigs, but which come from different organism,
    # are not stitched into a single sequence.

    ref_seq = 'A' * 100 + 'C' * 100

    contigs = [
        GenotypedContig(name='a',
                        seq='A' * 50 + 'C' * 20,
                        ref_name='testref-1',
                        group_ref='testref-1',
                        ref_seq=ref_seq,
                        match_fraction=0.5,
                        ),
        GenotypedContig(name='b',
                        seq='A' * 20 + 'C' * 50,
                        ref_name='testref-2',
                        group_ref='testref-2',
                        ref_seq=ref_seq,
                        match_fraction=0.5,
                        ),
        ]

    results = list(stitch_contigs(contigs))
    assert len(results) == 2

    assert { contig.seq for contig in contigs } \
        == { contig.seq for contig in results }

    assert len(visualizer().elements) > len(contigs)


def test_correct_processing_complex_nogaps(exact_aligner, visualizer):
    # Scenario: There are two reference organisms.
    # Each with 4 contigs.
    # For each, three overlapping contigs are stitched together, the non-overlapping is kept separate.

    ref_seq = 'A' * 100 + 'C' * 100 + 'T' * 100 + 'G' * 100

    contigs = [[
        GenotypedContig(name='a',
                        seq='A' * 50 + 'C' * 20,
                        ref_name=ref_name,
                        group_ref=ref_name,
                        ref_seq=ref_seq,
                        match_fraction=0.5,
                        ),
        GenotypedContig(name='b',
                        seq='A' * 20 + 'C' * 50,
                        ref_name=ref_name,
                        group_ref=ref_name,
                        ref_seq=ref_seq,
                        match_fraction=0.5,
                        ),
        GenotypedContig(name='c',
                        seq='C' * 70 + 'T' * 20,
                        ref_name=ref_name,
                        group_ref=ref_name,
                        ref_seq=ref_seq,
                        match_fraction=0.5,
                        ),
        GenotypedContig(name='d',
                        seq='T' * 20 + 'G' * 50,
                        ref_name=ref_name,
                        group_ref=ref_name,
                        ref_seq=ref_seq,
                        match_fraction=0.5,
                        ),
        ] for ref_name in ['testref-1', 'testref-2']]

    contigs = sum(contigs, start=[])

    results = list(stitch_contigs(contigs))
    assert len(results) == 4

    assert 170 == len(results[0].seq)
    assert results[0].seq == 'A' * 50 + 'C' * 100 + 'T' * 20
    assert results[0].group_ref == 'testref-1'

    assert 170 == len(results[1].seq)
    assert results[1].seq == 'A' * 50 + 'C' * 100 + 'T' * 20
    assert results[1].group_ref == 'testref-2'

    assert results[2].seq == contigs[3].seq
    assert results[3].seq == contigs[7].seq

    assert len(visualizer().elements) > len(contigs)


def test_stitching_when_one_contig_completely_covered_by_another(exact_aligner, visualizer):
    # Scenario: If one contig is completely covered by another contig,
    # the completely covered contig must be dropped.

    ref_seq = 'A' * 100 + 'C' * 100

    contigs = [
        GenotypedContig(name='a',
                        seq='A' * 20 + 'C' * 20,
                        ref_name='testref',
                        group_ref='testref',
                        ref_seq=ref_seq,
                        match_fraction=0.5,
                        ),
        GenotypedContig(name='b',
                        seq='A' * 50 + 'C' * 50,
                        ref_name='testref',
                        group_ref='testref',
                        ref_seq=ref_seq,
                        match_fraction=0.5,
                        ),
        ]

    results = list(stitch_contigs(contigs))
    assert len(results) == 1

    # Test to ensure that the final result contains the contig 'b' and
    # does not contain the completely covered contig 'a'.
    assert results[0].seq == contigs[1].seq

    assert len(visualizer().elements) > len(contigs)


def test_stitching_contig_with_big_noncovered_gap(exact_aligner, visualizer):
    # Scenario: One contig has a big gap, which is however not covered by anything else.

    ref_seq = 'A' * 100 + 'C' * 100 + 'T' * 100

    contigs = [
        GenotypedContig(name='a',
                        seq= 'A' * 50 + 'T' * 50, # mind the C gap
                        ref_name='testref',
                        group_ref='testref',
                        ref_seq=ref_seq,
                        match_fraction=0.5,
                        ),
        ]

    results = list(stitch_contigs(contigs))

    assert { contig.seq for contig in contigs } \
        == { contig.seq for contig in results }

    assert len(visualizer().elements) > len(contigs)


def test_stitching_contig_with_big_noncovered_gap_2(exact_aligner, visualizer):
    # Scenario: One contig has a big gap, which is however not covered by anything else.

    ref_seq = 'A' * 100 + 'C' * 100 + 'T' * 100 + 'G' * 100

    contigs = [
        GenotypedContig(name='a',
                        seq='A' * 50 + 'T' * 50, # mind the C gap
                        ref_name='testref',
                        group_ref='testref',
                        ref_seq=ref_seq,
                        match_fraction=0.5,
                        ),
        GenotypedContig(name='B',
                        seq='G' * 50,
                        ref_name='testref',
                        group_ref='testref',
                        ref_seq=ref_seq,
                        match_fraction=0.5,
                        ),
        ]

    results = list(stitch_contigs(contigs))

    assert { contig.seq for contig in contigs } \
        == { contig.seq for contig in results }

    assert len(visualizer().elements) > len(contigs)


def test_stitching_contig_with_big_covered_gap(exact_aligner, visualizer):
    # Scenario: If one contig has a big gap covered by another contig.

    ref_seq = 'G' * 100 + 'A' * 100 + 'C' * 100 + 'T' * 100 + 'G' * 100

    contigs = [
        GenotypedContig(name='a',
                        seq='G' * 50 + 'A' * 50 + 'T' * 100, # mind the gap
                        ref_name='testref',
                        group_ref='testref',
                        ref_seq=ref_seq,
                        match_fraction=0.5,
                        ),
        GenotypedContig(name='b',
                        seq='A' * 100 + 'C' * 100 + 'T' * 50,
                        ref_name='testref',
                        group_ref='testref',
                        ref_seq=ref_seq,
                        match_fraction=0.5,
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

    ref_seq = 'G' * 100 + 'A' * 9 + 'C' * 100 + 'T' * 100

    contigs = [
        GenotypedContig(name='a',
                        seq='G' * 100 + 'A' * 0 + 'C' * 100, # mind the gap
                        ref_name='testref',
                        group_ref='testref',
                        ref_seq=ref_seq,
                        match_fraction=0.5,
                        ),
        GenotypedContig(name='b',
                        seq='A' * 9 + 'C' * 50,
                        ref_name='testref',
                        group_ref='testref',
                        ref_seq=ref_seq,
                        match_fraction=0.5,
                        ),
        ]

    contigs = list(align_all_to_reference(contigs))
    assert len(contigs) == 2
    assert len(list(contigs[0].alignment.deletions())) == 1
    assert len(list(contigs[1].alignment.deletions())) == 0
    results = list(split_contigs_with_gaps(contigs))
    assert len(results) == 3

    assert len(visualizer().elements) > len(contigs)

    assert all(x.seq == x.lstrip_query().rstrip_query().seq for x in results)
    assert { contig.seq for contig in contigs } \
        == { contig.seq for contig in results }


def test_stitching_partial_align(exact_aligner, visualizer):
    # Scenario: A single contig has a sequence that partially aligns to the reference sequence.

    contigs = [
        GenotypedContig(name='a',
                        seq='T' * 10 + 'C' * 20 + 'A' * 10,
                        ref_name='testref',
                        group_ref='testref',
                        ref_seq='A' * 20 + 'C' * 20 + 'T' * 20,
                        match_fraction=0.3,
                        ),
        ]

    results = list(stitch_contigs(contigs))
    assert len(results) == len(contigs)
    for result in results:
        assert any(result.seq in contig.seq for contig in contigs)

    assert len(visualizer().elements) > len(contigs)

    assert all(x.seq != x.lstrip_query().rstrip_query().seq for x in results)

    assert { contig.seq for contig in contigs } \
        != { contig.lstrip_query().rstrip_query().seq for contig in results }


def test_partial_align_consensus(exact_aligner, visualizer):
    # Scenario: A single contig partially aligns to the reference sequence, and a consensus sequence is being stitched.

    contigs = [
        GenotypedContig(name='a',
                        seq='T' * 10 + 'C' * 20 + 'A' * 10,
                        ref_name='testref',
                        group_ref='testref',
                        ref_seq='A' * 20 + 'C' * 20 + 'T' * 20,
                        match_fraction=0.3,
                        ),
        ]

    results = list(stitch_consensus(contigs))
    assert len(results) == len(contigs)
    assert { contig.seq for contig in contigs } \
        == { contig.seq for contig in results }

    assert len(visualizer().elements) > len(contigs)


def test_stitching_partial_align_multiple_sequences(exact_aligner, visualizer):
    # Scenario: Multiple contigs have sequences that partially align to the same reference sequence.

    ref_seq='A' * 20 + 'C' * 20 + 'T' * 20

    contigs = [
        GenotypedContig(name='a',
                        seq='T' * 10 + 'C' * 20 + 'A' * 10,
                        ref_name='testref',
                        group_ref='testref',
                        ref_seq=ref_seq,
                        match_fraction=0.3,
                        ),
        GenotypedContig(name='b',
                        seq='C' * 20 + 'A' * 10 + 'G' * 10,
                        ref_name='testref',
                        group_ref='testref',
                        ref_seq=ref_seq,
                        match_fraction=0.3,
                        ),
        ]

    results = list(stitch_contigs(contigs))
    assert len(results) == 1
    for result in results:
        assert any(result.seq in contig.seq for contig in contigs)

    assert len(visualizer().elements) > len(contigs)

    assert { contig.seq for contig in contigs } \
        != { contig.lstrip_query().rstrip_query().seq for contig in results }


def test_partial_align_consensus_multiple_sequences(exact_aligner, visualizer):
    # Scenario: Multiple contigs partially align to the same reference sequence, and a consensus sequence is being stitched from them.

    ref_seq='A' * 20 + 'C' * 20 + 'T' * 20

    contigs = [
        GenotypedContig(name='a',
                        seq='A' * 20,
                        ref_name='testref',
                        group_ref='testref',
                        ref_seq=ref_seq,
                        match_fraction=0.3,
                        ),
        GenotypedContig(name='b',
                        seq='T' * 20,
                        ref_name='testref',
                        group_ref='testref',
                        ref_seq=ref_seq,
                        match_fraction=0.3,
                        ),
        ]

    results = list(stitch_consensus(contigs))
    assert len(results) == 1
    assert results[0].seq == contigs[0].seq + contigs[1].seq

    assert len(visualizer().elements) > len(contigs)


def test_partial_align_consensus_multiple_overlaping_sequences(exact_aligner, visualizer):
    # Scenario: Multiple contigs partially align to the same reference sequence, and a consensus sequence is being stitched from them.

    ref_seq='A' * 20 + 'C' * 20 + 'T' * 20

    contigs = [
        GenotypedContig(name='a',
                        seq='T' * 10 + 'A' * 5 + 'C' * 20 + 'A' * 10,
                        ref_name='testref',
                        group_ref='testref',
                        ref_seq=ref_seq,
                        match_fraction=0.3,
                        ),
        GenotypedContig(name='b',
                        seq='C' * 20 + 'T' * 5 + 'A' * 10 + 'G' * 10,
                        ref_name='testref',
                        group_ref='testref',
                        ref_seq=ref_seq,
                        match_fraction=0.3,
                        ),
        ]

    results = list(stitch_consensus(contigs))
    assert len(results) == 1
    assert results[0].seq == 'T' * 10 + 'A' * 5 + 'C' * 20 + 'T' * 5 + 'A' * 10 + 'G' * 10
    assert results[0].seq == contigs[0].seq[:-10] + contigs[1].seq[20:]

    assert len(visualizer().elements) > len(contigs)


def test_main_invocation(exact_aligner, tmp_path, hcv_db):
    pwd = os.path.dirname(__file__)
    contigs = os.path.join(pwd, "data", "exact_parts_contigs.csv")
    stitched_contigs = os.path.join(tmp_path, "stitched.csv")
    stitcher.main([contigs, stitched_contigs])

    assert os.path.exists(contigs)
    assert os.path.exists(stitched_contigs)

    # Check the contents of stitched_contigs
    with open(stitched_contigs, 'r') as stitched_file:
        stitched_data = stitched_file.read()

    expected_file_path = os.path.join(pwd, "data", "exact_parts_contigs_stitched.csv")
    with open(expected_file_path, 'r') as expected_file:
        expected_data = expected_file.read()

    assert stitched_data == expected_data, "The contents of the stitched contigs file do not match the expected contents."


def test_visualizer_simple(exact_aligner, tmp_path, hcv_db):
    pwd = os.path.dirname(__file__)
    contigs = os.path.join(pwd, "data", "exact_parts_contigs.csv")
    stitched_contigs = os.path.join(tmp_path, "stitched.csv")
    plot = os.path.join(tmp_path, "exact_parts_contigs.plot.svg")
    stitcher.main([contigs, stitched_contigs, "--debug", "--plot", plot])

    assert os.path.exists(contigs)
    assert os.path.exists(stitched_contigs)

    # Check the contents of stitched_contigs
    with open(stitched_contigs, 'r') as stitched_file:
        stitched_data = stitched_file.read()

    expected_file_path = os.path.join(pwd, "data", "exact_parts_contigs_stitched.csv")
    with open(expected_file_path, 'r') as expected_file:
        expected_data = expected_file.read()
        assert stitched_data == expected_data, "The contents of the stitched contigs file do not match the expected contents."

    # Check the contents of stitched_contigs
    expected_plot = os.path.join(pwd, "data", "exact_parts_contigs.plot.svg")
    with open(plot, 'r') as stitched_file, \
         open(expected_plot, 'r') as expected_file:
        stitched_data = stitched_file.read()
        expected_data = expected_file.read()
        assert stitched_data == expected_data, "The contents of the stitched plot file do not match the expected contents."


#  _   _       _ _     _            _
# | | | |_ __ (_) |_  | |_ ___  ___| |_ ___
# | | | | '_ \| | __| | __/ _ \/ __| __/ __|
# | |_| | | | | | |_  | ||  __/\__ \ |_\__ \
#  \___/|_| |_|_|\__|  \__\___||___/\__|___/
#

@pytest.mark.parametrize("intervals, expected", [
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
    ([(1, 3), (2, 4), (6, 8), (10, 11), (11, 12)],
     [(1, 4), (6, 8), (10, 12)]),

    # Intervals out of initial order
    ([(4, 6), (1, 2)],
     [(1, 2), (4, 6)]),

    # Overlapping intervals with out of order inputs
    ([(1, 4), (3, 5), (2, 3), (7, 10), (9, 12)],
     [(1, 5), (7, 12)]),

    # Large set of intervals with various overlaps
    ([(1, 4), (2, 6), (5, 8), (7, 8), (10, 15), (11, 12), (13, 14), (17, 18)],
     [(1, 8), (10, 15), (17, 18)]),

    # Intervals where end is less than start should return as is or be handled explicitly depending on implementation
    ([(5, 3), (1, 2)],
     [(1, 2), (5, 3)]),

    # Intervals that are exactly one after the other in sequence / Intervals that are completely disjoint
    ([(1, 2), (4, 5), (7, 8)],
     [(1, 2), (4, 5), (7, 8)]),

    # Overlapping intervals that merge into one large interval
    ([(2, 6), (4, 10), (5, 15), (14, 20)],
     [(2, 20)]),

    # Same interval repeated multiple times
    ([(1, 5), (1, 5), (1, 5)],
     [(1, 5)]),

    # Single point intervals
    ([(1, 1), (5, 5), (3, 3)],
     [(1, 1), (3, 3), (5, 5)]),

    ([(1, 1), (5, 5), (3, 3), (1, 1), (1, 1)],
     [(1, 1), (3, 3), (5, 5)]),

    ([(1, 1), (2, 3)],
     [(1, 3)]),

    # Intervals that start with negative numbers
    ([(-5, 0), (-2, 3), (1, 7), (9, 12)],
     [(-5, 7), (9, 12)]),
])
def test_merge_intervals(intervals, expected):
    assert merge_intervals(intervals) == expected


class MockAlignedContig:
    def __init__(self, ref_name, group_ref, r_st, r_ei, name="contig"):
        self.ref_name = ref_name
        self.group_ref = group_ref
        self.alignment = MockAlignment(r_st, r_ei)
        self.name = name

    def overlaps(self, other):
        return AlignedContig.overlaps(self, other)


class MockAlignment:
    def __init__(self, r_st, r_ei):
        self.r_st = r_st
        self.r_ei = r_ei


# Simple function to create mock AlignedContig objects for testing, including ref_name.
def create_mock_aligned_contig(ref_name, r_st, r_ei, name="contig"):
    return MockAlignedContig(ref_name, ref_name, r_st, r_ei, name)


@pytest.mark.parametrize("contigs, expected_covered_name", [
    # No contigs are completely covered.
    ([('ref1', 0, 100), ('ref1', 101, 200)], None),
    ([('ref1', 0, 50), ('ref1', 51, 100)], None),

    # A single contig is completely covered by one other contig.
    ([('ref1', 0, 100), ('ref1', 0, 200)], 'contig1'),
    ([('ref1', 50, 150), ('ref1', 0, 200)], 'contig1'),

    # A single contig completely covers another, but with different reference names.
    ([('ref1', 0, 50), ('ref2', 0, 100)], None),

    # Single coverage with exact match.
    ([('ref1', 0, 100), ('ref1', 0, 100)], 'contig1'),

    # A single contig is completely covered at the beginning by one and at the end by another contig.
    ([('ref1', 0, 50), ('ref1', 50, 100), ('ref1', 25, 75)], 'contig3'),

    # Contigs overlap but none are completely covered.
    ([('ref1', 0, 50), ('ref1', 40, 90), ('ref1', 80, 120)], None),

    # Multiple contigs with some covered completely by a single other contig.
    ([('ref1', 0, 200), ('ref1', 10, 30), ('ref1', 170, 190)], 'contig2'),

    # Multiple contigs with complex overlaps and one completely covered.
    ([('ref1', 30, 60), ('ref1', 0, 50), ('ref1', 20, 70), ('ref1', 60, 90)], 'contig1'),

    # Edge case where a contig starts where another ends.
    ([('ref1', 0, 50), ('ref1', 50, 100)], None),

    # Contigs are completely covered in a nested fashion.
    ([('ref1', 0, 200), ('ref1', 50, 150), ('ref1', 100, 125)], 'contig2'),

    # Contigs are adjacent and cover each other completely.
    ([('ref1', 0, 100), ('ref1', 101, 200), ('ref1', 0, 200)], 'contig1'),

    # Single large contig covers several smaller non-adjacent contigs.
    ([('ref1', 0, 500), ('ref1', 50, 100), ('ref1', 200, 250), ('ref1', 300, 350)], 'contig2'),

    # Single large contig covers several smaller adjacent contigs.
    ([('ref1', 50, 100), ('ref1', 70, 300), ('ref1', 101, 199), ('ref1', 200, 350)], 'contig2'),

    # Single small contig is covered by several larger contigs.
    ([('ref1', 0, 250), ('ref1', 200, 300), ('ref1', 600, 800), ('ref1', 250, 700)], 'contig2'),

    # Complex case with multiple contigs and complete coverage by combinations.
    ([('ref1', 0, 100), ('ref1', 30, 130), ('ref1', 60, 160), ('ref1', 90, 190), ('ref1', 120, 220)], 'contig2'),

    # Contigs with same start but different end, where one is covered.
    ([('ref1', 0, 100), ('ref1', 0, 50)], 'contig2'),

    # Contigs with same end but different start, where one is covered.
    ([('ref1', 50, 100), ('ref1', 0, 100)], 'contig1'),

    # Contig covered by two overlapping contigs that don't individually cover the whole range.
    ([('ref1', 0, 75), ('ref1', 25, 100), ('ref1', 0, 100)], 'contig1'),

    # Two contigs are covered completely by one large contig.
    ([('ref1', 0, 300), ('ref1', 50, 100), ('ref1', 200, 250)], 'contig2'),

    # No contigs at all.
    ([], None),
])
def test_find_covered(contigs, expected_covered_name):
    mock_contigs = [create_mock_aligned_contig(ref_name, r_st, r_ei, f'contig{i+1}')
                    for i, (ref_name, r_st, r_ei) in enumerate(contigs)]
    covered, covering = find_covered_contig(mock_contigs)
    if expected_covered_name is None:
        assert covered is None
    else:
        assert covered is not None
        assert covered.name == expected_covered_name


def test_concordance_same_length_inputs():
    with pytest.raises(ValueError):
        calculate_concordance('abc', 'ab')

def test_concordance_completely_different_strings():
    result = calculate_concordance('a'*30, 'b'*30)
    assert all(n == 0 for n in result)

def generate_random_string_pair(length):
    left = ''.join(random.choice('ACGT') for _ in range(length))
    right = ''.join(random.choice('ACGT') for _ in range(length))
    return left, right


@pytest.mark.parametrize(
    'left, right, expected',
    [("aaaaa", "aaaaa", [0.6, 0.68, 0.7, 0.68, 0.6]),
     ("abcdd", "abcdd", [0.6, 0.68, 0.7, 0.68, 0.6]),
     ("aaaaaaaa", "baaaaaab", [0.3, 0.62, 0.71, 0.75, 0.75, 0.71, 0.62, 0.3]),
     ("aaaaaaaa", "aaaaaaab", [0.64, 0.73, 0.79, 0.8, 0.79, 0.73, 0.64, 0.31]),
     ("aaaaaaaa", "aaaaaaab", [0.64, 0.73, 0.79, 0.8, 0.79, 0.73, 0.64, 0.31]),
     ("aaaaaaaa", "aaaaabbb", [0.6, 0.68, 0.7, 0.68, 0.6, 0.29, 0.19, 0.13]),
     ("aaaaaaaa", "aaabbaaa", [0.56, 0.63, 0.62, 0.39, 0.39, 0.62, 0.63, 0.56]),
     ("aaaaa", "bbbbb", [0] * 5),
     ]
)
def test_concordance_simple(left, right, expected):
    result = [round(float(x), 2) for x in calculate_concordance(left, right)]
    assert result == expected


@pytest.mark.parametrize(
    'left, right, expected',
    [("a" * 128, "a" * 128, 64),
     ("a" * 128, "a" * 64 + "b" * 64, 32),
     ("a" * 128, "a" * 64 + "ba" * 32, 32),
     ("a" * 128, "a" * 54 + "b" * 20 + "a" * 54, 28), # two peaks
     ("a" * 128, "a" * 63 + "b" * 2 + "a" * 63, 32), # two peaks
     ("a" * 1280, "b" * 640 + "a" * 640, round(1280 * 3 / 4)),
     ]
)
def test_concordance_simple_index(left, right, expected):
    concordance = calculate_concordance(left, right)
    concordance_d = list(disambiguate_concordance(concordance))
    index = max(range(len(concordance)), key=lambda i: concordance_d[i])
    if abs(index - expected) > 3:
        assert index == expected


def generate_test_cases(num_cases):
    with fixed_random_seed(42):
        length = random.randint(1, 80)
        return [generate_random_string_pair(length) for _ in range(num_cases)]

concordance_cases = generate_test_cases(num_cases=100)


@pytest.mark.parametrize('left, right', concordance_cases)
def test_concordance_output_range(left, right):
    result = calculate_concordance(left, right)
    assert all(0 <= n <= 1 for n in result), "All values in result should be between 0 and 1"


@pytest.mark.parametrize('left, right', concordance_cases)
def test_concordance_higher_if_more_matches_added(left, right):
    # Insert exact matches in the middle
    matching_sequence = 'A' * 30
    insert_position = len(left) // 2
    new_left = left[:insert_position] + matching_sequence + left[insert_position + len(matching_sequence):]
    new_right = right[:insert_position] + matching_sequence + right[insert_position + len(matching_sequence):]

    old_conc = calculate_concordance(left, right)
    new_conc = calculate_concordance(new_left, new_right)
    old_average = sum(old_conc) / len(old_conc)
    new_average = sum(new_conc) / len(new_conc)
    assert old_average <= new_average


@pytest.mark.parametrize('left, right', concordance_cases)
def test_concordance_higher_in_matching_areas(left, right):
    # Insert exact matches in the middle
    matching_sequence = 'A' * 30
    insert_position = len(left) // 2
    new_left = left[:insert_position] + matching_sequence + left[insert_position + len(matching_sequence):]
    new_right = right[:insert_position] + matching_sequence + right[insert_position + len(matching_sequence):]

    concordance_scores = calculate_concordance(new_left, new_right)

    # Check concordance in the matching area
    matching_area_concordance = concordance_scores[insert_position:insert_position + len(matching_sequence)]

    # Calculate average concordance inside and outside the matching area
    average_inside = sum(matching_area_concordance) / len(matching_sequence)
    average_outside = (sum(concordance_scores) - sum(matching_area_concordance)) / (len(concordance_scores) - len(matching_sequence))

    # Assert that the concordance is indeed higher in the matching area
    assert average_inside > average_outside, "Concordance in matching areas should be higher than in non-matching areas"
