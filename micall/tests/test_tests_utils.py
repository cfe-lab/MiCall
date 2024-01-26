
import pytest
from micall.tests.utils import MockAligner, MockAlignment

def test_basic_mapping():

    aligner = MockAligner('acgt' + 'a' * 20 + 'acgt')

    alignment = list(aligner.map('a' * 10))

    assert len(alignment) == 5

    alignment = alignment[0]

    assert isinstance(alignment, MockAlignment)
    assert alignment.mapq == 60
    assert alignment.strand == 1
    assert alignment.r_st == 4
    assert alignment.r_en == 14
    assert alignment.q_st == 0
    assert alignment.q_en == 10


def test_exact_match():
    aligner = MockAligner("abcdefg")
    alignments = list(aligner.map("abc"))
    assert len(alignments) == 1
    assert alignments[0].r_st == 0
    assert alignments[0].r_en == 3


def test_no_match():
    aligner = MockAligner("abcdefg")
    alignments = list(aligner.map("xyz"))
    assert len(alignments) == 0


def test_partial_match():
    aligner = MockAligner("abcdefg")
    alignments = list(aligner.map("abxyabc"))
    assert len(alignments) == 1
    assert alignments[0].r_st == 0
    assert alignments[0].r_en == 3


def test_multiple_matches():
    aligner = MockAligner("A" * 40)
    alignments = list(aligner.map("A" * 20))
    assert len(alignments) == 5
    assert alignments[0].r_st == 0
    assert alignments[0].r_en == 20
    assert alignments[1].r_st == 20
    assert alignments[1].r_en == 40


def test_multiple_matches_bigger_query():
    aligner = MockAligner("A" * 40)
    alignments = list(aligner.map("A" * 50))
    assert len(alignments) == 5
    assert alignments[0].r_st == 0
    assert alignments[0].r_en == 40
    assert alignments[1].r_st == 0
    assert alignments[1].r_en == 40


def test_empty_reference():
    aligner = MockAligner("A" * 0)
    alignments = list(aligner.map("A" * 20))
    assert len(alignments) == 0


def test_empty_query():
    aligner = MockAligner("A" * 40)
    alignments = list(aligner.map("A" * 0))
    assert len(alignments) == 0
