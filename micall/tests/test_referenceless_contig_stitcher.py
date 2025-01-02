import pytest
from micall.utils.referenceless_contig_stitcher import stitch_consensus
from micall.utils.contig_stitcher_contigs import Contig


@pytest.fixture(autouse=True)
def disable_acceptable_prob_check(monkeypatch):
    monkeypatch.setattr("micall.utils.referenceless_contig_stitcher.ACCEPTABLE_STITCHING_PROB", 1)


@pytest.mark.parametrize(
    "seqs, expected",
    [
        #
        # Singletons. Copied from `find_maximum_overlap`.
        #

        (('aaaaaxxxx', 'xxxbbbbb'), ('aaaaaxxxxbbbbb',)),
        (('aaaaaxxxx', 'xxxbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb'), ('aaaaaxxxxbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb',)),
        (('aaaaaxxxxcccccccccccccccccccc', 'xxxbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb'), ('aaaaaxxxxcccccccccccccccccccc', 'xxxbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb')),
        (('aaaaaxxxxx', 'bbbbbxxx'), ("bbbbbxxxxx",)),
        (('bbbbbxxx', 'aaaaaxxxxx'), ('aaaaaxxxxx',)),
        (('aaaaaxxx', 'bbbbbxxxxx'), ('bbbbbxxxxx',)),
        (('xxxaaaaa', 'xxxxxbbbbb'), ('xxxaaaaabb',)),
        (('xxxxxaaaaa', 'xxxbbbbb'), ('xxxbbbbbaa',)),
        (('xxxxaaaaa', 'bbbbbxxx'), ('bbbbbxxxxaaaaa',)),
        (('aaaaxxxxaaaaa', 'bxxx'), ('aaaaxxxxaaaaa',)),
        (('aaaa', 'bbbb'), ('aaaa', 'bbbb',)),
        (('aaaax', 'xbbbb'), ('aaaaxbbbb',)),

        #
        # Multiple.
        #
        (('aaax', 'xbbby', 'ycccz', 'zddd'), ('aaaxbbbyccczddd',)),
        (('aaa',), ('aaa',)),
        ((), ()),
    ],
)
def test_stitch_simple_cases(seqs, expected):
    contigs = [Contig(None, seq) for seq in seqs]
    consenses = tuple(contig.seq for contig in stitch_consensus(contigs))
    assert consenses == expected
