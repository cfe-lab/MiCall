import pytest
from micall.utils.referenceless_contig_stitcher import stitch_consensus
from micall.utils.contig_stitcher_contigs import Contig


@pytest.mark.parametrize(
    "seqs, expected",
    [
        #
        # Singletons. Copied from `find_maximum_overlap`.
        #

        (('aaaaaxxxx', 'xxxbbbbb'), ('aaaaaxxxxbbbbb',)),
        (('aaaaaxxxx', 'xxxbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb'), ('aaaaaxxxxbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb',)),
        (('aaaaaxxxxcccccccccccccccccccc', 'xxxbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb'), ('aaaaaxxxxccccccccccccccccccccbbbbbbbbbbb',)),
        (('aaaaaxxxxx', 'bbbbbxxx'), ("bbbbbxxx",)),
        (('bbbbbxxx', 'aaaaaxxxxx'), ('aaaaaxxxxx',)),
        (('aaaaaxxx', 'bbbbbxxxxx'), ('bbbbbxxxxx',)),
        (('xxxaaaaa', 'xxxxxbbbbb'), ('xxxxxbbbbb',)),
        (('xxxxxaaaaa', 'xxxbbbbb'), ('xxxbbbbb',)),
        (('xxxxaaaaa', 'bbbbbxxx'), ('bbbbbxxx',)),
        (('aaaaxxxxaaaaa', 'bxxx'), ('aaaaxxxxaaaaa',)),
        (('aaaa', 'bbbb'), ('aaaa', 'bbbb',)),
        (('aaaax', 'xbbbb'), ('aaaaxbbbb',)),

        #
        # Multiple.
        #
        (('aaax', 'xbbby', 'ycccz', 'zddd'), ('aaaxbbbyccczddd',)),
    ],
)
def test_stitch_simple_cases(seqs, expected):
    contigs = [Contig(None, seq) for seq in seqs]
    consenses = tuple(contig.seq for contig in stitch_consensus(contigs))
    assert consenses == expected
