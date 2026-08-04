from micall.utils.probe_finder import unpack_mixtures_and_reverse, ProbeFinder
from micall.utils.translation import reverse_and_complement


def test_unpack_mixtures_and_reverse():
    seq = 'ACTGRCT'  # R = A+G

    # {(seq, is_reversed)}
    expected_seqs = {('ACTGACT', False),
                     ('ACTGGCT', False),
                     # These two are reverse complements of the previous two.
                     ('AGTCAGT', True),
                     ('AGCCAGT', True)}

    seqs = unpack_mixtures_and_reverse(seq)

    assert seqs == expected_seqs


def test_probe_finder():
    target_seq = 'ATCGACCTAGCT'
    contig_seq = 'ATCGACCTAGCTAATTCCAGT'

    finder = ProbeFinder(contig_seq, target_seq)

    assert finder.contig_match == target_seq
    assert not finder.is_reversed
    assert finder.dist == 0
    assert finder.score == 60


def test_probe_finder_with_changes():
    target_seq = '    ATCGACCTAGCT'.strip()  # align for readability
    contig_seq = '    ATCGACCTGGCTAATTCCAGT'.strip()
    expected_match = 'ATCGACCTGGCT'
    # change                  ^

    finder = ProbeFinder(contig_seq, target_seq)

    assert finder.contig_match == expected_match
    assert finder.dist == 1
    assert finder.end_dist == 1


def test_probe_finder_reversed():
    target_seq = 'ATCGACCTAGCT'
    contig_seq = reverse_and_complement('ATCGACCTGGCTAATTCCAGT')
    expected_match = 'ATCGACCTGGCT'

    finder = ProbeFinder(contig_seq, target_seq)

    assert finder.contig_match == expected_match
    assert finder.is_reversed


def test_probe_finder_before_start():
    target_seq = ' ATCGACCTAGCT'.strip()
    contig_seq = '    GACCTAGCTAATTCCAGT'.strip()
    expected_match = 'GACCTAGCT'

    finder = ProbeFinder(contig_seq, target_seq)

    assert finder.contig_match == expected_match
    assert not finder.is_reversed
    assert finder.dist == 3
    assert finder.end_dist == 0


def test_probe_finder_after_end():
    target_seq = '         ATCGACCTAGCT'.strip()
    contig_seq = 'AATTCCAGTATCGTCCTA'
    expected_match = '     ATCGTCCTA'.strip()

    finder = ProbeFinder(contig_seq, target_seq)

    assert finder.contig_match == expected_match
    assert not finder.is_reversed
    assert finder.dist == 4
    assert finder.end_dist == 1
