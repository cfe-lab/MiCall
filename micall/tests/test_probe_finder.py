from micall.utils.probe_finder import unpack_mixtures


def test_unpack_mixtures():
    seq = 'ACTGRCT'  # R = A+G
    expected_seqs = {'ACTGACT',
                     'ACTGGCT'}

    seqs = unpack_mixtures(seq)

    assert seqs == expected_seqs
