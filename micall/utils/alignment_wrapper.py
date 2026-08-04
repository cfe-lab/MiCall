import typing

from gotoh import align_it


def align_nucs(seq1: str, seq2: str) -> typing.Tuple[str, str, int]:
    """ Align two sequences of nucleotides with default parameters.

    :return: the two aligned sequences, plus an alignment score
    """
    gap_open_penalty = 15
    gap_extend_penalty = 3
    use_terminal_gap_penalty = 1
    aligned1, aligned2, score = align_it(
        seq1,
        seq2,
        gap_open_penalty,
        gap_extend_penalty,
        use_terminal_gap_penalty)

    return aligned1, aligned2, score
