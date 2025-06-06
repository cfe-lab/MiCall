#! /usr/bin/env python

import argparse
from dataclasses import dataclass
import sys
import numpy as np
from typing import Sequence, Optional, Tuple, Iterable, Any
from itertools import chain
import scipy
from micall.utils.overlap_stitcher import calculate_overlap_score


@dataclass(frozen=False)
class OverlapFinder:
    alphabet: Tuple[object, ...]
    bit_arr1: np.ndarray[Any, np.dtype[Any]]
    bit_arr2: np.ndarray[Any, np.dtype[Any]]
    total: np.ndarray[Any, np.dtype[Any]]

    @staticmethod
    def make(alphabet: Iterable[str]) -> 'OverlapFinder':
        # FIXME: Automatically resize these arrays.
        total = np.zeros(100_000)
        bit_arr1 = np.zeros(100_000)
        bit_arr2 = np.zeros(100_000)
        alphabet_keys = {x: True for x in alphabet}.keys()
        alphabet_init = tuple(x.encode('utf-8') for x in alphabet_keys)
        return OverlapFinder(alphabet=alphabet_init,
                             total=total,
                             bit_arr1=bit_arr1,
                             bit_arr2=bit_arr2,
                             )


def choose_convolution_method(len1: int, len2: int) -> Any:
    if len1 * len2 > 10_000 * 10_000:
        method = scipy.signal.convolve
    else:
        method = np.convolve
    return method


def get_overlap_results(total: np.ndarray,
                        len_1: int, len_2: int,
                        ) -> Tuple[int, float]:
    len_total = len(total)

    max_overlap = min(len_1, len_2)
    def clip(arr):
        return np.clip(arr, 0, max_overlap + 1)

    overlap_sizes = np.zeros(len(total))
    left_size = len(total) // 2
    right_size = len(total) - left_size
    overlap_sizes[:left_size] = clip(1 + np.arange(left_size))
    overlap_sizes[-right_size:] = np.flip(clip(1 + np.arange(right_size)))

    assert len(overlap_sizes[:left_size]) + len(overlap_sizes[-right_size:]) == len(total)
    assert np.max(overlap_sizes[:left_size]) <= max_overlap + 1
    assert np.max(overlap_sizes[-right_size:]) <= max_overlap + 1
    assert np.min(overlap_sizes[:left_size]) == np.min(overlap_sizes[-right_size:]) == 1

    total = calculate_overlap_score(L=overlap_sizes, M=total)  # type: ignore

    # Return the shift value that yields maximum overlap
    max_value = np.max(total)
    if max_value <= 0:
        return (0, 0)

    max_indices = np.where(total == max_value)[0]
    left_max = max_indices[0]
    right_max = max_indices[-1]
    if left_max < (len_total - right_max):
        max_offset = left_max
    else:
        max_offset = right_max

    shift = max_offset - len_total
    return (int(shift), float(max_value))


def find_maximum_overlap(seq1: str,
                         seq2: str,
                         finder: Optional[OverlapFinder] = None,
                         ) -> Tuple[int, float]:
    """
    Calculate the offset at which two sequences (seq1 and seq2)
    overlap the most.

    This function uses a convolution-based approach to determine the
    degree of overlap between two sequences of objects. It considers
    each unique element in both sequences, calculates how well they
    align at different offsets, and returns the offset index that
    results in the maximum overlap.

    Parameters:
    seq1 (str): The first sequence of objects.
    seq2 (str): The second sequence of objects.

    Returns:
    int: The shift value for `seq1` relative to `seq2` to achieve the
    maximum overlap. Value of 0 aligns `seq1` before `seq2` with no
    overlap. Negative value shifts the start of `seq2` back, creating
    more and more overlap, until it shifts it past the start of `seq1`,
    when the overlap starts to decrease again.
    """

    if len(seq1) == 0 or len(seq2) == 0:
        raise ValueError(
            f"Expected non-empty sequences, but got {len(seq1)}, {len(seq2)}.")

    if finder is None:
        finder = OverlapFinder.make(chain(seq1, seq2))

    # Initialize an array to accumulate convolved results for
    # determining overlap
    len_total = len(seq1) + len(seq2) - 1

    # Slicing a NumPy array does not create a copy of the original array.
    # It creates a view.
    total = finder.total[:len_total]
    total.fill(0)

    bit_arr1 = finder.bit_arr1[:len(seq1)]
    bit_arr2 = finder.bit_arr2[:len(seq2)]

    np_arr1 = np.frombuffer(seq1.encode('utf-8'), dtype='S1')
    np_arr2 = np.flip(np.frombuffer(seq2.encode('utf-8'), dtype='S1'))

    method = choose_convolution_method(len(bit_arr1), len(bit_arr2))

    # Iterate over each unique element to determine overlap
    for element in finder.alphabet:
        bit_arr1.fill(0)
        bit_arr2.fill(0)
        bit_arr1[np_arr1 == element] = 1
        bit_arr2[np_arr2 == element] = 1

        # Compute the convolution of the two binary arrays
        convo = method(bit_arr1, bit_arr2, mode='full')

        # Add the convolution to the total results
        total += convo

    return get_overlap_results(total, len(bit_arr1), len(bit_arr2))


def show_maximum_overlap(arr1: Sequence[object],
                         arr2: Sequence[object],
                         shift: int,
                         ) -> str:

    arr1_line = '-' * (abs(shift) - len(arr1)) + ''.join(map(str, arr1))
    arr2_line = '-' * (len(arr1) - abs(shift)) + ''.join(map(str, arr2))
    max_len = max(len(arr1_line), len(arr2_line))
    lines = (x.ljust(max_len, '-')
             for x in [arr1_line, arr2_line])

    return '\n'.join(chain(lines, ['']))


def main(seq1: str, seq2: str) -> int:
    (shift, value) = find_maximum_overlap(seq1, seq2)
    ret = show_maximum_overlap(seq1, seq2, shift)
    print(ret, end='')
    return 0


def get_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser("Find maximum overlap between two sequences.")
    parser.add_argument('seq1')
    parser.add_argument('seq2')
    return parser


def cli_main(argv: Sequence[str]) -> int:
    parser = get_parser()
    args = parser.parse_args(argv)
    return main(args.seq1, args.seq2)


if __name__ == '__main__':
    sys.exit(cli_main(sys.argv[1:]))
