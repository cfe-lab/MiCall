#! /usr/bin/env python

import argparse
from dataclasses import dataclass
import sys
import numpy as np
from typing import Sequence, Optional, Tuple, Iterable, Any
from itertools import chain
import scipy


@dataclass(frozen=False)
class OverlapFinder:
    alphabet: Tuple[object, ...]
    bit_arr1: np.ndarray[Any, np.dtype[Any]]
    bit_arr2: np.ndarray[Any, np.dtype[Any]]
    total: np.ndarray[Any, np.dtype[Any]]

    @staticmethod
    def make(alphabet: Iterable[object]) -> 'OverlapFinder':
        # FIXME: Automatically resize these arrays.
        total = np.zeros(100_000)
        bit_arr1 = np.zeros(100_000)
        bit_arr2 = np.zeros(100_000)
        return OverlapFinder(tuple(alphabet),
                             total=total,
                             bit_arr1=bit_arr1,
                             bit_arr2=bit_arr2,
                             )


def find_maximum_overlap(arr1: Sequence[object],
                         arr2: Sequence[object],
                         finder: Optional[OverlapFinder] = None,
                         ) -> int:
    """
    Calculate the offset at which two sequences (arr1 and arr2)
    overlap the most.

    This function uses a convolution-based approach to determine the
    degree of overlap between two sequences of objects. It considers
    each unique element in both sequences, calculates how well they
    align at different offsets, and returns the offset index that
    results in the maximum overlap.

    Parameters:
    arr1 (Sequence[object]): The first sequence of objects.
    arr2 (Sequence[object]): The second sequence of objects.

    Returns:
    int: The shift value for `arr1` relative to `arr2` to achieve the
    maximum overlap. Value of 0 aligns `arr1` before `arr2` with no
    overlap. Negative value shifts the start of `arr2` back, creating
    more and more overlap, until it shifts it past the start of `arr1`,
    when the overlap starts to decrease again.
    """

    if len(arr1) == 0 or len(arr2) == 0:
        raise ValueError(
            f"Expected non-empty sequences, but got {len(arr1)}, {len(arr2)}.")

    if finder is None:
        # Create a set to store unique elements from both sequences
        alphabet_set = set()
        alphabet = []

        # Iterate over each element in both sequences to build a list of
        # unique elements
        for x in chain(arr1, arr2):
            if x not in alphabet_set:
                alphabet_set.add(x)
                alphabet.append(x)

        finder = OverlapFinder.make(alphabet)

    # Initialize an array to accumulate convolved results for
    # determining overlap
    len_total = len(arr1) + len(arr2) - 1

    # Slicing a NumPy array does not create a copy of the original array.
    # It creates a view.
    total = finder.total[:len_total]
    total.fill(0)

    bit_arr1 = finder.bit_arr1[:len(arr1)]
    bit_arr2 = finder.bit_arr2[:len(arr2)]

    np_arr1 = np.array(tuple(arr1))
    np_arr2 = np.array(tuple(reversed(arr2)))

    # Iterate over each unique element to determine overlap
    for element in finder.alphabet:
        bit_arr1.fill(0)
        bit_arr2.fill(0)
        bit_arr1[np_arr1 == element] = 1
        bit_arr2[np_arr2 == element] = 1

        # Compute the convolution of the two binary arrays
        convo = scipy.signal.convolve(bit_arr1, bit_arr2, mode='full')
        # Add the convolution to the total results
        total += convo

    # Return the shift value that yields maximum overlap
    max_value = np.max(total)
    if max_value <= 0:
        return 0

    max_indices = np.where(total == max_value)[0]
    left_max = max_indices[0]
    right_max = max_indices[-1]
    if left_max < (len_total - right_max):
        max_offset = left_max
    else:
        max_offset = right_max

    shift_value = max_offset - len_total
    return int(shift_value)


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
    shift = find_maximum_overlap(seq1, seq2)
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
