from typing import Iterable, Iterator, Optional, FrozenSet, Tuple, Sequence, TextIO
from dataclasses import dataclass
from fractions import Fraction
import csv

from micall.utils.contig_stitcher_contigs import Contig
from micall.utils.find_maximum_overlap import find_maximum_overlap


@dataclass(frozen=True)
class ContigsPath:
    # Contig representing all combined contigs in the path.
    whole: Contig

    # Id's of contigs that comprise this path.
    parts_ids: FrozenSet[int]

    # Lower is better. This is an estimated probability that
    # all the components in this path came together by accident.
    score: Fraction

    def has_contig(self, contig: Contig) -> bool:
        return contig.id in self.parts_ids

    @property
    def is_empty(self) -> bool:
        return len(self.parts_ids) == 0

    @staticmethod
    def empty() -> 'ContigsPath':
        return ContigsPath(Contig(None, ''), frozenset(), Fraction(1))


@dataclass(frozen=True)
class Overlap:

    # A negative integer.
    # Represents the shift value for left contig relative to right
    # contig to achieve the maximum overlap.
    shift: int


def get_overlap(left: Contig, right: Contig) -> Optional[Overlap]:
    shift = find_maximum_overlap(left.seq, right.seq)
    if shift == 0:
        return None

    return Overlap(shift)


def combine_contigs(a: Contig,
                    b: Contig,
                    overlap: Overlap,
                    ) -> Tuple[Contig, Fraction]:
    # FIXME: This is a trivial implementation. It often doesn't work.
    #        It is easy to improve.

    if abs(overlap.shift) < len(a.seq):
        a_goes_first = True
    else:
        a_goes_first = False

    (left, right) = (a, b) if a_goes_first else (b, a)
    adjusted_left = left.seq[:-abs(overlap.shift)]
    adjusted_right = right.seq[abs(overlap.shift):]
    assert adjusted_left is not None
    assert adjusted_right is not None

    result_seq = left.seq + adjusted_right
    result_contig = Contig(None, result_seq)

    result_probability = Fraction(1, 2)

    return (result_contig, result_probability)


def try_combine_contigs(a: Contig, b: Contig,
                        ) -> Optional[Tuple[Contig, Fraction]]:
    # TODO: Memoize this function.
    #       Two-layer caching seems most optimal:
    #       first by key=contig.id, then by key=contig.seq.

    overlap = get_overlap(a, b)
    if overlap is None:
        return None

    return combine_contigs(a, b, overlap)


def extend_by_1(path: ContigsPath, candidate: Contig) -> Iterator[ContigsPath]:
    if path.has_contig(candidate):
        return

    if path.is_empty:
        (combined, prob) = (candidate, Fraction(1))
    else:
        combination = try_combine_contigs(path.whole, candidate)
        if combination is None:
            return
        (combined, prob) = combination

    score = path.score * prob
    new_elements = path.parts_ids.union([candidate.id])
    new_path = ContigsPath(combined, new_elements, score)
    yield new_path


def calc_extension(contigs: Sequence[Contig],
                   path: ContigsPath,
                   ) -> Iterator[ContigsPath]:

    for contig in contigs:
        yield from extend_by_1(path, contig)


def calc_multiple_extensions(paths: Iterable[ContigsPath],
                             contigs: Sequence[Contig],
                             ) -> Iterator[ContigsPath]:
    for path in paths:
        yield from calc_extension(contigs, path)


def calculate_all_paths(contigs: Sequence[Contig]) -> Iterator[ContigsPath]:
    initial = ContigsPath.empty()
    paths: Iterable[ContigsPath] = [initial]
    while paths:
        paths = tuple(calc_multiple_extensions(paths, contigs))
        yield from paths


def find_most_probable_path(contigs: Sequence[Contig]) -> ContigsPath:
    paths = calculate_all_paths(contigs)
    return min(paths, key=lambda path: path.score)


def stitch_consensus(contigs: Iterable[Contig]) -> Iterable[Contig]:
    remaining = tuple(contigs)
    while remaining:
        most_probable = find_most_probable_path(remaining)
        yield most_probable.whole
        remaining = tuple(contig for contig in remaining
                          if not most_probable.has_contig(contig))


def read_referenceless_contigs(input_csv: TextIO) -> Iterable[Contig]:
    for row in csv.DictReader(input_csv):
        seq = row['contig']
        yield Contig(name=None, seq=seq)
