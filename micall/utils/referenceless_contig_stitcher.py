from typing import Iterable, Iterator, Optional, FrozenSet, Tuple, Sequence, TextIO
from dataclasses import dataclass
from fractions import Fraction
from Bio import SeqIO, Seq
from Bio.SeqRecord import SeqRecord
from micall.utils.contig_stitcher_context import StitcherContext
from gotoh import align_it

from micall.utils.contig_stitcher_contigs import Contig
from micall.utils.find_maximum_overlap import find_maximum_overlap


ACCEPTABLE_STITCHING_PROB = Fraction(1, 100)


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
        return ContigsPath(Contig.empty(), frozenset(), Fraction(1))


@dataclass(frozen=True)
class Overlap:

    # A negative integer.
    # Represents the shift value for left contig relative to right
    # contig to achieve the maximum overlap.
    shift: int


def get_overlap(left: Contig, right: Contig) -> Optional[Overlap]:
    if len(left.seq) == 0 or len(right.seq) == 0:
        return None

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
        shift = overlap.shift
        left = a
        right = b
    else:
        shift = (len(b.seq) + len(a.seq) - abs(overlap.shift)) * -1
        left = b
        right = a

    left_initial_overlap = left.seq[len(left.seq) - abs(shift):]
    left_initial_remainder = left.seq[:len(left.seq) - abs(shift)]
    right_initial_overlap = right.seq[:abs(shift)]
    right_initial_remainder = right.seq[abs(shift):]

    gap_open_penalty = 15
    gap_extend_penalty = 3
    use_terminal_gap_penalty = 1
    aligned_left, aligned_right, alignment_score = align_it(
        str(left_initial_overlap),
        str(right_initial_overlap),
        gap_open_penalty,
        gap_extend_penalty,
        use_terminal_gap_penalty)

    resulting_matches = sum(1 for x, y
                            in zip(aligned_left, aligned_right)
                            if x == y and x != '-')

    result_seq = left_initial_remainder + left_initial_overlap + right_initial_remainder
    result_contig = Contig(None, result_seq)

    # FIXME: Calculate a more accurate probability.
    result_probability = Fraction(1.0) - Fraction(resulting_matches, len(left_initial_overlap) + 2)
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
        if prob > ACCEPTABLE_STITCHING_PROB:
            # FIXME: Adjust the threold to something more based.
            # FIXME: Print a warning that this happened.
            return

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


def write_contigs(output_fasta: TextIO, contigs: Iterable[Contig]):
    with StitcherContext.fresh():
        records = (SeqRecord(Seq.Seq(contig.seq),
                             description='',
                             id=contig.unique_name,
                             name=contig.unique_name)
                   for contig in contigs)
        SeqIO.write(records, output_fasta, "fasta")


def read_contigs(input_fasta: TextIO) -> Iterable[Contig]:
    for record in SeqIO.parse(input_fasta, "fasta"):
        yield Contig(name=record.name, seq=record.seq)


def referenceless_contig_stitcher(input_fasta: TextIO,
                                  output_fasta: Optional[TextIO],
                                  ) -> int:
    contigs = list(read_contigs(input_fasta))

    if output_fasta is not None:
        contigs = list(stitch_consensus(contigs))

    if output_fasta is not None:
        write_contigs(output_fasta, contigs)

    return len(contigs)