"""
Module for handling CIGAR strings and related alignment formats.
"""

from math import ceil, floor
import re
from typing import Container, Tuple, Iterable, Optional, Set
from dataclasses import dataclass
from functools import cached_property
from itertools import chain, dropwhile

from micall.utils.consensus_aligner import CigarActions


CIGAR_OP_MAPPING = {
    'M': CigarActions.MATCH,
    'I': CigarActions.INSERT,
    'D': CigarActions.DELETE,
    'N': CigarActions.SKIPPED,
    'S': CigarActions.SOFT_CLIPPED,
    'H': CigarActions.HARD_CLIPPED,
    'P': CigarActions.PADDING,
    '=': CigarActions.SEQ_MATCH,
    'X': CigarActions.MISMATCH,
}


def parse_cigar_operation(operation: str) -> CigarActions:
    if operation in CIGAR_OP_MAPPING:
        return CIGAR_OP_MAPPING[operation]
    else:
        raise ValueError(f"Unexpected CIGAR action: {operation}.")


def cigar_operation_to_str(op: CigarActions) -> str:
    return [k for (k, v) in CIGAR_OP_MAPPING.items() if v == op][0]


class CoordinateMapping:
    def __init__(self):
        self.query_to_ref_d = {}
        self.ref_to_query_d = {}
        self.ref_to_op_d = {}
        self.query_to_op_d = {}


    def extend(self,
               ref_index: Optional[int],
               query_index: Optional[int],
               op_index: Optional[int]):

        if ref_index is not None and query_index is not None:
            self.ref_to_query_d[ref_index] = query_index
            self.query_to_ref_d[query_index] = ref_index

        if op_index is not None:
            if ref_index is not None:
                self.ref_to_op_d[ref_index] = op_index
            if query_index is not None:
                self.query_to_op_d[query_index] = op_index


    def mapped_reference_coordinates(self) -> Set[int]:
        return set(self.ref_to_query_d.keys())


    def all_reference_coordinates(self) -> Set[int]:
        return set(self.ref_to_op_d.keys())


    def mapped_query_coordinates(self) -> Set[int]:
        return set(self.query_to_ref_d.keys())


    def all_query_coordinates(self) -> Set[int]:
        return set(self.query_to_op_d.keys())


    def ref_to_query(self, index) -> Optional[int]:
        return self.ref_to_query_d.get(index, None)


    def query_to_ref(self, index) -> Optional[int]:
        return self.query_to_ref_d.get(index, None)


    @staticmethod
    def _find_closest_key(mapping: dict, index: int) -> int:
        return min(mapping, key=lambda k: abs(mapping[k] - index))


    def ref_to_closest_query(self, index) -> int:
        return CoordinateMapping._find_closest_key(self.query_to_ref_d, index)


    def query_to_closest_ref(self, index) -> int:
        return CoordinateMapping._find_closest_key(self.ref_to_query_d, index)


    def ref_to_leftsup_query(self, index) -> Optional[int]:
        left_neihbourhood = (k for (k, v) in self.query_to_ref_d.items() if v <= index)
        return max(left_neihbourhood, default=None)


    def ref_to_rightinf_query(self, index) -> Optional[int]:
        right_neihbourhood = (k for (k, v) in self.query_to_ref_d.items() if index <= v)
        return min(right_neihbourhood, default=None)


    def ref_or_query_to_op(self, ref_index: int, query_index: int, conflict):
        r = self.ref_to_op_d.get(ref_index, None)
        q = self.query_to_op_d.get(query_index, None)
        if r is not None and q is not None:
            return conflict(r, q)

        return r if q is None else q


    def translate(self, reference_delta: int, query_delta: int) -> 'CoordinateMapping':
        ret = CoordinateMapping()

        ret.ref_to_query_d = {k + reference_delta: v + query_delta for (k, v) in self.ref_to_query_d.items()}
        ret.query_to_ref_d = {k + query_delta: v + reference_delta for (k, v) in self.query_to_ref_d.items()}
        ret.ref_to_op_d = {k + reference_delta: v for (k, v) in self.ref_to_op_d.items()}
        ret.query_to_op_d = {k + query_delta: v for (k, v) in self.query_to_op_d.items()}

        return ret


class Cigar(list):
    """
    A CIGAR string represents a read alignment against a reference sequence.
    It is a run-length encoded sequence of alignment operations listed below:

       M: Alignment match (can be a sequence match or mismatch)
       D: Deletion from the reference
       I: Insertion to the reference
       S: Soft clip on the read (ignored region, not aligned but present in the read)
       H: Hard clip on the read (ignored region, not present in the read)
       N: Skipped region from the reference
       P: Padding (silent deletion from padded reference, not applicable for our case)
       =: Sequence match
       X: Sequence mismatch

    CIGAR strings are defined in the SAM specification
    (https://samtools.github.io/hts-specs/SAMv1.pdf).
    """


    def __init__(self, cigar_lst: Iterable[Tuple[int, CigarActions]]):
        super().__init__([])
        for x in cigar_lst: self.append(x)


    @staticmethod
    def coerce(obj):
        if isinstance(obj, Cigar):
            return obj

        if isinstance(obj, str):
            return Cigar.parse(obj)

        if isinstance(obj, list):
            return Cigar(obj)

        raise TypeError(f"Cannot coerce {obj!r} to CIGAR string.")


    @staticmethod
    def parse(string):
        data = []
        while string:
            match = re.match(r'([0-9]+)([^0-9])', string)
            if match:
                num, operation = match.groups()
                data.append([int(num), parse_cigar_operation(operation)])
                string = string[match.end():]
            else:
                raise ValueError(f"Invalid CIGAR string. Invalid part: {string[:20]}")

        return Cigar(data)


    def append(self, item: Tuple[int, CigarActions]):
        # Type checking
        if not isinstance(item, list) and not isinstance(item, tuple):
            raise ValueError(f"Invalid CIGAR list: {item!r} is not a tuple.")
        if len(item) != 2:
            raise ValueError(f"Invalid CIGAR list: {item!r} is has a bad length.")

        num, operation = item
        if isinstance(operation, int):
            operation = CigarActions(operation)
        if not isinstance(num, int) or not isinstance(operation, CigarActions):
            raise ValueError(f"Invalid CIGAR list: {item!r} is not a number/operation tuple.")
        if num < 0:
            raise ValueError(f"Invalid CIGAR list: number of operations is negative.")

        # Normalization
        if num == 0:
            return

        if self:
            last_num, last_operation = self[-1]
            if operation == last_operation:
                self[-1] = (last_num + num, operation)
                return

        super().append((num, operation))


    def iterate_operations(self) -> Iterable[CigarActions]:
        for num, operation in self:
            for _ in range(num):
                yield operation


    def iterate_operations_with_pointers(self) -> Iterable[Tuple[CigarActions, Optional[int], Optional[int]]]:
        ref_pointer = 0
        query_pointer = 0

        for operation in self.iterate_operations():
            if operation in (CigarActions.MATCH, CigarActions.SEQ_MATCH, CigarActions.MISMATCH):
                yield (operation, ref_pointer, query_pointer)
                query_pointer += 1
                ref_pointer += 1

            elif operation in (CigarActions.INSERT, CigarActions.SOFT_CLIPPED):
                yield (operation, None, query_pointer)
                query_pointer += 1

            elif operation in (CigarActions.DELETE, CigarActions.SKIPPED):
                yield (operation, ref_pointer, None)
                ref_pointer += 1

            else:
                yield (operation, None, None)


    @cached_property
    def query_length(self):
        return max((query_pointer + 1 if query_pointer is not None else 0 for (_, _, query_pointer)
                    in self.iterate_operations_with_pointers()),
                   default=0)


    @cached_property
    def ref_length(self):
        return max((ref_pointer + 1 if ref_pointer is not None else 0 for (_, ref_pointer, _)
                    in self.iterate_operations_with_pointers()),
                   default=0)


    def slice_operations(self, start_inclusive, end_noninclusive) -> 'Cigar':
        return Cigar([(1, op) for op in self.iterate_operations()]
                     [start_inclusive:end_noninclusive])


    def lstrip(self, actions: Container[CigarActions]) -> 'Cigar':
        """ Return a copy of the Cigar with leading actions removed. """
        return Cigar(dropwhile(lambda tupl: tupl[1] in actions, self))


    def rstrip(self, actions: Container[CigarActions]) -> 'Cigar':
        """ Return a copy of the Cigar with trailing actions removed. """
        return Cigar(reversed(list(dropwhile(lambda tupl: tupl[1] in actions, reversed(self)))))


    @cached_property
    def coordinate_mapping(self) -> CoordinateMapping:
        """
        Convert a CIGAR string to coordinate mapping representing a reference-to-query and query-to-reference coordinate mappings.
        TODO(vitalik): describe the domains and holes.

        :param cigar: a CIGAR string.

        :return: Lists of integers representing the mappings of coordinates from the reference
                sequence to the query sequence, and back.
        """

        mapping = CoordinateMapping()

        for op_pointer, (operation, ref_pointer, query_pointer) in enumerate(self.iterate_operations_with_pointers()):
            mapping.extend(ref_pointer,
                           query_pointer,
                           op_pointer)

        return mapping


    def to_msa(self, reference_seq, query_seq) -> Tuple[str, str]:
        reference_msa = ''
        query_msa = ''

        for operation, ref_pointer, query_pointer in self.iterate_operations_with_pointers():
            if ref_pointer is None and query_pointer is None:
                continue

            try:
                if ref_pointer is not None:
                    reference_msa += reference_seq[ref_pointer]
                else:
                    reference_msa += '-'

                if query_pointer is not None:
                    query_msa += query_seq[query_pointer]
                else:
                    query_msa += '-'

            except IndexError:
                raise ValueError("CIGAR string corresponds to a larger match than either reference or query.")

        return reference_msa, query_msa


    def __repr__(self):
        return f'Cigar({str(self)!r})'


    def __str__(self):
        """ Inverse of Cigar.parse """
        return ''.join('{}{}'.format(num, cigar_operation_to_str(op)) for num, op in self)


@dataclass
class CigarHit:
    cigar: Cigar
    r_st: int
    r_ei: int # inclusive
    q_st: int
    q_ei: int # inclusive


    def __post_init__(self):
        self.cigar = Cigar.coerce(self.cigar)

        if self.ref_length != self.cigar.ref_length:
            raise ValueError(f"CIGAR string maps {self.cigar.ref_length}"
                             f" reference positions, but CIGAR hit range is {self.ref_length}")

        if self.query_length != self.cigar.query_length:
            raise ValueError(f"CIGAR string maps {self.cigar.query_length}"
                             f" query positions, but CIGAR hit range is {self.query_length}")


    @property
    def ref_length(self):
        return self.r_ei + 1 - self.r_st


    @property
    def query_length(self):
        return self.q_ei + 1 - self.q_st


    @staticmethod
    def from_default_alignment(r_st, r_ei, q_st, q_ei):
        ref_length = r_ei - r_st + 1
        query_length = q_ei - q_st + 1
        cigar = Cigar.coerce([[ref_length, CigarActions.DELETE],
                              [query_length, CigarActions.INSERT]])

        return CigarHit(cigar, r_st=r_st, r_ei=r_ei, q_st=q_st, q_ei=q_ei)


    def overlaps(self, other) -> bool:
        """
        Checks if this CIGAR hit overlaps with the other CIGAR hit,
        in either reference or query space.
        NOTE: only applicable if these hits come from the same reference and query.
        """

        def intervals_overlap(x, y):
            """ Check if two intervals [x0, x1] and [y0, y1] overlap. """
            return x[0] <= y[1] and x[1] >= y[0]

        return intervals_overlap((self.r_st, self.r_ei), (other.r_st, other.r_ei)) \
            or intervals_overlap((self.q_st, self.q_ei), (other.q_st, other.q_ei))


    def gaps(self) -> Iterable['CigarHit']:
        # TODO(vitalik): memoize whatever possible.

        covered_coordinates = self.coordinate_mapping.mapped_reference_coordinates()
        all_coordinates = self.coordinate_mapping.all_reference_coordinates()

        def make_gap(r_st, r_en):
            r_ei = r_en - 1
            left, midright = self.cut_reference(r_st - 0.5)
            middle, right = midright.cut_reference(r_ei + 0.5)
            return middle

        gap_start = None
        for coord in all_coordinates:
            if coord in covered_coordinates:
                if gap_start is not None:
                    yield make_gap(gap_start, coord)
                    gap_start = None
            else:
                if gap_start is None:
                    gap_start = coord

        if gap_start is not None:
            yield make_gap(gap_start, coord)


    def __add__(self, other):
        """
        Inserts deletions/insertions between self and other,
        then ajusts boundaries appropriately.
        """

        if self.overlaps(other):
            raise ValueError("Cannot combine overlapping CIGAR hits")

        cigar = self.cigar \
            + CigarHit.from_default_alignment(self.r_ei + 1, other.r_st - 1, self.q_ei + 1, other.q_st - 1).cigar \
            + other.cigar

        return CigarHit(cigar=cigar,
                        r_st=self.r_st,
                        r_ei=other.r_ei,
                        q_st=self.q_st,
                        q_ei=other.q_ei,
                        )


    @property
    def epsilon(self):
        return 1 / (self.query_length + self.ref_length + 100)


    def _slice(self, r_st, r_ei, q_st, q_ei) -> 'CigarHit':
        mapping = self.coordinate_mapping

        o_st = mapping.ref_or_query_to_op(r_st, q_st, min)
        o_ei = mapping.ref_or_query_to_op(r_ei, q_ei, max)
        if o_st is None or o_ei is None:
            cigar = Cigar([])
        else:
            cigar = self.cigar.slice_operations(o_st, o_ei + 1)

        return CigarHit(cigar=cigar,
                        r_st = r_st,
                        r_ei = r_ei,
                        q_st = q_st,
                        q_ei = q_ei,
                        )


    def _ref_cut_to_query_cut(self, cut_point: float):
        mapping = self.coordinate_mapping

        left_query_cut_point = mapping.ref_to_leftsup_query(floor(cut_point))
        right_query_cut_point = mapping.ref_to_rightinf_query(ceil(cut_point))

        if left_query_cut_point is None:
            return self.q_st - 0.1
        if right_query_cut_point is None:
            return self.q_ei + 0.1

        lerp = lambda start, end, t: (1 - t) * start + t * end
        query_cut_point = lerp(left_query_cut_point, right_query_cut_point,
                               cut_point - floor(cut_point))

        if float(query_cut_point).is_integer():
            # Disambiguate to the right.
            query_cut_point += self.epsilon

        return query_cut_point


    def cut_reference(self, cut_point: float) -> 'CigarHit':
        """
        Splits alignment in two parts such that cut_point is in between.
        Guarantees that the two parts do not share any elements,
        and that no element is lost.
        """

        if float(cut_point).is_integer():
            raise ValueError("Cut accepts fractions, not integers")

        if self.ref_length == 0 or \
           not (self.r_st - 1 < cut_point < self.r_ei + 1):
            raise IndexError("Cut point out of reference bounds")

        query_cut_point = self._ref_cut_to_query_cut(cut_point)
        assert (self.q_st - 1 <= query_cut_point <= self.q_ei + 1)

        left = self._slice(self.r_st, floor(cut_point),
                           self.q_st, floor(query_cut_point))
        right = self._slice(ceil(cut_point), self.r_ei,
                            ceil(query_cut_point), self.q_ei)

        return left, right


    def lstrip(self, actions: Container[CigarActions]) -> 'CigarHit':
        """ Return a copy of the CigarHit with leading actions removed. """

        cigar = self.cigar.lstrip(actions)
        reference_delta = cigar.ref_length - self.cigar.ref_length
        query_delta = cigar.query_length - self.cigar.query_length
        return CigarHit(cigar,
                        r_st=self.r_st,
                        r_ei=self.r_ei + reference_delta,
                        q_st=self.q_st,
                        q_ei=self.q_ei + query_delta)


    def rstrip(self, actions: Container[CigarActions]) -> 'CigarHit':
        """ Return a copy of the CigarHit with trailing actions removed. """

        cigar = self.cigar.rstrip(actions)
        reference_delta = cigar.ref_length - self.cigar.ref_length
        query_delta = cigar.query_length - self.cigar.query_length
        return CigarHit(cigar,
                        r_st=self.r_st,
                        r_ei=self.r_ei + reference_delta,
                        q_st=self.q_st,
                        q_ei=self.q_ei + query_delta)


    @cached_property
    def coordinate_mapping(self) -> CoordinateMapping:
        return self.cigar.coordinate_mapping.translate(self.r_st, self.q_st)


    def to_msa(self, reference_seq: str, query_seq: str) -> Tuple[str, str]:
        return self.cigar.to_msa(reference_seq[self.r_st:], query_seq[self.q_st:])


    def translate(self, reference_delta: int, query_delta: int) -> 'CigarHit':
        return CigarHit(cigar=self.cigar,
                        r_st=self.r_st + reference_delta,
                        r_ei=self.r_ei + reference_delta,
                        q_st=self.q_st + query_delta,
                        q_ei=self.q_ei + query_delta)


    def __repr__(self):
        return f'CigarHit({str(self.cigar)!r}, r_st={self.r_st!r}, r_ei={self.r_ei!r}, q_st={self.q_st!r}, q_ei={self.q_ei!r})'


def connect_cigar_hits(cigar_hits: Iterable[CigarHit]) -> CigarHit:
    """
    This function exists to deal with the fact that mappy does not always
    connect big gaps, and returns surrounding parts as two separate alignment hits.

    For those cases we simply connect all the parts that do not overlap.

    Order of cigar_hits matters because we ignore alignments
    that overlap with previously found alignments.
    """

    if not len(cigar_hits) > 0:
        raise ValueError("Expected a non-empty list of cigar hits")

    accumulator = []

    # Collect non-overlaping parts.
    # Earlier matches have priority over ones that come after.
    for hit in cigar_hits:
        if any(earlier.overlaps(hit) for earlier in accumulator):
            continue

        accumulator.append(hit)

    # Sort by interval start positions.
    sorted_parts = sorted(accumulator, key=lambda p: p.r_st)

    # Collect all intervals back together, connecting them with CigarActions.DELETE.
    return sum(sorted_parts[1:], start=sorted_parts[0])
