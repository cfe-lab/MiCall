from typing import Iterable, Iterator, Optional, Tuple, \
    Sequence, TextIO, MutableMapping, Literal, Union
from dataclasses import dataclass
from Bio import SeqIO, Seq
from Bio.SeqRecord import SeqRecord
import logging
from sortedcontainers import SortedList
import itertools

from micall.utils.contig_stitcher_context import context, ReferencelessStitcherContext
from micall.utils.overlap_stitcher import align_queries, \
    calculate_concordance, sort_concordance_indexes, calc_overlap_pvalue
import micall.utils.referenceless_contig_stitcher_events as events
from micall.utils.referenceless_contig_with_aligner import ContigWithAligner
from micall.utils.referenceless_contig_path import ContigsPath
from micall.utils.referenceless_score import Score, SCORE_EPSILON, SCORE_NOTHING
from micall.utils.contig_stitcher_contigs import Contig


logger = logging.getLogger(__name__)


ACCEPTABLE_STITCHING_SCORE: Score = calc_overlap_pvalue(L=15, M=15)
MAX_ALTERNATIVES = 30


def log(e: events.EventType) -> None:
    context.get().emit(e)
    logger.debug("%s", e)


@dataclass
class Pool:
    # TODO: Factor out a SortedRing structure out of here.
    paths: SortedList[ContigsPath]
    existing: MutableMapping[str, ContigsPath]
    size: int
    capacity: int
    smallest_score: Score

    @staticmethod
    def empty() -> 'Pool':
        return Pool(SortedList(), {}, 0, 999999, ACCEPTABLE_STITCHING_SCORE)

    @property
    def min_acceptable_score(self) -> Score:
        return self.smallest_score

    def resize(self, new_capacity: int) -> None:
        if new_capacity < self.size:
            stop = self.size - new_capacity
            del self.paths[:stop]
            self.size = new_capacity

        self.capacity = new_capacity

        if self.size > 0:
            smallest_path = self.paths[0]
            self.smallest_score = max(ACCEPTABLE_STITCHING_SCORE,
                                      smallest_path.probability)
        else:
            self.smallest_score = ACCEPTABLE_STITCHING_SCORE

    def add(self, path: ContigsPath) -> bool:
        key = path.whole.seq
        alternative = self.existing.get(key)
        if alternative is not None and alternative.score() >= path.score():
            return False

        if self.size > 0:
            smallest_path = self.paths[0]
            if self.size >= self.capacity:
                if smallest_path.score() >= path.score():
                    return False

                del self.paths[0]
                self.size -= 1

            if path.score() < smallest_path.score():
                smallest_path = path

            self.smallest_score = max(ACCEPTABLE_STITCHING_SCORE,
                                      smallest_path.probability)
        else:
            self.smallest_score = max(ACCEPTABLE_STITCHING_SCORE,
                                      path.score())

        self.size += 1
        self.paths.add(path)
        self.existing[key] = path
        return True


@dataclass(frozen=True)
class Overlap:
    # A negative integer.
    # Represents the shift value for left contig relative to right
    # contig to achieve the maximum overlap.
    shift: int
    size: int


ContigId = int
GET_OVERLAP_CACHE: MutableMapping[
    Tuple[ContigId, ContigId],
    Optional[Overlap]] = {}


def get_overlap(left: ContigWithAligner, right: ContigWithAligner) -> Optional[Overlap]:
    if len(left.seq) == 0 or len(right.seq) == 0:
        return None

    key = (left.id, right.id)
    cached: Union[Overlap, None, Literal[0]] = GET_OVERLAP_CACHE.get(key, 0)
    if cached != 0:
        return cached

    (shift, value) = left.find_maximum_overlap(right)
    if shift == 0:
        ret = None
        GET_OVERLAP_CACHE[key] = ret
        return ret

    if abs(shift) <= len(left.seq):
        size = min(abs(shift), len(right.seq))
    else:
        size = min(len(left.seq), len(left.seq) + len(right.seq) + shift)

    assert size > 0, f"{shift}, {len(left.seq)}, {len(right.seq)}"
    assert size <= len(left.seq), f"{shift}, {size}, {len(left.seq)}, {len(right.seq)}"
    assert size <= len(right.seq), f"{shift}, {size}, {len(left.seq)}, {len(right.seq)}"

    ret = Overlap(shift=shift, size=size)
    GET_OVERLAP_CACHE[key] = ret
    return ret


def combine_probability(current: Score, new: Score) -> Score:
    return current + new


def get_minimum_score(current: Score, minimum: Score) -> Score:
    """
    Returns minimum probability that is still acceptable.
    """

    return minimum - current


ALIGN_CACHE: MutableMapping[Tuple[str, str], Tuple[str, str]] = {}


def align_overlaps(left_overlap: str, right_overlap: str) -> Tuple[str, str]:
    key = (left_overlap, right_overlap)
    existing = ALIGN_CACHE.get(key)
    if existing is not None:
        return existing

    result = align_queries(left_overlap, right_overlap)
    ALIGN_CACHE[key] = result
    return result


CutoffsCacheResult = Optional[Tuple[int, int]]
CutoffsCache = MutableMapping[Tuple[ContigId, ContigId], CutoffsCacheResult]
CUTOFFS_CACHE: CutoffsCache = {}


def find_overlap_cutoffs(minimum_score: Score,
                         left: ContigWithAligner,
                         right: ContigWithAligner,
                         left_initial_overlap: str,
                         right_initial_overlap: str,
                         ) -> CutoffsCacheResult:

    key = (left.id, right.id)

    existing: Union[CutoffsCacheResult, Literal[False]] \
        = CUTOFFS_CACHE.get(key, False)

    if existing is None:
        return None

    elif existing is not False:
        (left_cutoff, right_cutoff) = existing

    elif len(left.seq) < len(right.seq):
        left_overlap_alignments = left.map_overlap(minimum_score, True, right_initial_overlap)
        left_cutoff = min((al.r_st + shift for al, shift in left_overlap_alignments), default=-1)
        if left_cutoff < 0:
            CUTOFFS_CACHE[key] = None
            return None

        right_overlap_alignments = right.map_overlap(minimum_score, False, left_initial_overlap)
        right_cutoff = max((al.r_en + shift for al, shift in right_overlap_alignments), default=-1)
        if right_cutoff < 0:
            CUTOFFS_CACHE[key] = None
            return None
    else:
        right_overlap_alignments = right.map_overlap(minimum_score, False, left_initial_overlap)
        right_cutoff = max((al.r_en + shift for al, shift in right_overlap_alignments), default=-1)
        if right_cutoff < 0:
            CUTOFFS_CACHE[key] = None
            return None

        left_overlap_alignments = left.map_overlap(minimum_score, True, right_initial_overlap)
        left_cutoff = min((al.r_st + shift for al, shift in left_overlap_alignments), default=-1)
        if left_cutoff < 0:
            CUTOFFS_CACHE[key] = None
            return None

    ret = (left_cutoff, right_cutoff)
    CUTOFFS_CACHE[key] = ret
    return ret


def try_combine_contigs(current_prob: Score,
                        pool: Pool,
                        a: ContigWithAligner, b: ContigWithAligner,
                        ) -> Optional[Tuple[ContigWithAligner, Score]]:

    minimum_score = get_minimum_score(current_prob, pool.min_acceptable_score)

    maximum_overlap_size = min(len(a.seq), len(b.seq)) - 1
    maximum_number_of_matches = maximum_overlap_size
    maximum_result_probability = calc_overlap_pvalue(L=maximum_overlap_size, M=maximum_number_of_matches)
    if maximum_result_probability < minimum_score:
        return None

    overlap = get_overlap(a, b)
    if overlap is None:
        return None

    optimistic_overlap_size = overlap.size
    optimistic_number_of_matches = optimistic_overlap_size
    optimistic_result_probability = calc_overlap_pvalue(L=optimistic_overlap_size, M=optimistic_number_of_matches)
    if optimistic_result_probability < minimum_score:
        return None

    if abs(overlap.shift) < len(a.seq):
        shift = overlap.shift
        left = a
        right = b
    else:
        shift = (len(b.seq) + len(a.seq) - abs(overlap.shift)) * -1
        left = b
        right = a

    left_initial_overlap = left.seq[len(left.seq) - abs(shift):(len(left.seq) - abs(shift) + len(right.seq))]
    right_initial_overlap = right.seq[:abs(shift)]

    assert len(right_initial_overlap) == overlap.size
    assert len(left_initial_overlap) == overlap.size, f"{len(left_initial_overlap)} == {overlap.size}"

    assert calc_overlap_pvalue(L=len(left_initial_overlap), M=len(left_initial_overlap)) >= minimum_score
    assert len(right_initial_overlap) == overlap.size, f"{len(right_initial_overlap)} == {overlap.size} in {len(right.seq), {shift}}"
    assert len(left_initial_overlap) == len(right_initial_overlap), f"{len(left_initial_overlap)} == {len(right_initial_overlap)}"

    cutoffs = find_overlap_cutoffs(minimum_score,
                                   left, right,
                                   left_initial_overlap,
                                   right_initial_overlap,
                                   )
    if cutoffs is None:
        return None

    (left_cutoff, right_cutoff) = cutoffs
    left_overlap = left.seq[left_cutoff:(left_cutoff + len(right.seq))]
    left_remainder = left.seq[:left_cutoff]
    right_overlap = right.seq[:right_cutoff]
    right_remainder = right.seq[right_cutoff:]
    aligned_left, aligned_right = align_overlaps(left_overlap, right_overlap)

    number_of_matches = sum(1 for x, y
                            in zip(aligned_left, aligned_right)
                            if x == y and x != '-')
    result_probability = calc_overlap_pvalue(L=len(left_overlap), M=number_of_matches)
    if result_probability < minimum_score:
        return None

    is_covered = len(right.seq) < abs(shift)
    if is_covered:
        log(events.Covered(left.unique_name, right.unique_name, cutoffs=cutoffs))
        return (left, SCORE_EPSILON)

    else:
        concordance = calculate_concordance(aligned_left, aligned_right)
        max_concordance_index = next(iter(sort_concordance_indexes(concordance)))
        left_overlap_chunk = ''.join(x for x in aligned_left[:max_concordance_index] if x != '-')
        right_overlap_chunk = ''.join(x for x in aligned_right[max_concordance_index:] if x != '-')

        result_seq = left_remainder + left_overlap_chunk + right_overlap_chunk + right_remainder
        result_contig = ContigWithAligner(None, result_seq)

        log(events.FoundOverlap(left_contig=left.unique_name,
                                right_contig=right.unique_name,
                                result_contig=result_contig.unique_name,
                                overlap_size=len(aligned_left),
                                cutoffs=cutoffs,
                                ))
        return (result_contig, result_probability)


def extend_by_1(pool: Pool,
                path: ContigsPath,
                candidate: ContigWithAligner,
                ) -> bool:
    if path.has_contig(candidate):
        return False

    combination = try_combine_contigs(path.probability, pool, path.whole, candidate)
    if combination is None:
        return False

    (combined, prob) = combination
    probability = combine_probability(path.probability, prob)
    new_elements = path.parts_ids.union([candidate.id])
    new_path = ContigsPath(combined, new_elements, probability)
    return pool.add(new_path)


def calc_extension(pool: Pool,
                   contigs: Sequence[ContigWithAligner],
                   path: ContigsPath,
                   ) -> bool:
    ret = False
    for contig in contigs:
        ret = extend_by_1(pool, path, contig) or ret
    return ret


def calc_multiple_extensions(pool: Pool,
                             paths: Iterable[ContigsPath],
                             contigs: Sequence[ContigWithAligner],
                             ) -> bool:
    ret = False
    for path in paths:
        ret = calc_extension(pool, contigs, path) or ret
    return ret


def calculate_all_paths(paths: Sequence[ContigsPath],
                        contigs: Sequence[ContigWithAligner],
                        ) -> Iterable[ContigsPath]:
    pool = Pool.empty()
    pool.resize(MAX_ALTERNATIVES)

    for path in sorted(paths):
        if pool.size >= pool.capacity:
            break
        pool.add(path)

    log(events.CalculatingAll())
    for cycle in itertools.count(1):
        log(events.CycleStart(cycle, pool.size))

        if not calc_multiple_extensions(pool, pool.paths, contigs):
            break

        log(events.CycleEnd(cycle, pool.size, pool))

    return pool.paths


def find_most_probable_path(seeds: Sequence[ContigsPath],
                            contigs: Sequence[ContigWithAligner],
                            ) -> ContigsPath:
    paths = calculate_all_paths(seeds, contigs)
    return max(paths, key=lambda path: path.score())


def contig_size_fun(contig: Contig) -> int:
    return -len(contig.seq)


def find_best_candidates(first: ContigWithAligner,
                         contigs: Iterable[ContigWithAligner],
                         ) -> Iterator[ContigsPath]:

    initial = ContigsPath.singleton(first)
    yield initial

    pool = Pool.empty()
    pool.resize(1)

    for candidate in contigs:
        if first.id >= candidate.id:
            continue

        extend_by_1(pool=pool,
                    path=initial,
                    candidate=candidate,
                    )

    yield from pool.paths


def stitch_consensus_overlaps(contigs: Iterable[ContigWithAligner]) -> Iterator[ContigWithAligner]:
    remaining = tuple(sorted(contigs, key=contig_size_fun))

    log(events.InitializingSeeds())
    seeds = tuple(sorted(map(ContigsPath.singleton, remaining),
                         key=lambda path: len(path.whole.seq),
                         reverse=True))
    log(events.Starting(len(seeds)))

    while remaining:
        most_probable = find_most_probable_path(seeds, remaining)
        log(events.Constructed(most_probable))
        yield most_probable.whole
        remaining = tuple(contig for contig in remaining
                          if not most_probable.has_contig(contig))
        seeds = tuple(path for path in seeds
                      if most_probable.parts_ids.isdisjoint(path.parts_ids))
        log(events.Remove(len(most_probable.parts_ids), len(remaining)))
        if len(most_probable.parts_ids) == 1:
            log(events.GiveUp())
            yield from remaining
            return


def try_combine_1(contigs: Iterable[ContigWithAligner],
                  ) -> Optional[Tuple[ContigWithAligner,
                                      ContigWithAligner,
                                      ContigWithAligner]]:
    for first in contigs:
        for second in contigs:
            if first.id >= second.id:
                continue

            pool = Pool.empty()
            result = try_combine_contigs(current_prob=SCORE_NOTHING,
                                         pool=pool,
                                         a=first, b=second,
                                         )
            if result is not None:
                (combined, prob) = result
                return (first, second, combined)

    return None


def o2_loop(contigs: Iterable[ContigWithAligner],
            ) -> Iterator[ContigWithAligner]:
    buf = {contig.id: contig for contig in contigs}

    while True:
        combination = try_combine_1(buf.values())
        if combination is None:
            break

        (first, second, combined) = combination
        del buf[first.id]
        del buf[second.id]
        buf[combined.id] = combined

    yield from sorted(buf.values(), key=contig_size_fun)


def stitch_consensus(contigs: Iterable[ContigWithAligner],
                     ) -> Iterator[ContigWithAligner]:
    contigs = tuple(stitch_consensus_overlaps(contigs))
    log(events.InitiallyProduced(len(contigs)))
    contigs = o2_loop(contigs)
    yield from contigs


def write_contigs(output_fasta: TextIO, contigs: Iterable[ContigWithAligner]):
    records = (SeqRecord(Seq.Seq(contig.seq),
                         description='',
                         id=contig.unique_name,
                         name=contig.unique_name)
               for contig in contigs)
    SeqIO.write(records, output_fasta, "fasta")


def read_contigs(input_fasta: TextIO) -> Iterable[ContigWithAligner]:
    for record in SeqIO.parse(input_fasta, "fasta"):
        yield ContigWithAligner(name=record.name, seq=str(record.seq))


def referenceless_contig_stitcher_with_ctx(
        input_fasta: TextIO,
        output_fasta: Optional[TextIO],
) -> int:
    contigs = tuple(read_contigs(input_fasta))
    log(events.Loaded(len(contigs)))

    if output_fasta is not None:
        contigs = tuple(stitch_consensus(contigs))
        log(events.Outputting(len(contigs)))

    if output_fasta is not None:
        write_contigs(output_fasta, contigs)

    return len(contigs)


def referenceless_contig_stitcher(input_fasta: TextIO,
                                  output_fasta: Optional[TextIO],
                                  ) -> int:
    with ReferencelessStitcherContext.fresh():
        return referenceless_contig_stitcher_with_ctx(input_fasta, output_fasta)
