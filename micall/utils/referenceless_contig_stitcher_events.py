from typing import Union, Tuple, Optional
from dataclasses import dataclass

from micall.utils.referenceless_contig_path import ContigsPath
from micall.utils.referenceless_score import Score


@dataclass(frozen=True)
class GiveUp:
    def __str__(self) -> str:
        return "Giving up on attempts to stitch more overlaps since most probable is a singleton."


@dataclass(frozen=True)
class Remove:
    n_removed: int
    n_remaining: int

    def __str__(self) -> str:
        return f"Removed {self.n_removed} components from the working list, having {self.n_remaining} still to process."


@dataclass(frozen=True)
class CalculatingAll:
    def __str__(self) -> str:
        return "Calculating all paths..."


@dataclass(frozen=True)
class CycleStart:
    i_cycle: int
    n_paths: int

    def __str__(self) -> str:
        return f"Cycle {self.i_cycle} started with {self.n_paths} paths."


class CycleEnd:
    def __init__(self, i_cycle: int, n_paths: int, pool: object) -> None:
        self.i_cycle: int = i_cycle
        self.n_paths: int = n_paths
        self.pool = pool

    def __str__(self) -> str:
        fittest = self.pool.ring[-1]  # type: ignore
        length = len(fittest.whole.seq)  # type: ignore
        parts = len(fittest.contigs_ids)  # type: ignore

        return f"Cycle {self.i_cycle} finished with {self.n_paths} new paths. " \
            f"The fittest has length {length} and consists of {parts} parts."


@dataclass(frozen=True)
class InitializingSeeds:
    def __str__(self) -> str:
        return "Initializing initial seeds..."


@dataclass(frozen=True)
class Starting:
    n_seeds: int

    def __str__(self) -> str:
        return f"Starting with {self.n_seeds} initial seeds."


@dataclass(frozen=True)
class Constructed:
    path: ContigsPath

    def __str__(self) -> str:
        return f"Constructed a path of length {len(self.path.whole.seq)}."


@dataclass(frozen=True)
class Loaded:
    n_contigs: int

    def __str__(self) -> str:
        return f"Loaded {self.n_contigs} contigs."


@dataclass(frozen=True)
class Outputting:
    n_contigs: int

    def __str__(self) -> str:
        return f"Outputting {self.n_contigs} contigs."


@dataclass(frozen=True)
class InitiallyProduced:
    n_contigs: int

    def __str__(self) -> str:
        return f"Initial overlap stitching produced {self.n_contigs} contigs."


@dataclass(frozen=True)
class Covered:
    left_contig: str
    right_contig: str

    def __str__(self) -> str:
        return f"Contig {self.left_contig} covers contig {self.right_contig} completely."


@dataclass(frozen=True)
class CombinedContigs:
    left_contig: str
    right_contig: str
    result_contig: str
    overlap_size: int

    def __str__(self) -> str:
        return f"Found a significant overlap of size {self.overlap_size}" \
            f" between contigs {self.left_contig} and {self.right_contig}," \
            f" resulting in contig {self.result_contig}."


@dataclass(frozen=True)
class CalculatedCutoffs:
    left_contig: str
    right_contig: str
    overlap_size: int
    cutoffs: Optional[Tuple[int, int]]

    def __str__(self) -> str:
        return f"Calculated cutoff for an overlap of size {self.overlap_size}" \
            f" between contigs {self.left_contig} and {self.right_contig}" \
            f" to be {self.cutoffs}."


@dataclass(frozen=True)
class DeterminedOverlap:
    left_contig: str
    right_contig: str
    aligned_size: int
    number_of_matches: int
    relative_score: Score

    def __str__(self) -> str:
        return f"Overlap between contigs {self.left_contig} and {self.right_contig}" \
            f" has aligned size {self.aligned_size}, {self.number_of_matches} matches," \
            f" and the score of {self.relative_score}."


EventType = Union[GiveUp, Remove, CalculatingAll, CycleStart,
                  CycleEnd, InitializingSeeds, Starting, Constructed,
                  Loaded, Outputting, InitiallyProduced, Covered,
                  CombinedContigs, CalculatedCutoffs, DeterminedOverlap]
