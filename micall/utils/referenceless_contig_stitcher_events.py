from typing import Union
from dataclasses import dataclass

from micall.utils.referenceless_contig_path import ContigsPath


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
        return f"Cycle {self.i_cycle} started with {self.n_paths} paths"


class CycleEnd:
    def __init__(self, i_cycle: int, n_paths: int, pool: object) -> None:
        self.i_cycle: int = i_cycle
        self.n_paths: int = n_paths
        fittest = pool.paths[-1]  # type: ignore
        self.length = len(fittest.whole.seq)  # type: ignore
        self.parts = len(fittest.parts_ids)  # type: ignore

    def __str__(self) -> str:
        return f"Cycle {self.i_cycle} finished with {self.n_paths} new paths. " \
            f"The fittest has length {self.length} and consists of {self.parts} parts."


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


EventType = Union[GiveUp, Remove, CalculatingAll, CycleStart,
                  CycleEnd, InitializingSeeds, Starting, Constructed]
