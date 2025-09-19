from abc import ABC, abstractmethod
from typing import List, TypeVar, Generic, Iterator, Self, Dict, Set, Tuple, Any
from contextvars import ContextVar
from contextlib import contextmanager
from copy import deepcopy

import micall.utils.referencefull_contig_stitcher_events as full_events
import micall.utils.referenceless_contig_stitcher_events as less_events
import micall.utils.registry as registry


T = TypeVar('T')


class GenericStitcherContext(ABC, Generic[T]):
    def __init__(self):
        self.events: List[T] = []

    def emit(self, event: T) -> None:
        self.events.append(event)

    @classmethod
    @abstractmethod
    def _context(cls) -> ContextVar[Self]: ...

    @classmethod
    def get(cls) -> Self:
        return cls._context().get()

    @classmethod
    def set(cls, value: Self) -> None:
        cls._context().set(value)

    @classmethod
    @contextmanager
    def fresh(cls) -> Iterator[Self]:
        ctx = cls()
        class_context = cls._context()
        token = class_context.set(ctx)
        try:
            with registry.fresh():
                yield ctx
        finally:
            class_context.reset(token)

    @classmethod
    @contextmanager
    def stage(cls) -> Iterator[Self]:
        class_context = cls._context()
        try:
            existing = class_context.get()
        except BaseException:
            with cls.fresh() as ret:
                yield ret
                return

        ctx = deepcopy(existing)
        token = class_context.set(ctx)
        try:
            with registry.stage():
                yield ctx
        finally:
            class_context.reset(token)


class ReferencefullStitcherContext(GenericStitcherContext[full_events.EventType]):
    @staticmethod
    def _context() -> ContextVar['ReferencefullStitcherContext']:
        return _referencefull_context


class ReferencelessStitcherContext(GenericStitcherContext[less_events.EventType]):
    def __init__(self) -> None:
        # debug flags
        self.is_debug2: bool = False
        # overlap detection caches: key=(left_id,right_id) -> Overlap (typed as Any to avoid circular import)
        self.get_overlap_cache: Dict[Tuple[int, int], Any] = {}
        self.get_overlap_negative: Set[Tuple[int, int]] = set()
        # alignment cache for overlap windows: key=(left_overlap,right_overlap) -> (aligned_left, aligned_right)
        self.align_cache: Dict[Tuple[str, str], Tuple[str, str]] = {}
        # cutoffs caches: key=(left_id,right_id) -> (left_cutoff,right_cutoff)
        self.cutoffs_cache: Dict[Tuple[int, int], Tuple[int, int]] = {}
        self.cutoffs_negative: Set[Tuple[int, int]] = set()
        super().__init__()

    @staticmethod
    def _context() -> ContextVar['ReferencelessStitcherContext']:
        return _referenceless_context


_referencefull_context: ContextVar[ReferencefullStitcherContext] = ContextVar("ReferencefullStitcherContext")
_referenceless_context: ContextVar[ReferencelessStitcherContext] = ContextVar("ReferencelessStitcherContext")
