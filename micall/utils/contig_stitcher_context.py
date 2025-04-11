
from abc import ABC, abstractmethod
from dataclasses import dataclass
from typing import List, TypeVar, Generic, Iterator
from contextvars import ContextVar
from contextlib import contextmanager
from copy import deepcopy

import micall.utils.referencefull_contig_stitcher_events as full_events
import micall.utils.referenceless_contig_stitcher_events as less_events
import micall.utils.registry as registry


T = TypeVar('T')


@dataclass(frozen=True)
class GenericStitcherContext(ABC, Generic[T]):
    events: List[T]

    def emit(self, event: T) -> None:
        self.events.append(event)

    @staticmethod
    @abstractmethod
    def make() -> 'GenericStitcherContext': ...

    @staticmethod
    @abstractmethod
    def _context() -> ContextVar: ...

    @classmethod
    @contextmanager
    def fresh(cls) -> Iterator['GenericStitcherContext']:
        ctx = cls.make()
        class_context = cls._context()
        token = class_context.set(ctx)
        try:
            with registry.fresh():
                yield ctx
        finally:
            class_context.reset(token)

    @classmethod
    @contextmanager
    def stage(cls) -> Iterator['GenericStitcherContext']:
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
    def make() -> 'ReferencefullStitcherContext':
        return ReferencefullStitcherContext([])

    @staticmethod
    def get() -> 'ReferencefullStitcherContext':
        return _referencefull_context.get()

    @staticmethod
    def set(value: 'ReferencefullStitcherContext') -> None:
        _referencefull_context.set(value)

    @staticmethod
    def _context() -> ContextVar['ReferencefullStitcherContext']:
        return _referencefull_context


class ReferencelessStitcherContext(GenericStitcherContext[less_events.EventType]):
    @staticmethod
    def make() -> 'ReferencelessStitcherContext':
        return ReferencelessStitcherContext([])

    @staticmethod
    def get() -> 'ReferencelessStitcherContext':
        return _referenceless_context.get()

    @staticmethod
    def set(value: 'ReferencelessStitcherContext') -> None:
        _referenceless_context.set(value)

    @staticmethod
    def _context() -> ContextVar['ReferencelessStitcherContext']:
        return _referenceless_context


_referencefull_context: ContextVar[ReferencefullStitcherContext] = ContextVar("ReferencefullStitcherContext")
_referenceless_context: ContextVar[ReferencelessStitcherContext] = ContextVar("ReferencelessStitcherContext")
