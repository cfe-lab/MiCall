
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

    @classmethod
    @contextmanager
    def fresh(cls) -> Iterator['GenericStitcherContext']:
        ctx = cls.make()
        token = context.set(ctx)
        try:
            with registry.stage():
                yield ctx
        finally:
            context.reset(token)

    @staticmethod
    @contextmanager
    def stage() -> Iterator['GenericStitcherContext']:
        ctx = deepcopy(context.get())
        token = context.set(ctx)
        try:
            with registry.fresh():
                yield ctx
        finally:
            context.reset(token)


class ReferencefullStitcherContext(GenericStitcherContext[full_events.EventType]):
    @staticmethod
    def make() -> 'GenericStitcherContext':
        return ReferencefullStitcherContext([])


class ReferencelessStitcherContext(GenericStitcherContext[less_events.EventType]):
    @staticmethod
    def make() -> 'GenericStitcherContext':
        return ReferencelessStitcherContext([])


context: ContextVar[GenericStitcherContext] = ContextVar("GenericStitcherContext")
