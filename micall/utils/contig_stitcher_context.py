from typing import List, TypeVar, Generic
from contextvars import ContextVar
from contextlib import contextmanager
from dataclasses import dataclass
from copy import deepcopy

import micall.utils.referencefull_contig_stitcher_events as full_events
import micall.utils.registry as registry


T = TypeVar('T')


@dataclass(frozen=True)
class GenericStitcherContext(Generic[T]):
    events: List[T]

    def emit(self, event: T) -> None:
        self.events.append(event)

    @staticmethod
    def make() -> 'GenericStitcherContext':
        return GenericStitcherContext(events=[])

    @staticmethod
    @contextmanager
    def fresh():
        ctx = GenericStitcherContext.make()
        token = context.set(ctx)
        try:
            with registry.ensure():
                yield ctx
        finally:
            context.reset(token)

    @staticmethod
    @contextmanager
    def stage():
        ctx = deepcopy(context.get())
        token = context.set(ctx)
        try:
            with registry.fresh():
                yield ctx
        finally:
            context.reset(token)


ReferencefullStitcherContext = GenericStitcherContext[full_events.EventType]
ReferencelessStitcherContext = GenericStitcherContext[object]

context: ContextVar[GenericStitcherContext] = ContextVar("GenericStitcherContext")
