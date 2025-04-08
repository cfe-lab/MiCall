from typing import List, Dict, TypeVar, Generic
from contextvars import ContextVar
from contextlib import contextmanager
from dataclasses import dataclass
from copy import deepcopy

import micall.utils.referencefull_contig_stitcher_events as full_events


T = TypeVar('T')


@dataclass(frozen=True)
class GenericStitcherContext(Generic[T]):
    uniq_dict: Dict[object, Dict[object, int]]
    events: List[T]

    def register(self, key: object, value: object) -> int:
        if value not in self.uniq_dict:
            self.uniq_dict[value] = {}

        existing = self.uniq_dict[value]
        if key not in existing:
            existing[key] = len(existing) + 1

        return existing[key]

    def emit(self, event: T) -> None:
        self.events.append(event)

    @staticmethod
    def make() -> 'GenericStitcherContext':
        return GenericStitcherContext(events=[], uniq_dict={})

    @staticmethod
    @contextmanager
    def fresh():
        ctx = GenericStitcherContext.make()
        token = context.set(ctx)
        try:
            yield ctx
        finally:
            context.reset(token)

    @staticmethod
    @contextmanager
    def stage():
        ctx = deepcopy(context.get())
        token = context.set(ctx)
        try:
            yield ctx
        finally:
            context.reset(token)


ReferencefullStitcherContext = GenericStitcherContext[full_events.EventType]
ReferencelessStitcherContext = GenericStitcherContext[object]

context: ContextVar[GenericStitcherContext] = ContextVar("GenericStitcherContext")
