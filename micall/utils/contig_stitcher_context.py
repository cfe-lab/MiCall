from typing import List, Dict
from contextvars import ContextVar
from contextlib import contextmanager
from dataclasses import dataclass
from copy import deepcopy

import micall.utils.contig_stitcher_events as st_events


@dataclass
class StitcherContext:
    uniq_dict: Dict[object, Dict[object, int]]
    events: List[st_events.EventType]

    def register(self, key: object, value: object) -> int:
        if value not in self.uniq_dict:
            self.uniq_dict[value] = {}

        existing = self.uniq_dict[value]
        if key not in existing:
            existing[key] = len(existing) + 1

        return existing[key]

    def emit(self, event: st_events.EventType) -> None:
        self.events.append(event)

    @staticmethod
    def make() -> 'StitcherContext':
        return StitcherContext(events=[], uniq_dict={})

    @staticmethod
    @contextmanager
    def fresh():
        ctx = StitcherContext.make()
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


context: ContextVar[StitcherContext] = ContextVar("StitcherContext")
