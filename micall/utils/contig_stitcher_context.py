from typing import Iterable, Optional, Tuple, List, Dict, Union, Literal, TypeVar, Callable, Set
from contextvars import ContextVar
from contextlib import contextmanager

import micall.utils.contig_stitcher_events as events


class StitcherContext:
    def __init__(self) -> None:
        self.name_generator_state: int = 0
        self.nameset: Set[str] = set()
        self.events: List[events.EventType] = []

    def generate_new_name(self) -> str:
        while True:
            self.name_generator_state += 1
            name = f"c{self.name_generator_state}"
            if name not in self.nameset:
                self.nameset.add(name)
                return name

    def emit(self, event: events.EventType) -> None:
        self.events.append(event)


    @staticmethod
    @contextmanager
    def fresh():
        ctx = StitcherContext()
        token = context.set(ctx)
        try:
            yield ctx
        finally:
            context.reset(token)


context: ContextVar[StitcherContext] = ContextVar("StitcherContext")
