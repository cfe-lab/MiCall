
from typing import MutableMapping
from contextlib import contextmanager
from contextvars import ContextVar


class Registry:
    def __init__(self) -> None:
        self._uniq_dict: MutableMapping[object, MutableMapping[object, int]] \
            = {}

    def add(self, key: object, value: object) -> int:
        if value not in self._uniq_dict:
            self._uniq_dict[value] = {}

        existing = self._uniq_dict[value]
        if key not in existing:
            existing[key] = len(existing) + 1

        return existing[key]


@contextmanager
def fresh():
    ctx = Registry()
    token = context.set(ctx)
    try:
        yield ctx
    finally:
        context.reset(token)


context: ContextVar[Registry] = ContextVar("Registry")
