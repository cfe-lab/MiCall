
from typing import MutableMapping
from contextlib import contextmanager
from contextvars import ContextVar
from copy import deepcopy


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
    token = _context.set(ctx)
    try:
        yield ctx
    finally:
        _context.reset(token)


@contextmanager
def stage():
    try:
        existing = get()
    except BaseException:
        with fresh() as ret:
            yield ret
            return

    ctx = deepcopy(existing)
    token = _context.set(ctx)
    try:
        yield ctx
    finally:
        _context.reset(token)


def get() -> Registry:
    return _context.get()


def set(value: Registry) -> None:
    _context.set(value)


_context: ContextVar[Registry] = ContextVar("Registry")
