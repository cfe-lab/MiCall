
import logging
from typing import List, Tuple, Iterable, Callable


class InMemoryLogHandler(logging.Handler):
    def __init__(self, name: str, *args, **kwargs) -> None:
        super().__init__(*args, **kwargs)
        self.name: str = name
        self.logs: List[logging.LogRecord] = []
        self.callbacks = []


    def emit(self, record: logging.LogRecord):
        self.logs.append(record)
        for callback in self.callbacks:
            callback(record)


    def addCallback(self, callback: Callable[[logging.LogRecord], object]):
        self.callbacks.append(callback)


def add_structured_handler(logger: logging.Logger):
    memory_handler = InMemoryLogHandler(logger.name)
    logger.addHandler(memory_handler)
    return memory_handler
