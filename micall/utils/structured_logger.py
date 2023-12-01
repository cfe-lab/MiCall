
import logging
from typing import List, Tuple, Iterable


LoggerName = str
_structured_logs: List[Tuple[LoggerName, logging.LogRecord]] = []


class InMemoryLogHandler(logging.Handler):
    def __init__(self, name: str, *args, **kwargs) -> None:
        super().__init__(*args, **kwargs)
        self.name: str = name

    def emit(self, record: logging.LogRecord):
        _structured_logs.append((self.name, record))


def register_structured_logger(logger: logging.Logger):
    memory_handler = InMemoryLogHandler(logger.name)
    logger.addHandler(memory_handler)


def iterate_messages() -> Iterable[Tuple[LoggerName, logging.LogRecord]]:
    yield from _structured_logs
