"""
This module provides methods for creating files atomically.
"""

from contextlib import contextmanager
from pathlib import Path
from typing import TextIO, Iterator
import random
import string


def random_string(length: int) -> str:
    chars = string.ascii_letters + string.digits
    return ''.join(random.choices(chars, k=length))


@contextmanager
def new_atomic_text_file(path: Path) -> Iterator[TextIO]:
    random_name_part = random_string(9)
    temporary_name = f".newatomicfile-{random_name_part}.{path.name}~"
    temporary_path = path.parent / temporary_name
    with temporary_path.open("w+t") as writer:
        try:
            yield writer
        except:
            temporary_path.unlink()
            raise
    temporary_path.rename(path)
