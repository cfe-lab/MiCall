"""
This module provides methods for creating directories
that only materialize if they are fully initialized.
"""

import random
import shutil
import string
from collections.abc import Iterator
from contextlib import contextmanager

from micall.utils.dir_path import DirPath


def random_string(length: int) -> str:
    chars = string.ascii_letters + string.digits
    return ''.join(random.choices(chars, k=length))


@contextmanager
def new_atomic_directory(path: DirPath) -> Iterator[DirPath]:
    random_name_part = random_string(9)
    temporary_name = f".newatomicdir-{random_name_part}.{path.name}~"
    temporary_path = path.parent / temporary_name
    temporary_path.mkdir(exist_ok=False, parents=True)

    try:
        yield temporary_path
    except:
        shutil.rmtree(temporary_path)
        raise

    if path.exists():
        shutil.rmtree(path)

    temporary_path.rename(path)
