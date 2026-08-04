from pathlib import Path
import re

from micall.utils.dir_path import DirPath


def find_file(directory: DirPath, pattern: str) -> Path:
    for subdir in directory.iterdir():
        name = subdir.name
        if re.findall(pattern, name):
            return directory / name
    raise ValueError(f"Cannot find file {pattern!r}"
                     f" in directory {str(directory)!r}.")
