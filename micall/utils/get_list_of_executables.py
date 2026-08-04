#! /usr/bin/env python3

import argparse
import sys
from typing import Sequence, Iterator
import os
from pathlib import Path
import re


def is_executable_script(content: str) -> bool:
    if content.startswith("#!"):
        return True

    if re.findall(r'__name__.+__main__', content):
        return True

    return False


def iterate_executables() -> Iterator[Path]:
    script_path: Path = Path(__file__).resolve()
    base_dir = script_path.parent.parent.parent

    # Only scan specific source directories, not the entire project
    source_dirs = ["micall"]

    for source_dir in source_dirs:
        source_path = base_dir / source_dir
        if not source_path.exists():
            continue

        # Iterate over all files in the source directory.
        for root, _, files in os.walk(source_path):
            for file in files:
                # Process only files with a .py extension.
                if not file.endswith(".py"):
                    continue

                file_path = Path(root) / file
                content = file_path.read_text()

                if is_executable_script(content):
                    relative = file_path.relative_to(base_dir)
                    yield relative


def main(argv: Sequence[str]) -> int:
    """
    Main function to list the script files.

    Args:
        argv: A list of command-line arguments.

    Returns:
        An exit status code (0 for success).
    """

    parser = argparse.ArgumentParser(description="List executable Python scripts.")
    parser.parse_args(argv)

    for path in iterate_executables():
        print(path)

    return 0


if __name__ == "__main__":
    exit(main(sys.argv[1:]))
