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

    if re.findall(r'__name__\s*==\s*[\'"]__main__', content):
        return True

    if 'import argparse' in content:
        return True

    if 'from argparse' in content:
        return True

    return False


def iterate_executables() -> Iterator[Path]:
    script_path: Path = Path(__file__).resolve()
    base_dir = script_path.parent.parent.parent

    # Iterate over all files in the base directory.
    for root, _, files in os.walk(base_dir):
        for file in files:

            # Process only files with a .py extension.
            if not file.endswith('.py'):
                continue

            file_path = Path(root) / file
            content = file_path.read_text()

            if is_executable_script(content):
                yield file_path


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
