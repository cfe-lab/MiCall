#! /usr/bin/env python3

import argparse
import sys
from typing import Sequence
import os
from pathlib import Path
import re


def dir_path(string: str) -> Path:
    if os.path.isdir(string):
        return Path(string)
    else:
        raise ValueError("Path %r is not a directory.", string)


def parse_arguments(argv: Sequence[str]) -> argparse.Namespace:
    """
    Parse command-line arguments.

    Args:
        argv: A list of command-line arguments.

    Returns:
        A Namespace object containing parsed arguments.
    """

    script_path: Path = Path(__file__).resolve()
    base_dir = script_path.parent.parent.parent

    parser = argparse.ArgumentParser(
        description="List executable Python scripts.")
    parser.add_argument(
        '-d', '--directory', type=dir_path, default=base_dir,
        help='The root directory to search for Python files (default: parent of the current script\'s directory).'
    )

    return parser.parse_args(argv)


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


def main(argv: Sequence[str]) -> int:
    """
    Main function to list the script files.

    Args:
        argv: A list of command-line arguments.

    Returns:
        An exit status code (0 for success).
    """

    args = parse_arguments(argv)
    base_dir: Path = args.directory

    # Iterate over all files in the specified directory
    for root, _, files in os.walk(base_dir):
        for file in files:
            # Process only files with a .py extension
            if not file.endswith('.py'):
                continue

            file_path = os.path.join(root, file)
            try:
                with open(file_path, 'r') as f:
                    content = f.read()
            except (IOError, OSError) as e:
                print(f"Error reading {file_path}: {e}", file=sys.stderr)
                continue

            if is_executable_script(content):
                print(file_path)

    return 0


if __name__ == "__main__":
    exit(main(sys.argv[1:]))
