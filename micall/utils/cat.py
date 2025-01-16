#! /usr/bin/env python

"""
This script concatenates files and writes their contents to an output file.

It mimics the behavior of the Unix `cat` command, with the additional
feature of specifying an output file where the concatenated content is
saved. It no longer supports reading from the standard input.

The script uses efficient techniques to handle potentially large files
by leveraging the `shutil.copyfileobj()` method to stream file
contents.
"""

import argparse
import shutil
import sys
from pathlib import Path
from typing import Iterable, Sequence


def cat(inputs: Iterable[Path], output: Path) -> int:
    """
    Concatenates the contents of input inputs and writes the result to
    an output file.
    """

    with open(output, 'w') as out:
        try:
            for file in inputs:
                try:
                    with open(file) as f:
                        shutil.copyfileobj(f, out)
                except FileNotFoundError:
                    print(f"The file {str(file)!r} was not found."
                          " Please check the file path and try again.",
                          file=sys.stderr)
                    return 1
                except IOError as e:
                    print(f"IO Error while processing the file {str(file)!r}: {e}",
                          file=sys.stderr)
                    return 1

            return 0
        except BaseException as e:
            print(f"An unexpected error occurred: {e}", file=sys.stderr)
            return 1


def main(argv: Sequence[str]) -> int:
    parser = argparse.ArgumentParser(
        description="Concatenate inputs and write to an output file.")

    parser.add_argument(
        'inputs',
        metavar='FILE',
        type=Path,
        nargs='+',
        help='Input inputs to concatenate.'
    )

    parser.add_argument(
        'output',
        metavar='OUTPUT_FILE',
        type=Path,
        help='Output file where the concatenated result will be saved.'
    )

    args = parser.parse_args(argv)
    return cat(args.inputs, args.output)


def entry() -> None:
    sys.exit(main(sys.argv[1:]))


if __name__ == '__main__':
    entry()
