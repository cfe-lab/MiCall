
import sys
from importlib.metadata import version

import micall.core.contig_stitcher as stitcher
import release_test_publish as pub

assert stitcher
assert sys
assert pub


def get_version() -> str:
    if __package__ is None:
        return "development"
    else:
        return str(version(__package__))


def cli():
    print(get_version())
    print('Bye!')
