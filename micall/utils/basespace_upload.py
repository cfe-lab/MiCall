#!/usr/bin/env python3.6

""" Upload samples to Illumina BaseSpace.

This script depends on the BaseSpace command-line interface.
https://help.basespace.illumina.com/articles/descriptive/basespace-cli
"""

import os
from argparse import ArgumentParser
from itertools import groupby
from subprocess import check_call


def parse_args():
    parser = ArgumentParser(description='Upload a set of files to BaseSpace.')
    parser.add_argument('name', help='a name for the new project')
    parser.add_argument('files', nargs='+', help='files to upload')

    return parser.parse_args()


def upload(project_name, filenames):
    for _, group_files in groupby(
            sorted(filenames),
            key=lambda s: os.path.basename(s).split('_')[0]):
        upload_args = ['bs',
                       'upload',
                       'sample',
                       '-p', project_name]
        upload_args.extend(group_files)
        check_call(upload_args)


def main():
    args = parse_args()
    upload(args.name, args.files)
    print('Done.')


if __name__ == '__main__':
    main()
