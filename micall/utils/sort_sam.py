#! /usr/bin/env python3
import os
from argparse import ArgumentParser, FileType

import subprocess


def parse_args():
    parser = ArgumentParser(description='Sort SAM file before viewing.')
    parser.add_argument('sam', help='SAM file to sort')

    return parser.parse_args()


def main():
    args = parse_args()
    # samtools view -Sb example.sam -o example.bam
    sam_name = args.sam
    sam_root, _ = os.path.splitext(sam_name)
    subprocess.check_call(
        ['samtools', 'view', '-Sb', sam_name, '-o', sam_root + '.bam'])
    # samtools sort example.bam -o example.sorted.bam
    subprocess.check_call(
        ['samtools', 'sort', sam_root + '.bam', '-o', sam_root + '.sorted.bam'])
    # samtools view -h -o example.sorted.sam example.sorted.bam
    subprocess.check_call(['samtools',
                           'view',
                           '-h',
                           '-o',
                           sam_root + '.sorted.sam',
                           sam_root + '.sorted.bam'])

if __name__ == '__main__':
    main()
