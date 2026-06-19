#! /usr/bin/env python3
import csv
import json
import os
from argparse import ArgumentParser, FileType, ArgumentDefaultsHelpFormatter

import subprocess

import typing

from micall.utils.externals import ProjectsFile


def parse_args():
    parser = ArgumentParser(description='Sort SAM file before viewing.',
                            formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('sam', help='SAM file to sort, or prelim.csv')
    with ProjectsFile().path() as projects_file_path:
        parser.add_argument('--projects',
                            type=FileType(),
                            help='JSON file with project definitions',
                            default=str(projects_file_path))
    return parser.parse_args()


def convert_from_csv(csv_name: str,
                     sam_name: str,
                     ref_name: str,
                     projects_file: typing.TextIO):
    with open(csv_name) as csv_file, \
            open(sam_name, 'w') as sam_file, \
            open(ref_name, 'w') as ref_file:
        writer = csv.writer(sam_file, delimiter='\t', lineterminator=os.linesep)
        writer.writerow(['@HD', 'VN:1.0', 'SO:unsorted'])

        reader = csv.reader(csv_file)
        header = next(reader)
        region_index = 2
        assert header[region_index] == 'rname'

        regions = {row[region_index] for row in reader}

        projects = json.load(projects_file)
        for region, info in projects['regions'].items():
            if region in regions:
                reference = ''.join(info['reference'])
                reference_length = len(reference)
                row = ['@SQ', f'SN:{region}', f'LN:{reference_length}']
                writer.writerow(row)
                print(f'>{region}', file=ref_file)
                for line in info['reference']:
                    print(line, file=ref_file)

        csv_file.seek(0)
        reader = csv.reader(csv_file)
        next(reader)
        writer.writerows(reader)


def main():
    args = parse_args()
    sam_name = args.sam
    sam_root, sam_ext = os.path.splitext(sam_name)
    if sam_ext == '.csv':
        csv_name = sam_name
        sam_name = sam_root + '.sam'
        ref_name = sam_root + '_ref.fasta'
        convert_from_csv(csv_name, sam_name, ref_name, args.projects)
    # samtools view -Sb example.sam -o example.bam
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


main()
