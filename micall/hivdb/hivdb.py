#! /usr/bin/env python3.4
import os
from argparse import ArgumentParser, FileType
from csv import DictReader, DictWriter
from itertools import groupby
from operator import itemgetter

from micall.hivdb.asi_algorithm import AsiAlgorithm
from ..core.aln2counts import AMINO_ALPHABET

MIN_FRACTION = 0.05  # prevalence of mutations to report
REPORTED_REGIONS = ('PR', 'RT', 'INT')
RULES_PATH = os.path.join(os.path.dirname(__file__), 'HIVDB_8.3.xml')


def parse_args():
    parser = ArgumentParser(
        description='Make resistance calls and list mutations from amino counts.')
    parser.add_argument('aminos_csv', type=FileType(), help='amino counts')
    parser.add_argument('resistance_csv',
                        type=FileType('w'),
                        help='resistance calls')
    # parser.add_argument('mutations_csv',
    #                     type=FileType('w'),
    #                     help='Relevant mutations present above a threshold')
    return parser.parse_args()


def read_aminos(amino_csv, min_fraction, reported_regions=None):
    coverage_columns = list(AMINO_ALPHABET) + ['del']
    report_names = coverage_columns[:]
    report_names[-1] = 'd'
    for region, rows in groupby(DictReader(amino_csv),
                                itemgetter('region')):
        if reported_regions is not None and region not in reported_regions:
            continue
        aminos = []
        for row in rows:
            counts = list(map(int, (row[f] for f in coverage_columns)))
            coverage = sum(counts)
            min_count = max(1, coverage * min_fraction)  # needs at least 1
            pos_aminos = [report_names[i]
                          for i, count in enumerate(counts)
                          if count >= min_count and report_names[i] != '*']
            ins_count = int(row['ins'])
            if ins_count >= min_count:
                pos_aminos.append('i')
            aminos.append(pos_aminos)
        yield region, aminos


def write_resistance(resistance_csv, aminos):
    writer = DictWriter(resistance_csv,
                        ['region', 'drug', 'level_name', 'level', 'score'],
                        lineterminator=os.linesep)
    writer.writeheader()
    asi = AsiAlgorithm(RULES_PATH)
    for region, amino_seq in aminos:
        result = asi.interpret(amino_seq, region)
        for drug_result in result.drugs:
            writer.writerow(dict(region=region,
                                 drug=drug_result.code,
                                 level_name=drug_result.level_name,
                                 level=drug_result.level,
                                 score=drug_result.score))


def hivdb(amino_csv, resistance_csv):
    aminos = read_aminos(amino_csv, MIN_FRACTION, REPORTED_REGIONS)
    write_resistance(resistance_csv, aminos)


def main():
    args = parse_args()
    hivdb(args.amino_csv, args.resistance_csv)

if __name__ == '__main__':
    main()
