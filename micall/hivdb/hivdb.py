#! /usr/bin/env python3.4
import os
from argparse import ArgumentParser, FileType
from csv import DictReader, DictWriter
from itertools import groupby
from operator import itemgetter

from micall.hivdb.asi_algorithm import AsiAlgorithm
from micall.core.aln2counts import AMINO_ALPHABET

MIN_FRACTION = 0.05  # prevalence of mutations to report
MIN_COVERAGE_SCORE = 4
REPORTED_REGIONS = {'PR': 'PR', 'RT': 'RT', 'INT': 'IN'}
RULES_PATH = os.path.join(os.path.dirname(__file__), 'HIVDB_8.3.xml')


def parse_args():
    parser = ArgumentParser(
        description='Make resistance calls and list mutations from amino counts.')
    parser.add_argument('aminos_csv', type=FileType(), help='amino counts')
    parser.add_argument('coverage_scores_csv', type=FileType())
    parser.add_argument('resistance_csv',
                        type=FileType('w'),
                        help='resistance calls')
    parser.add_argument('mutations_csv',
                        type=FileType('w'),
                        help='Relevant mutations present above a threshold')
    return parser.parse_args()


def select_reported_regions(choices, reported_regions):
    split_choices = set()
    for choice in choices:
        split_choices.update(choice.split('_'))
    return {region: translation
            for region, translation in reported_regions.items()
            if region in split_choices}

def find_good_regions(original_regions, coverage_scores_csv):
    good_regions = {region: [name, False] for region, name in original_regions.items()}
    for row in DictReader(coverage_scores_csv):
        score = int(row['on.score'])
        region_code = row['region']
        region = good_regions.get(region_code)
        if region is not None:
            region[1] |= score >= MIN_COVERAGE_SCORE
    return good_regions


def read_aminos(amino_csv, min_fraction, reported_regions=None):
    coverage_columns = list(AMINO_ALPHABET) + ['del']
    report_names = coverage_columns[:]
    report_names[-1] = 'd'
    missing_regions = set()
    if reported_regions:
        missing_regions.update(reported_regions.keys())
    for region, rows in groupby(DictReader(amino_csv),
                                itemgetter('region')):
        if reported_regions is None:
            translated_region = region
        else:
            missing_regions.discard(region)
            translated_region, is_reported = reported_regions.get(region,
                                                                  (None, None))
            if translated_region is None:
                continue
            if not is_reported:
                yield translated_region, None
                continue
        aminos = []
        for row in rows:
            counts = list(map(int, (row[f] for f in coverage_columns)))
            coverage = int(row['coverage'])
            min_count = max(1, coverage * min_fraction)  # needs at least 1
            pos_aminos = {report_names[i]: count/coverage
                          for i, count in enumerate(counts)
                          if count >= min_count and report_names[i] != '*'}
            ins_count = int(row['ins'])
            if ins_count >= min_count:
                pos_aminos['i'] = ins_count / coverage
            aminos.append(pos_aminos)
        yield translated_region, aminos
    for region in missing_regions:
        if reported_regions is None:
            translated_region = region
        else:
            translated_region, _ = reported_regions.get(region, (None, None))
        yield translated_region, None


def write_insufficient_data(resistance_writer, region, asi):
    drug_classes = asi.gene_def[region]
    for drug_class in drug_classes:
        for drug_code in asi.drug_class[drug_class]:
            drug_name = asi.drugs[drug_code][0]
            resistance_writer.writerow(dict(region=region,
                                            drug_class=drug_class,
                                            drug=drug_code,
                                            drug_name=drug_name,
                                            level_name='Insufficient data available',
                                            level=0,
                                            score=0.0))


def write_resistance(aminos, resistance_csv, mutations_csv):
    resistance_writer = DictWriter(
        resistance_csv,
        ['region', 'drug_class', 'drug', 'drug_name', 'level', 'level_name', 'score'],
        lineterminator=os.linesep)
    resistance_writer.writeheader()
    mutations_writer = DictWriter(mutations_csv,
                                  ['drug_class', 'mutation', 'prevalence'],
                                  lineterminator=os.linesep)
    mutations_writer.writeheader()
    asi = AsiAlgorithm(RULES_PATH)
    for region, amino_seq in aminos:
        if amino_seq is None:
            write_insufficient_data(resistance_writer, region, asi)
            continue
        result = asi.interpret(amino_seq, region)
        for drug_result in result.drugs:
            resistance_writer.writerow(dict(region=region,
                                            drug_class=drug_result.drug_class,
                                            drug=drug_result.code,
                                            drug_name=drug_result.name,
                                            level_name=drug_result.level_name,
                                            level=drug_result.level,
                                            score=drug_result.score))
        for drug_class, class_mutations in result.mutations.items():
            for mutation in class_mutations:
                amino = mutation[-1]
                pos = int(mutation[1:-1])
                pos_aminos = amino_seq[pos-1]
                prevalence = pos_aminos[amino]
                mutations_writer.writerow(dict(drug_class=drug_class,
                                               mutation=mutation,
                                               prevalence=prevalence))


def hivdb(amino_csv,
          coverage_scores_csv,
          resistance_csv,
          mutations_csv,
          region_choices=None):
    if region_choices is None:
        selected_regions = REPORTED_REGIONS
    else:
        selected_regions = select_reported_regions(region_choices, REPORTED_REGIONS)
    good_regions = find_good_regions(selected_regions, coverage_scores_csv)
    aminos = read_aminos(amino_csv, MIN_FRACTION, good_regions)
    write_resistance(aminos, resistance_csv, mutations_csv)


def main():
    args = parse_args()
    hivdb(args.aminos_csv,
          args.coverage_scores_csv,
          args.resistance_csv,
          args.mutations_csv)

if __name__ == '__main__':
    main()
