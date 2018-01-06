#! /usr/bin/env python3.4
import json
import os
from argparse import ArgumentParser, FileType
from collections import namedtuple
from csv import DictReader, DictWriter
from itertools import groupby
from operator import itemgetter

from micall.hivdb.asi_algorithm import AsiAlgorithm
from micall.core.aln2counts import AMINO_ALPHABET

MIN_FRACTION = 0.05  # prevalence of mutations to report
MIN_COVERAGE_SCORE = 4
REPORTED_REGIONS = {'PR', 'RT', 'IN', 'NS3', 'NS5a', 'NS5b'}
HIV_RULES_PATH = os.path.join(os.path.dirname(__file__), 'HIVDB_8.3.xml')
HCV_RULES_PATH = os.path.join(os.path.dirname(__file__), 'hcv_rules.json')

AminoList = namedtuple('AminoList', 'region aminos seed')


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
    return split_choices & reported_regions


def get_reported_region(reference):
    if reference == 'INT':
        return 'IN'
    for region in ('NS3', 'NS5a', 'NS5b'):
        if reference.startswith('HCV') and reference.endswith(region):
            return region
    return reference


def get_genotype(seed):
    if seed is None:
        return None
    parts = seed.split('-')
    virus = parts[0]
    if virus != 'HCV':
        return None
    full_genotype = parts[1]
    return full_genotype[0]


def find_good_regions(original_regions, coverage_scores_csv):
    good_regions = {}
    for row in DictReader(coverage_scores_csv):
        region_code = row['region']
        entry = good_regions.get(region_code)
        if entry is None:
            reported_region = get_reported_region(region_code)
            if reported_region not in original_regions:
                continue
            entry = [reported_region, False]
            good_regions[region_code] = entry
        score = int(row['on.score'])
        entry[1] |= score >= MIN_COVERAGE_SCORE
    return good_regions


def read_aminos(amino_csv, min_fraction, reported_regions=None):
    coverage_columns = list(AMINO_ALPHABET) + ['del']
    report_names = coverage_columns[:]
    report_names[-1] = 'd'
    missing_regions = set()
    if reported_regions:
        missing_regions.update(reported_regions.keys())
    for (region, seed), rows in groupby(DictReader(amino_csv),
                                itemgetter('region', 'seed')):
        if reported_regions is not None:
            missing_regions.discard(region)
            translated_region, is_reported = reported_regions.get(region,
                                                                  (None, None))
            if translated_region is None:
                continue
            if not is_reported:
                yield AminoList(region, None, None)
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
        yield AminoList(region, aminos, seed)
    for region in missing_regions:
        yield AminoList(region, None, None)


def write_insufficient_data(resistance_writer, region, asi):
    reported_region = get_reported_region(region)
    drug_classes = asi.gene_def[region]
    for drug_class in drug_classes:
        for drug_code in asi.drug_class[drug_class]:
            drug_name = asi.drugs[drug_code][0]
            resistance_writer.writerow(dict(region=reported_region,
                                            drug_class=drug_class,
                                            drug=drug_code,
                                            drug_name=drug_name,
                                            level_name='Insufficient data available',
                                            level=0,
                                            score=0.0))


def write_resistance(aminos, resistance_csv, mutations_csv):
    """ Calculate resistance scores and write them to files.

    :param list[AminoList] aminos: region is the coordinate
        reference name that this gene region was mapped to, and prevalance is a
        float between 0.0 and 1.0
    :param resistance_csv: open file to write resistance calls to, grouped by
        genotype, region, drug_class
    :param mutations_csv: open file to write mutations to, grouped by genotype,
        drug_class
    """
    resistance_writer = DictWriter(
        resistance_csv,
        ['region',
         'drug_class',
         'drug',
         'drug_name',
         'level',
         'level_name',
         'score',
         'genotype'],
        lineterminator=os.linesep)
    resistance_writer.writeheader()
    mutations_writer = DictWriter(mutations_csv,
                                  ['drug_class',
                                   'mutation',
                                   'prevalence',
                                   'genotype'],
                                  lineterminator=os.linesep)
    mutations_writer.writeheader()
    algorithms = load_asi()
    for region, amino_seq, seed in aminos:
        asi = algorithms.get(region)
        if asi is None:
            continue
        reported_region = get_reported_region(region)
        genotype = get_genotype(seed)
        if amino_seq is None:
            write_insufficient_data(resistance_writer, region, asi)
            continue
        result = asi.interpret(amino_seq, region)
        for drug_result in result.drugs:
            resistance_writer.writerow(dict(region=reported_region,
                                            drug_class=drug_result.drug_class,
                                            drug=drug_result.code,
                                            drug_name=drug_result.name,
                                            level_name=drug_result.level_name,
                                            level=drug_result.level,
                                            score=drug_result.score,
                                            genotype=genotype))
        for drug_class, class_mutations in result.mutations.items():
            for mutation in class_mutations:
                amino = mutation[-1]
                pos = int(mutation[1:-1])
                pos_aminos = amino_seq[pos-1]
                prevalence = pos_aminos[amino]
                mutations_writer.writerow(dict(drug_class=drug_class,
                                               mutation=mutation,
                                               prevalence=prevalence,
                                               genotype=genotype))


def load_asi():
    asi = AsiAlgorithm(HIV_RULES_PATH)
    algorithms = {'PR': asi,
                  'RT': asi,
                  'INT': asi}
    with open(HCV_RULES_PATH) as f:
        hcv_rules = json.load(f)
    references = {genotype['reference']
                  for rule in hcv_rules
                  for genotype in rule['genotypes']}
    for reference in references:
        algorithms[reference] = AsiAlgorithm(rules_json=hcv_rules,
                                             reference=reference)

    return algorithms


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
