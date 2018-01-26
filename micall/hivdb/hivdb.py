#! /usr/bin/env python3.4
import os
from argparse import ArgumentParser, FileType
from collections import namedtuple
from csv import DictReader, DictWriter
from itertools import groupby
from operator import itemgetter, attrgetter

import yaml
from pyvdrm.vcf import Mutation

from micall.hivdb.asi_algorithm import AsiAlgorithm, ResistanceLevels
from micall.core.aln2counts import AMINO_ALPHABET

MIN_FRACTION = 0.05  # prevalence of mutations to report
MIN_COVERAGE = 100
REPORTED_REGIONS = {'PR', 'RT', 'IN', 'NS3', 'NS5a', 'NS5b'}
HIV_RULES_PATH = os.path.join(os.path.dirname(__file__), 'HIVDB_8.3.xml')
HCV_RULES_PATH = os.path.join(os.path.dirname(__file__), 'hcv_rules.yaml')

AminoList = namedtuple('AminoList', 'region aminos genotype')


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
    full_genotype = parts[1].upper()
    if full_genotype.startswith('1'):
        if full_genotype == '1B':
            return full_genotype
        return '1A'
    if full_genotype == '6E':
        return full_genotype
    return full_genotype[0]


def read_aminos(amino_csv, min_fraction, reported_regions=None, min_coverage=0):
    coverage_columns = list(AMINO_ALPHABET) + ['del']
    report_names = coverage_columns[:]
    report_names[-1] = 'd'
    for (region, seed), rows in groupby(DictReader(amino_csv),
                                        itemgetter('region', 'seed')):
        genotype = get_genotype(seed)
        if reported_regions is not None:
            translated_region = get_reported_region(region)
            if translated_region not in reported_regions:
                continue
        aminos = []
        total_coverage = 0
        if region.endswith('-NS3'):
            max_pos = 181
        elif region.endswith('-NS5a'):
            max_pos = 101
        elif region.endswith('-NS5b'):
            max_pos = 561
        else:
            max_pos = None
        # TODO: Remove last_coverage check after validating against Ruby version.
        last_coverage = 0
        for row in rows:
            counts = list(map(int, (row[f] for f in coverage_columns)))
            coverage = int(row['coverage'])
            pos = int(row['refseq.aa.pos'])
            if max_pos and pos <= max_pos:
                total_coverage += coverage
            if pos == max_pos:
                last_coverage = coverage
            if coverage == 0 or coverage < min_coverage:
                pos_aminos = {}
            else:
                min_count = max(1, coverage * min_fraction)  # needs at least 1
                pos_aminos = {report_names[i]: count/coverage
                              for i, count in enumerate(counts)
                              if count >= min_count and report_names[i] != '*'}
                ins_count = int(row['ins'])
                if ins_count >= min_count:
                    pos_aminos['i'] = ins_count / coverage
            aminos.append(pos_aminos)
        if max_pos is None or min(last_coverage,
                                  total_coverage // max_pos) >= min_coverage:
            # Need original region to look up wild type.
            yield AminoList(region, aminos, genotype)


def get_algorithm_regions(algorithm):
    return ('INT' if region == 'IN' else region
            for region in algorithm.gene_def)


def create_empty_aminos(region, genotype, algorithms):
    algorithm = algorithms[genotype]
    std_name = 'IN' if region == 'INT' else region
    std_length = len(algorithm.stds[std_name])
    return AminoList(region,
                     [{}] * std_length,
                     genotype)


def filter_aminos(all_aminos, algorithms):
    all_aminos = list(all_aminos)
    good_aminos = [amino_list
                   for amino_list in all_aminos
                   if any(amino_list.aminos)]
    good_genotypes = {amino_list.genotype for amino_list in good_aminos}
    good_regions = {(amino_list.genotype, amino_list.region)
                    for amino_list in good_aminos}
    expected_regions = {(genotype, region)
                        for genotype in good_genotypes
                        for region in get_algorithm_regions(algorithms[genotype])}
    missing_regions = sorted(expected_regions - good_regions)
    good_aminos += [create_empty_aminos(region, genotype, algorithms)
                    for genotype, region in missing_regions]
    return good_aminos


def write_resistance(aminos, resistance_csv, mutations_csv, algorithms=None):
    """ Calculate resistance scores and write them to files.

    :param list[AminoList] aminos: region is the coordinate
        reference name that this gene region was mapped to, and prevalance is a
        float between 0.0 and 1.0
    :param resistance_csv: open file to write resistance calls to, grouped by
        genotype, region, drug_class
    :param mutations_csv: open file to write mutations to, grouped by genotype,
        drug_class
    :param dict algorithms: {region: AsiAlgorithm}
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
    if algorithms is None:
        algorithms = load_asi()
    for genotype, genotype_aminos in groupby(aminos, attrgetter('genotype')):
        region_results = []
        for region, amino_seq, _ in genotype_aminos:
            asi = algorithms.get(genotype)
            if asi is None:
                continue
            result = asi.interpret(amino_seq, region)
            region_results.append((region, amino_seq, result))
        if all(drug_result.level == ResistanceLevels.FAIL.level
               for region, amino_seq, result in region_results
               for drug_result in result.drugs):
            continue
        for region, amino_seq, result in region_results:
            reported_region = get_reported_region(region)
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
                mutations = [Mutation(m) for m in class_mutations]
                mutations.sort()
                for mutation in mutations:
                    amino = mutation.variant
                    pos = mutation.pos
                    pos_aminos = amino_seq[pos-1]
                    prevalence = pos_aminos[amino]
                    mutations_writer.writerow(dict(drug_class=drug_class,
                                                   mutation=mutation,
                                                   prevalence=prevalence,
                                                   genotype=genotype))


def load_asi():
    asi = AsiAlgorithm(HIV_RULES_PATH)
    algorithms = {None: asi}
    with open(HCV_RULES_PATH) as f:
        hcv_rules = yaml.safe_load(f)
    genotypes = {genotype['genotype']
                 for rule in hcv_rules
                 for genotype in rule['genotypes']}
    for genotype in genotypes:
        algorithms[genotype] = AsiAlgorithm(rules_yaml=hcv_rules,
                                            genotype=genotype)

    return algorithms


def hivdb(amino_csv,
          resistance_csv,
          mutations_csv,
          region_choices=None):
    if region_choices is None:
        selected_regions = REPORTED_REGIONS
    else:
        selected_regions = select_reported_regions(region_choices, REPORTED_REGIONS)
    aminos = read_aminos(amino_csv, MIN_FRACTION, selected_regions, MIN_COVERAGE)
    algorithms = load_asi()
    filtered_aminos = filter_aminos(aminos, algorithms)
    write_resistance(filtered_aminos, resistance_csv, mutations_csv, algorithms)


def main():
    args = parse_args()
    hivdb(args.aminos_csv, args.resistance_csv, args.mutations_csv)


if __name__ == '__main__':
    main()
