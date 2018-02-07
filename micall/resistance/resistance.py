#! /usr/bin/env python3.4
import os
from argparse import ArgumentParser, FileType
from collections import namedtuple, defaultdict
from csv import DictReader, DictWriter
from itertools import groupby, chain, zip_longest
from operator import itemgetter, attrgetter

import yaml
from pyvdrm.vcf import Mutation

from micall.resistance.asi_algorithm import AsiAlgorithm, ResistanceLevels
from micall.core.aln2counts import AMINO_ALPHABET
from micall.utils.sample_sheet_parser import sample_sheet_parser

MIN_FRACTION = 0.05  # prevalence of mutations to report
MIN_COVERAGE = 100
REPORTED_REGIONS = {'PR', 'RT', 'IN', 'NS3', 'NS5a', 'NS5b'}
HIV_RULES_PATH = os.path.join(os.path.dirname(__file__), 'HIVDB_8.3.xml')
HCV_RULES_PATH = os.path.join(os.path.dirname(__file__), 'hcv_rules.yaml')

AminoList = namedtuple('AminoList', 'region aminos genotype')

SampleGroup = namedtuple('SampleGroup', 'enum names')


class LowCoverageError(Exception):
    pass


def parse_args():
    parser = ArgumentParser(
        description='Make resistance calls and list mutations from amino counts.')
    parser.add_argument('aminos_csv', type=FileType(), help='amino counts')
    parser.add_argument('midi_aminos_csv',
                        type=FileType(),
                        help='amino counts for HCV MIDI region or the same as aminos_csv')
    parser.add_argument('resistance_csv',
                        type=FileType('w'),
                        help='resistance calls')
    parser.add_argument('mutations_csv',
                        type=FileType('w'),
                        help='Relevant mutations present above a threshold')
    parser.add_argument('resistance_fail_csv',
                        type=FileType('w'),
                        help='regions that failed to make a resistance call')
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


def create_fail_writer(fail_csv):
    writer = DictWriter(fail_csv,
                        ['seed', 'region', 'reason'],
                        lineterminator=os.linesep)
    writer.writeheader()
    return writer


def check_coverage(region, rows, start_pos=1, end_pos=None):
    if end_pos is None:
        if region.endswith('-NS3'):
            end_pos = 181
        elif region.endswith('-NS5a'):
            end_pos = 101
        elif region.endswith('-NS5b'):
            end_pos = 228
        else:
            return
    start_coverage = 0
    end_coverage = 0
    total_coverage = 0
    for row in rows:
        pos = int(row['refseq.aa.pos'])
        if start_pos <= pos <= end_pos:
            coverage = int(row['coverage'])
            total_coverage += coverage
            if pos == start_pos:
                start_coverage = coverage
            if pos == end_pos:
                end_coverage = coverage
    average_coverage = total_coverage / (end_pos - start_pos + 1)
    if average_coverage < MIN_COVERAGE:
        raise LowCoverageError('low average coverage')
    if (end_coverage < MIN_COVERAGE or
            (start_pos == 1 and start_coverage < MIN_COVERAGE)):
        raise LowCoverageError('not enough high-coverage amino acids')


def combine_aminos(amino_csv, midi_amino_csv, fail_writer):
    midi_rows = defaultdict(list)  # {seed: [row]}
    if midi_amino_csv.name != amino_csv.name:
        for (seed, region), rows in groupby(DictReader(midi_amino_csv),
                                            itemgetter('seed', 'region')):
            if not region.endswith('-NS5b'):
                continue
            rows = list(rows)
            try:
                check_coverage(region, rows, start_pos=231, end_pos=561)
            except LowCoverageError as ex:
                fail_writer.writerow(dict(seed=seed,
                                          region=region,
                                          reason='MIDI: ' + ex.args[0]))
                continue
            midi_rows[seed] = [row
                               for row in rows
                               if 226 < int(row['refseq.aa.pos'])]
    for (seed, region), rows in groupby(DictReader(amino_csv),
                                        itemgetter('seed', 'region')):
        rows = list(rows)
        try:
            check_coverage(region, rows)
        except LowCoverageError as ex:
            fail_writer.writerow(dict(seed=seed,
                                      region=region,
                                      reason=ex.args[0]))
            rows = []
        if region.endswith('-NS5b'):
            region_midi_rows = midi_rows[seed]
            rows = combine_midi_rows(rows, region_midi_rows)
        yield from rows


def combine_midi_rows(main_rows, midi_rows):
    main_row_map = {int(row['refseq.aa.pos']): row
                    for row in main_rows}
    midi_row_map = {int(row['refseq.aa.pos']): row
                    for row in midi_rows}
    positions = sorted({pos
                        for pos in chain(main_row_map.keys(),
                                         midi_row_map.keys())})
    for pos in positions:
        main_row = main_row_map.get(pos)
        midi_row = midi_row_map.get(pos)
        if midi_row is None:
            if pos <= 336:
                yield main_row
        elif main_row is None:
            yield midi_row
        elif (pos <= 336 and
              int(main_row['coverage']) > int(midi_row['coverage'])):
            yield main_row
        else:
            yield midi_row
            

def read_aminos(amino_rows, min_fraction, reported_regions=None, min_coverage=0):
    coverage_columns = list(AMINO_ALPHABET) + ['del']
    report_names = coverage_columns[:]
    report_names[-1] = 'd'
    for (region, seed), rows in groupby(amino_rows,
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
        last_covered_pos = None
        for row in rows:
            counts = list(map(int, (row[f] for f in coverage_columns)))
            coverage = int(row['coverage'])
            pos = int(row['refseq.aa.pos'])
            while pos > len(aminos) + 1:
                aminos.append({})
            if max_pos and pos <= max_pos:
                total_coverage += coverage
            if pos == max_pos:
                last_coverage = coverage
            if coverage == 0 or coverage < min_coverage:
                pos_aminos = {}
            else:
                last_covered_pos = pos
                min_count = max(1, coverage * min_fraction)  # needs at least 1
                pos_aminos = {report_names[i]: count/coverage
                              for i, count in enumerate(counts)
                              if count >= min_count and report_names[i] != '*'}
                ins_count = int(row['ins'])
                if ins_count >= min_count:
                    pos_aminos['i'] = ins_count / coverage
            aminos.append(pos_aminos)
        if (region.endswith('-NS5b') and
                last_covered_pos is not None and
                last_covered_pos < 400):
            # Override last_coverage check when MIDI is missing.
            last_coverage = min_coverage
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
            if region == 'INT':
                region = 'IN'
            result = interpret(asi, amino_seq, region)
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
                    prevalence = pos_aminos.get(amino, 0)
                    mutations_writer.writerow(dict(drug_class=drug_class,
                                                   mutation=mutation,
                                                   prevalence=prevalence,
                                                   genotype=genotype))


def interpret(asi, amino_seq, region):
    ref_seq = asi.stds[region]
    # TODO: Make this more general instead of only applying to missing MIDI.
    is_missing_midi = region.endswith('-NS5b') and len(amino_seq) < len(ref_seq)

    if is_missing_midi:
        amino_seq += [{}] * (len(ref_seq) - len(amino_seq))
    result = asi.interpret(amino_seq, region)

    if not is_missing_midi:
        for drug_result in result.drugs:
            if drug_result.level == ResistanceLevels.FAIL.level:
                break
        else:
            # No missing coverage, we're done.
            return result

    # At least some data found, assume wild type in gaps.
    new_amino_seq = []
    for wild_type, old_aminos in zip_longest(ref_seq, amino_seq):
        if old_aminos:
            new_amino_seq.append(old_aminos)
        else:
            new_amino_seq.append({wild_type: 1.0})
    amino_seq = new_amino_seq

    new_result = asi.interpret(amino_seq, region)
    new_drug_results = {drug_result.code: drug_result
                        for drug_result in new_result.drugs}
    if is_missing_midi:
        allowed_levels = (ResistanceLevels.FAIL.level,
                          ResistanceLevels.RESISTANCE_LIKELY.level,
                          ResistanceLevels.NOT_INDICATED.level)
    else:
        allowed_levels = (ResistanceLevels.FAIL.level,
                          ResistanceLevels.UNKNOWN_MUTATIONS.level,
                          ResistanceLevels.RESISTANCE_POSSIBLE.level,
                          ResistanceLevels.RESISTANCE_LIKELY.level,
                          ResistanceLevels.NOT_INDICATED.level)

    for drug_result in result.drugs:
        if drug_result.level == ResistanceLevels.FAIL.level or is_missing_midi:
            new_drug_result = new_drug_results[drug_result.code]
            if new_drug_result.level in allowed_levels:
                drug_result.level = new_drug_result.level
                drug_result.level_name = new_drug_result.level_name
                drug_result.score = new_drug_result.score
            else:
                drug_result.level = ResistanceLevels.FAIL.level
                drug_result.level_name = ResistanceLevels.FAIL.name
                drug_result.score = 0.0
    result.mutations = new_result.mutations
    return result


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


def find_groups(file_names, sample_sheet_path, included_projects=None):
    with open(sample_sheet_path) as sample_sheet_file:
        run_info = sample_sheet_parser(sample_sheet_file)

    midi_files = {row['sample']: row['filename']
                  for row in run_info['DataSplit']
                  if row['project'] == 'MidHCV'}
    wide_names = {row['filename']: row['sample']
                  for row in run_info['DataSplit']
                  if (row['project'] != 'MidHCV' and
                      (included_projects is None or
                       row['project'] in included_projects))}
    for file_name in file_names:
        trimmed_file_name = '_'.join(file_name.split('_')[:2])
        sample_name = wide_names.get(trimmed_file_name)
        if sample_name is None:
            # Project was not included.
            continue
        midi_file = midi_files.get(sample_name + 'MIDI')
        yield SampleGroup(sample_name, (file_name, midi_file))


def report_resistance(amino_csv,
                      midi_amino_csv,
                      resistance_csv,
                      mutations_csv,
                      fail_csv,
                      region_choices=None):
    if region_choices is None:
        selected_regions = REPORTED_REGIONS
    else:
        selected_regions = select_reported_regions(region_choices, REPORTED_REGIONS)
    fail_writer = create_fail_writer(fail_csv)
    amino_rows = combine_aminos(amino_csv, midi_amino_csv, fail_writer)
    aminos = read_aminos(amino_rows, MIN_FRACTION, selected_regions, MIN_COVERAGE)
    algorithms = load_asi()
    filtered_aminos = filter_aminos(aminos, algorithms)
    write_resistance(filtered_aminos, resistance_csv, mutations_csv, algorithms)


def main():
    args = parse_args()
    report_resistance(args.aminos_csv,
                      args.midi_aminos_csv,
                      args.resistance_csv,
                      args.mutations_csv,
                      args.resistance_fail_csv)


if __name__ == '__main__':
    main()
