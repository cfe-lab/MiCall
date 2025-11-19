import os
import typing
from argparse import ArgumentParser, FileType
from collections import namedtuple
from csv import DictReader, DictWriter
from itertools import groupby, chain, zip_longest
from operator import itemgetter, attrgetter

import yaml
from pyvdrm.vcf import Mutation

from micall.core.project_config import ProjectConfig
from micall.data.landmark_reader import LandmarkReader
from micall.resistance.asi_algorithm import AsiAlgorithm, HcvResistanceLevels, HivResistanceLevels
from micall.core.aln2counts import AMINO_ALPHABET

MIN_FRACTION = 0.05  # prevalence of mutations to report
MIN_COVERAGE = 100
REPORTED_REGIONS = {'PR', 'RT', 'IN', 'NS3', 'NS5a', 'NS5b'}
LAST_MAIN_POS = 336  # last position to take from main file if midi is missing
NS3_END_POS = 181  # positions to check for coverage
NS5A_END_POS = 101
NS5B_MAIN_END_POS = 228
MIDI_START_POS = 231
MIDI_END_POS = 561

# Rules configuration - remember to update version numbers in genreport.yaml.
HIVDB_VERSION = '9.8'
HIV_RULES_PATH = os.path.join(os.path.dirname(__file__), f'HIVDB_{HIVDB_VERSION}.xml')
HCV_RULES_PATH = os.path.join(os.path.dirname(__file__), 'hcv_rules.yaml')

NOTHING_MAPPED_MESSAGE = 'nothing mapped'

AminoList = namedtuple('AminoList', 'region aminos genotype seed is_report_needed')
AminoList.__new__.__defaults__ = (False,)  # default for is_report_needed


class LowCoverageError(Exception):
    pass


def parse_args():
    parser = ArgumentParser(
        description='Make resistance calls and list mutations from amino counts.')
    parser.add_argument('main_amino_csv',
                        type=FileType(),
                        help='CSV containing amino frequencies from main sample')
    parser.add_argument('midi_amino_csv',
                        type=FileType(),
                        help='CSV containing amino frequencies from MIDI sample')
    parser.add_argument('main_nuc_csv',
                        type=FileType(),
                        help='CSV containing nucleotide frequencies from main sample')
    parser.add_argument('resistance_csv',
                        type=FileType('w'),
                        help='resistance calls')
    parser.add_argument('mutations_csv',
                        type=FileType('w'),
                        help='Relevant mutations present above a threshold')
    parser.add_argument('nuc_mutations_csv',
                        type=FileType('w'),
                        help='Relevant nucleotide mutations present above a threshold')
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
                        ['seed', 'region', 'coord_region', 'genotype', 'reason'],
                        lineterminator=os.linesep)
    writer.writeheader()
    return writer


def check_coverage(region, rows, start_pos=1, end_pos=None):
    if end_pos is None:
        if region.endswith('-NS3'):
            end_pos = NS3_END_POS
        elif region.endswith('-NS5a'):
            end_pos = NS5A_END_POS
        elif region.endswith('-NS5b'):
            end_pos = NS5B_MAIN_END_POS
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


def combine_aminos(amino_csv, midi_amino_csv, failures: dict):
    """ Combine amino rows from the two amplified regions.

    :param amino_csv: an open file with the amino counts from the whole genome
        region.
    :param midi_amino_csv: an open file with the amino counts from the MIDI
        region.
    :param failures: {(seed, region, is_midi): message} messages for any regions
        that did not have enough coverage. is_midi is True if the region was
        in the midi_amino_csv file.
    """
    is_midi = True
    midi_rows = {}  # {genotype: [row]}
    midi_start = MIDI_START_POS
    midi_end = MIDI_END_POS
    is_midi_separate = midi_amino_csv.name != amino_csv.name
    if is_midi_separate:
        for (seed, region), rows in groupby(DictReader(midi_amino_csv),
                                            itemgetter('seed', 'region')):
            if not region.endswith('-NS5b'):
                continue
            rows = list(rows)
            try:
                check_coverage(region, rows, start_pos=midi_start, end_pos=midi_end)
            except LowCoverageError as ex:
                failures[(seed, region, is_midi)] = ex.args[0]
                continue
            midi_rows[get_genotype(seed)] = [row
                                             for row in rows
                                             if (NS5B_MAIN_END_POS - 2) < int(row['refseq.aa.pos'])]
    is_midi = False
    for (seed, region), rows in groupby(DictReader(amino_csv),
                                        itemgetter('seed', 'region')):
        low_coverage_message = None
        all_rows = main_rows = list(rows)
        try:
            check_coverage(region, main_rows)
        except LowCoverageError as ex:
            main_rows = []
            low_coverage_message = ex.args[0]
        if region.endswith('-NS5b'):
            genotype = get_genotype(seed)
            region_midi_rows = midi_rows.pop(genotype, None)
            if region_midi_rows is None:
                region_midi_rows = all_rows
                if is_midi_separate:
                    failures.setdefault((seed, region, True),
                                        NOTHING_MAPPED_MESSAGE)
                try:
                    check_coverage(region,
                                   region_midi_rows,
                                   start_pos=midi_start,
                                   end_pos=midi_end)
                    low_coverage_message = None
                except LowCoverageError:
                    region_midi_rows = []
            main_rows = combine_midi_rows(main_rows,
                                          region_midi_rows,
                                          seed)
        if low_coverage_message:
            failures[(seed, region, is_midi)] = low_coverage_message
        yield from main_rows

    # Check for MIDI regions that had no match.
    for genotype, rows in sorted(midi_rows.items()):
        region = rows[0]['region']
        seed = rows[0]['seed']
        failures.setdefault((seed, region, False), NOTHING_MAPPED_MESSAGE)
        yield from rows


def write_failure(fail_writer, seed, region, reason):
    reported_region = get_reported_region(region)
    genotype = get_genotype(seed)
    fail_writer.writerow(dict(seed=seed,
                              region=reported_region,
                              coord_region=region,
                              genotype=genotype,
                              reason=reason))


def combine_midi_rows(main_rows, midi_rows, seed):
    main_row_map = {int(row['refseq.aa.pos']): row
                    for row in main_rows}
    midi_row_map = {int(row['refseq.aa.pos']): row
                    for row in midi_rows}
    positions = sorted({pos
                        for pos in chain(main_row_map.keys(),
                                         midi_row_map.keys())})

    for midi_row in midi_row_map.values():
        midi_row['seed'] = seed  # Override MIDI seed with main seed.
    for pos in positions:
        main_row = main_row_map.get(pos)
        midi_row = midi_row_map.get(pos)
        if midi_row is None:
            if pos <= LAST_MAIN_POS:
                yield main_row
        elif main_row is None:
            yield midi_row
        elif (pos <= LAST_MAIN_POS and
              int(main_row['coverage']) > int(midi_row['coverage'])):
            yield main_row
        else:
            yield midi_row


def read_aminos(amino_rows,
                min_fraction,
                reported_regions=None,
                min_coverage=0,
                algorithms=None):
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
        for row in rows:
            counts = list(map(int, (row[f] for f in coverage_columns)))
            coverage = int(row['coverage'])
            pos = int(row['refseq.aa.pos'])
            while pos > len(aminos) + 1:
                aminos.append({})
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
        if algorithms is None:
            is_report_required = False
        else:
            std_length = get_std_length(region, genotype, algorithms)
            while len(aminos) < std_length:
                aminos.append({})
            asi_algorithm = algorithms.get(genotype)
            key_positions = asi_algorithm.get_gene_positions(region)
            if not region.endswith('NS5b'):
                is_report_required = all(aminos[pos-1] for pos in key_positions)
            else:
                whole_genome_positions = {pos for pos in key_positions if pos < MIDI_START_POS}
                midi_positions = key_positions - whole_genome_positions
                is_report_required = (
                    all(pos <= len(aminos) and aminos[pos-1]
                        for pos in whole_genome_positions) or
                    all(pos <= len(aminos) and aminos[pos-1]
                        for pos in midi_positions))
        # Need original region to look up wild type.
        yield AminoList(region, aminos, genotype, seed, is_report_required)


def get_algorithm_regions(algorithm):
    regions = []
    for region in algorithm.gene_def:
        if region == 'IN':
            regions.append('INT')
        elif region == 'CA':
            continue
        else:
            regions.append(region)
    return regions


def create_empty_aminos(region, genotype, seed, algorithms):
    std_length = get_std_length(region, genotype, algorithms)
    return AminoList(region,
                     [{}] * std_length,
                     genotype,
                     seed,
                     False)


def get_std_length(region, genotype, algorithms):
    algorithm = algorithms[genotype]
    std_name = 'IN' if region == 'INT' else region
    std_length = len(algorithm.stds[std_name])
    return std_length


def filter_aminos(all_aminos, algorithms):
    all_aminos = list(all_aminos)
    genotypes_needed = {(amino_list.genotype, amino_list.seed)
                        for amino_list in all_aminos
                        if amino_list.is_report_needed}
    selected_aminos = [amino_list
                       for amino_list in all_aminos
                       if (amino_list.genotype, amino_list.seed) in genotypes_needed]
    selected_regions = {(amino_list.genotype or '', amino_list.seed, amino_list.region)
                        for amino_list in selected_aminos}
    expected_regions = {(genotype or '', seed, region)
                        for genotype, seed in genotypes_needed
                        for region in get_algorithm_regions(algorithms[genotype])}
    missing_regions = sorted(expected_regions - selected_regions)
    selected_aminos += [create_empty_aminos(region, genotype or None, seed, algorithms)
                        for genotype, seed, region in missing_regions]
    selected_aminos.sort()
    return selected_aminos


def get_position_consensus(position_aminos):
    if not position_aminos:
        return '-'
    consensus = ''.join(position_aminos)
    if len(consensus) > 1:
        return '[' + ''.join(sorted(consensus)) + ']'
    return consensus


def write_consensus(resistance_consensus_writer, amino_list, alg_version):
    if resistance_consensus_writer is None:
        return
    consensus = ''.join(get_position_consensus(position_aminos)
                        for position_aminos in amino_list.aminos).rstrip('-')
    if not consensus:
        return
    stripped_consensus = consensus.lstrip('-')
    offset = len(consensus) - len(stripped_consensus)
    reported_region = get_reported_region(amino_list.region)
    resistance_consensus_writer.writerow(dict(seed=amino_list.seed,
                                              region=reported_region,
                                              coord_region=amino_list.region,
                                              version=alg_version,
                                              offset=offset,
                                              sequence=stripped_consensus))


def write_resistance(aminos,
                     resistance_csv,
                     mutations_csv,
                     algorithms=None,
                     resistance_consensus_csv=None):
    """ Calculate resistance scores and write them to files.

    :param list[AminoList] aminos: region is the coordinate
        reference name that this gene region was mapped to, and prevalance is a
        float between 0.0 and 1.0
    :param resistance_csv: open file to write resistance calls to, grouped by
        genotype, region, drug_class
    :param mutations_csv: open file to write mutations to, grouped by genotype,
        drug_class
    :param dict algorithms: {region: AsiAlgorithm}
    :param resistance_consensus_csv: open file to write resistance consensus to
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
         'genotype',
         'seed',
         'coord_region',
         'version'],
        lineterminator=os.linesep)
    resistance_writer.writeheader()
    mutations_writer = DictWriter(mutations_csv,
                                  ['drug_class',
                                   'mutation',
                                   'prevalence',
                                   'genotype',
                                   'region',
                                   'seed',
                                   'coord_region',
                                   'version'],
                                  lineterminator=os.linesep)
    mutations_writer.writeheader()
    if resistance_consensus_csv is None:
        resistance_consensus_writer = None
    else:
        resistance_consensus_writer = create_consensus_writer(
            resistance_consensus_csv)
    if algorithms is None:
        algorithms = load_asi()
    for genotype, genotype_aminos in groupby(aminos, attrgetter('genotype')):
        region_results = []
        for amino_list in genotype_aminos:
            asi = algorithms.get(genotype)
            write_consensus(resistance_consensus_writer,
                            amino_list,
                            asi.alg_version)
            if asi is None:
                continue
            region = amino_list.region
            if region == 'INT':
                region = 'IN'
            result = interpret(asi, amino_list.aminos, region)
            region_results.append((region, amino_list, asi.alg_version, result))
        if all(drug_result.level == HcvResistanceLevels.FAIL.level
               for region, amino_seq, alg_version, result in region_results
               for drug_result in result.drugs):
            continue
        for region, amino_list, alg_version, result in region_results:
            amino_seq = amino_list.aminos
            reported_region = get_reported_region(region)
            for drug_result in result.drugs:
                resistance_writer.writerow(dict(region=reported_region,
                                                drug_class=drug_result.drug_class,
                                                drug=drug_result.code,
                                                drug_name=drug_result.name,
                                                level_name=drug_result.level_name,
                                                level=drug_result.level,
                                                score=drug_result.score,
                                                genotype=genotype,
                                                seed=amino_list.seed,
                                                coord_region=amino_list.region,
                                                version=alg_version))
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
                                                   genotype=genotype,
                                                   region=reported_region,
                                                   seed=amino_list.seed,
                                                   coord_region=amino_list.region,
                                                   version=alg_version))


def write_nuc_mutations(nuc_csv: typing.TextIO,
                        nuc_mutations_csv: typing.TextIO):
    nuc_rows = DictReader(nuc_csv)
    mutations_writer = DictWriter(nuc_mutations_csv,
                                  ['seed',
                                   'region',
                                   'wt',
                                   'refseq_nuc_pos',
                                   'var',
                                   'prevalence',
                                   'ref_genome_pos'],
                                  lineterminator=os.linesep)
    mutations_writer.writeheader()
    for seed, seed_rows in groupby(nuc_rows, itemgetter('seed')):
        if seed != 'SARS-CoV-2-seed':
            continue
        landmark_reader = LandmarkReader.load()
        projects = ProjectConfig.loadDefault()
        for region_name, region_rows in groupby(seed_rows, itemgetter('region')):
            region = landmark_reader.get_gene(seed, region_name, drop_stop_codon=False)
            seed_seq = projects.getReference(seed)
            ref_seq = seed_seq[region['start']-1:region['end']]
            for row in region_rows:
                nuc_pos = int(row['refseq.nuc.pos'])
                wild_type = ref_seq[nuc_pos-1]
                coverage = int(row['coverage'])
                if coverage == 0:
                    continue
                for nuc in 'ACGT':
                    if nuc == wild_type:
                        continue
                    nuc_count = int(row[nuc])
                    prevalence = nuc_count / coverage
                    if prevalence >= 0.05:
                        mutations_writer.writerow(dict(seed=seed,
                                                       region=region_name,
                                                       wt=wild_type,
                                                       refseq_nuc_pos=nuc_pos,
                                                       var=nuc,
                                                       prevalence=prevalence,
                                                       ref_genome_pos=row['genome.pos']))


def create_consensus_writer(resistance_consensus_csv):
    resistance_consensus_writer = DictWriter(resistance_consensus_csv,
                                             ['seed',
                                              'region',
                                              'coord_region',
                                              'version',
                                              'offset',
                                              'sequence'],
                                             lineterminator=os.linesep)
    resistance_consensus_writer.writeheader()
    return resistance_consensus_writer


def interpret(asi, amino_seq, region):
    ref_seq = asi.stds[region]
    # TODO: Make this more general instead of only applying to missing MIDI.
    is_missing_midi = False
    if region.endswith('-NS5b'):
        pos = 0
        # starting at the end of the amino list, find the first amino that is not empty
        for pos, amino in enumerate(reversed(amino_seq)):
            if amino != {}:
                break
        if len(amino_seq) - pos <= LAST_MAIN_POS:
            is_missing_midi = True

    if is_missing_midi:
        amino_seq += [{}] * (len(ref_seq) - len(amino_seq))
    result = asi.interpret(amino_seq, region)

    if not is_missing_midi:
        for drug_result in result.drugs:
            if drug_result.level == HcvResistanceLevels.FAIL.level:
                break
        else:
            # No missing coverage, we're done.
            return result

    # At least some data found, assume wild type in gaps.
    new_amino_seq = []
    wild_type_seq = list(ref_seq)
    if region == 'IN':
        wild_type_seq[231] = 'D'  # Our reference is different from HIVdb's.
    for wild_type, old_aminos in zip_longest(wild_type_seq, amino_seq):
        if old_aminos:
            new_amino_seq.append(old_aminos)
        else:
            new_amino_seq.append({wild_type: 1.0})
    amino_seq = new_amino_seq

    new_result = asi.interpret(amino_seq, region)
    new_drug_results = {drug_result.code: drug_result
                        for drug_result in new_result.drugs}
    if is_missing_midi:
        allowed_levels = (HcvResistanceLevels.FAIL.level,
                          HcvResistanceLevels.RESISTANCE_LIKELY.level,
                          HcvResistanceLevels.NOT_INDICATED.level)
    elif asi.alg_name == 'HCV_RULES':
        allowed_levels = (HcvResistanceLevels.FAIL.level,
                          HcvResistanceLevels.UNKNOWN_MUTATIONS.level,
                          HcvResistanceLevels.RESISTANCE_POSSIBLE.level,
                          HcvResistanceLevels.RESISTANCE_LIKELY.level,
                          HcvResistanceLevels.NOT_INDICATED.level)
    else:
        allowed_levels = (HivResistanceLevels.FAIL.level,
                          HivResistanceLevels.LOW.level,
                          HivResistanceLevels.INTERMEDIATE.level,
                          HivResistanceLevels.HIGH.level)

    for drug_result in result.drugs:
        if drug_result.level == HcvResistanceLevels.FAIL.level or is_missing_midi:
            new_drug_result = new_drug_results[drug_result.code]
            if new_drug_result.level in allowed_levels:
                drug_result.level = new_drug_result.level
                drug_result.level_name = new_drug_result.level_name
                drug_result.score = new_drug_result.score
            else:
                drug_result.level = HcvResistanceLevels.FAIL.level
                drug_result.level_name = HcvResistanceLevels.FAIL.name
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
        if len(genotype) > 1 and genotype[0] in genotypes:
            backup_genotype = genotype[0]
        else:
            backup_genotype = None
        asi = AsiAlgorithm(rules_yaml=hcv_rules,
                           genotype=genotype,
                           backup_genotype=backup_genotype)
        algorithms[genotype] = asi

    return algorithms


def write_failures(failures, filtered_aminos, fail_csv):
    fail_writer = create_fail_writer(fail_csv)
    reported_keys = set()  # {(seed, region)}
    for amino_list in filtered_aminos:
        reported_keys.add((amino_list.seed, amino_list.region))
        if not any(amino_list.aminos):
            failure_key = (amino_list.seed, amino_list.region, False)
            failures.setdefault(failure_key, NOTHING_MAPPED_MESSAGE)
    for (seed, region, is_midi), message in sorted(failures.items()):
        if not reported_keys or (seed, region) in reported_keys:
            if is_midi:
                message = 'MIDI: ' + message
            write_failure(fail_writer, seed, region, message)


def report_resistance(amino_csv,
                      midi_amino_csv,
                      nuc_csv,
                      resistance_csv,
                      mutations_csv,
                      nuc_mutations_csv,
                      fail_csv,
                      region_choices=None,
                      resistance_consensus_csv=None):
    if region_choices is None:
        selected_regions = REPORTED_REGIONS
    else:
        selected_regions = select_reported_regions(region_choices, REPORTED_REGIONS)
    failures = {}
    amino_rows = combine_aminos(amino_csv, midi_amino_csv, failures)
    algorithms = load_asi()
    aminos = read_aminos(amino_rows,
                         MIN_FRACTION,
                         selected_regions,
                         MIN_COVERAGE,
                         algorithms)
    filtered_aminos = filter_aminos(aminos, algorithms)
    write_failures(failures, filtered_aminos, fail_csv)
    write_resistance(filtered_aminos,
                     resistance_csv,
                     mutations_csv,
                     algorithms,
                     resistance_consensus_csv)
    write_nuc_mutations(nuc_csv, nuc_mutations_csv)


def main():
    args = parse_args()
    report_resistance(args.main_amino_csv,
                      args.midi_amino_csv,
                      args.main_nuc_csv,
                      args.resistance_csv,
                      args.mutations_csv,
                      args.nuc_mutations_csv,
                      args.resistance_fail_csv)


if __name__ == '__main__':
    main()
