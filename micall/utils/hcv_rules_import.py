import re
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter, FileType
from collections import namedtuple, defaultdict
from functools import partial

import os
from operator import itemgetter

import yaml
from pyvdrm.vcf import MutationSet
from openpyxl import load_workbook

from micall.core.project_config import ProjectConfig

READY_TABS = ('NS3_GT1a', 'NS3_GT1b')
PHENOTYPE_SCORES = {'likely susceptible': 0,
                    'resistance possible': 4,
                    'resistance likely': 8,
                    'effect unknown': 'effect unknown'}
DATA_CHANGES = {'NS3_GT1a': {'A65': ('V132I', 'I132I'),
                             'A104': ('I/V170F', 'I170F'),
                             'A105': ('I/V170T', 'I170T')},
                'NS3_GT1b': {'A138': ('NS4A V23A', None)}}
RuleSet = namedtuple(
    'RuleSet',
    ['region',
     'genotype',
     'drug_name',
     'first_column',
     'phenotype_column',
     'last_column',
     'mutations'])
Reference = namedtuple('Reference', 'name sequence')


def parse_args():
    parser = ArgumentParser(description='Read HCV rules from a spreadsheet.',
                            formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('spreadsheet',
                        type=partial(load_workbook, read_only=False))
    parser.add_argument('rules_yaml',
                        type=FileType('w'),
                        help='YAML file to write new rules to')
    return parser.parse_args()


class SplitDumper(yaml.SafeDumper):
    def write_plain(self, text, split=True):
        delimiter = ','
        if split:
            pieces = text.split(delimiter)
        else:
            pieces = [text]
        buffer = ''
        for i, piece in enumerate(pieces):
            if i > 0:
                buffer += delimiter
            if self.column-1 + len(buffer) + len(piece) <= self.best_width:
                buffer += piece
            else:
                super(SplitDumper, self).write_plain(buffer, split)
                self.write_indent()
                buffer = piece
        super(SplitDumper, self).write_plain(buffer)


def main():
    args = parse_args()
    references = load_references()
    all_rule_sets = []
    make_data_changes(args.spreadsheet)
    for ws in args.spreadsheet:
        if ws.title not in READY_TABS:
            continue
        print(ws.title)
        if not ws.title.startswith('NS'):
            print('Skipped.')
            continue
        rule_sets = read_rule_sets(ws, references)
        all_rule_sets.extend(rule_sets)
    write_rules(all_rule_sets, references, args.rules_yaml)


def make_data_changes(wb):
    for sheet_name, changes in DATA_CHANGES.items():
        ws = wb[sheet_name]
        for coordinate, (old, new) in changes.items():
            cell = ws[coordinate]
            assert cell.value == old, (sheet_name, coordinate, cell.value)
            ws[coordinate] = new


def read_rule_sets(ws, references):
    header_rows = []
    monitored_drugs = set()
    rule_sets = None
    previous_row = None
    for row in ws.rows:
        if rule_sets is None:
            if row[0].value and row[0].value.startswith('WT'):
                rule_sets = find_rule_sets(ws, header_rows)
                find_phenotype_scores(row, rule_sets)
            else:
                header_rows.append(row)
        elif row[0].value is not None:
            find_phenotype_scores(row, rule_sets)
        else:
            # In footer rows.
            if previous_row is not None:
                for col, cell in enumerate(previous_row, 1):
                    label = str(cell.value).lower()
                    if label.startswith('positions monitored'):
                        drug_name = monitor_positions(ws,
                                                      rule_sets,
                                                      references,
                                                      col,
                                                      row[col-1].value)
                        monitored_drugs.add(drug_name)
            previous_row = row
    if rule_sets is None:
        raise ValueError('No mutation started with WT in {}.'.format(ws.title))
    all_drugs = {rule_set.drug_name for rule_set in rule_sets}
    missing_drugs = all_drugs - monitored_drugs
    if missing_drugs:
        raise ValueError("No 'monitored positions' label for {} in {}.".format(
            ', '.join(sorted(missing_drugs)),
            ws.title))
    return rule_sets


def monitor_positions(worksheet, rule_sets, references, col, positions):
    rule_set = None
    for rule_set in rule_sets:
        if rule_set.first_column <= col <= rule_set.last_column:
            if positions == 'None':
                break
            if positions is None:
                raise ValueError(
                    "No list of monitored positions for {} in {}.".format(
                        rule_set.drug_name,
                        worksheet.title))
            reference = references[(rule_set.genotype.upper(), rule_set.region)]
            try:
                seen_variants = {int(pos): set()
                                 for pos in positions.split(',')}
            except ValueError:
                raise ValueError(
                    'Invalid monitored position for {} in {}: {!r}.'.format(
                        rule_set.drug_name,
                        worksheet.title,
                        positions))
            for mutation in rule_set.mutations:
                match = re.match(r'^([A-Z])(\d+)([A-Z]+)$', mutation)
                assert match, mutation
                wild_type = match.group(1)
                pos = int(match.group(2))
                assert wild_type == reference.sequence[pos-1], (reference.name,
                                                                mutation)
                variants = match.group(3)
                pos_variants = seen_variants.get(pos)
                if pos_variants is not None:
                    pos_variants.update(variants)
            for pos, pos_variants in seen_variants.items():
                wild_type = reference.sequence[pos - 1]
                mutation = '{}{}!{}{}'.format(wild_type,
                                              pos,
                                              ''.join(sorted(pos_variants)),
                                              wild_type)
                rule_set.mutations[mutation] = 'Effect unknown'
            break
    return rule_set.drug_name


def find_phenotype_scores(row, rule_sets):
    mutation = row[0].value
    if '+' in mutation:
        return
    for rule_set in rule_sets:
        phenotype_cell = row[rule_set.phenotype_column-1]
        phenotype = phenotype_cell.value
        if phenotype is not None:
            try:
                score = PHENOTYPE_SCORES[phenotype.strip(' ?*').lower()]
            except Exception:
                raise ValueError('Unknown phenotype for {} {}: {!r}.'.format(
                    rule_set.drug_name,
                    mutation,
                    phenotype))
            if score:
                rule_set.mutations[mutation] = score


def create_rule_set(sheet_name,
                    drug_name,
                    first_column,
                    phenotype_column,
                    last_column):
    region, genotype = sheet_name.split('_GT')
    short_name = get_short_drug_name(drug_name)
    return RuleSet(region,
                   genotype,
                   short_name,
                   first_column,
                   phenotype_column,
                   last_column,
                   {})


def load_drug_codes():
    config_path = os.path.normpath(os.path.join(__file__,
                                                '..',
                                                '..',
                                                'resistance',
                                                'genreport.yaml'))
    with open(config_path) as config_file:
        report_config = yaml.safe_load(config_file)
        drug_codes = {}
        for section in report_config:
            known_regions = section.get('known_regions', [])
            if 'NS3' not in known_regions:
                continue
            for drug_class in section['known_drugs'].values():
                for code, name in drug_class:
                    short_name = get_short_drug_name(name)
                    drug_codes[short_name] = code
    return drug_codes


def load_references():
    projects = ProjectConfig.loadDefault()
    references = {}  # {(genotype, region): Reference}
    for ref_name, sequence in projects.getAllReferences().items():
        match = re.match(r'HCV(.*?)-.*-([^-]+)$', ref_name)
        if match:
            genotype = match.group(1)
            region = match.group(2)
            if region in ('NS3', 'NS5a', 'NS5b'):
                reference = Reference(ref_name, sequence)
                references[(genotype, region)] = reference
                if genotype == '6':
                    references[('6E', region)] = reference
    return references


def get_short_drug_name(name):
    return name.split()[0]


def find_rule_sets(ws, header_rows):
    drug_row_num = len(header_rows) - 1
    drug_ranges = []
    for cell_range in ws.merged_cells.ranges:
        if cell_range.max_row == drug_row_num:
            drug_ranges.append(cell_range)
    rule_sets = []
    drug_row = header_rows[-2]
    for cell in drug_row:
        if cell.value:
            for drug_range in drug_ranges:
                if cell.coordinate in drug_range:
                    for col in range(drug_range.min_col,
                                     drug_range.max_col + 1):
                        header = header_rows[-1][col-1].value
                        if header == 'Phenotype':
                            rule_sets.append(create_rule_set(ws.title,
                                                             cell.value,
                                                             drug_range.min_col,
                                                             col,
                                                             drug_range.max_col))
                            break
                    else:
                        template = 'No Phenotype column between columns ' \
                                   '{} and {} of {}.'
                        raise ValueError(
                            template.format(drug_range.min_col,
                                            drug_range.max_col,
                                            ws.title))
                    break
    return rule_sets


def format_score(score):
    try:
        if '"' in score:
            raise ValueError('Text score contained quotes: {}.'.format(score))
        return '"{}"'.format(score)
    except TypeError:
        return score


def write_rules(rule_sets, references, rules_file):
    drug_codes = load_drug_codes()
    drug_summaries = {}
    for rule_set in rule_sets:
        drug_name = rule_set.drug_name
        try:
            drug_summary = drug_summaries[drug_name]
        except KeyError:
            drug_code = drug_codes[drug_name]
            drug_summary = drug_summaries[drug_name] = dict(name=drug_name,
                                                            code=drug_code,
                                                            genotypes=[])
        positions = defaultdict(dict)  # {pos: {score: MutationSet}}
        for mutation, score in rule_set.mutations.items():
            if '+' in mutation or ' ' in mutation:
                continue
            try:
                new_mutation_set = MutationSet(mutation)
            except ValueError:
                raise ValueError('Unable to parse mutation: ' + mutation)

            pos_scores = positions[new_mutation_set.pos]
            old_mutation_set = pos_scores.get(score)
            if old_mutation_set is not None:
                new_mutation_set = MutationSet(
                    wildtype=old_mutation_set.wildtype,
                    pos=old_mutation_set.pos,
                    mutations=(old_mutation_set.mutations |
                               new_mutation_set.mutations))
            pos_scores[score] = new_mutation_set
        score_terms = sorted((mutation_set, format_score(score))
                             for pos, pos_scores in positions.items()
                             for score, mutation_set in pos_scores.items())
        score_formula = 'SCORE FROM ( {} )'.format(', '.join(
            '{} => {}'.format(mutation_set, score)
            for mutation_set, score in score_terms))
        genotype = rule_set.genotype.upper()
        reference = references[(genotype, rule_set.region)]
        reference_name = reference.name
        drug_summary['genotypes'].append(dict(genotype=genotype,
                                              region=rule_set.region,
                                              reference=reference_name,
                                              rules=score_formula))
    drugs = sorted(drug_summaries.values(), key=itemgetter('code'))
    yaml.dump(drugs, rules_file, default_flow_style=False, Dumper=SplitDumper)


if __name__ == '__main__':

    main()
