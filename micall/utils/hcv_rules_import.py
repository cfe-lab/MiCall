from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter, FileType, Namespace
from collections import namedtuple, defaultdict
from functools import partial
from itertools import product, groupby
from operator import itemgetter, attrgetter
import os
import re
import sys

import yaml
from pyvdrm.hcvr import HCVR
from pyvdrm.vcf import MutationSet, VariantCalls
from openpyxl import load_workbook

from micall.core.project_config import ProjectConfig

READY_TABS = ('NS3_GT1a',
              'NS3_GT1b',
              'NS3_GT2',
              'NS3_GT3',
              'NS3_GT4',
              'NS3_GT5',
              'NS3_GT6',
              'NS5A_GT1a',
              'NS5A_GT1b',
              'NS5A_GT2',
              'NS5A_GT3',
              'NS5A_GT4',
              'NS5A_GT5',
              'NS5A_GT6')
PHENOTYPE_SCORES = {'likely susceptible': 0,
                    'resistance possible': 4,
                    'resistance likely': 8,
                    'effect unknown': 'effect unknown'}
HCV_REGIONS = ('NS3', 'NS5a', 'NS5b')
RuleSet = namedtuple(
    'RuleSet',
    ['region',
     'genotype',
     'drug_name',
     'first_column',
     'phenotype_column',
     'last_column',
     'mutations',  # {condition: score}
     'fold_shifts'])  # {condition: shift_text}
RuleSet.__new__.__defaults__ = (None, )
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


class WorksheetReader:
    def __init__(self, worksheets, *footer_readers):
        self.worksheets = worksheets
        self.missing_wild_types = []
        self.missing_monitored_positions = []
        self.footer_readers = footer_readers

    def __iter__(self):
        for ws in self.worksheets:
            yield from self._iter_worksheet(ws)

    def _iter_worksheet(self, ws):
        wild_type_row, footer_row = self._find_row_limits(ws)
        if wild_type_row is None:
            self.missing_wild_types.append(ws.title)
            return
        header_row = wild_type_row - 1
        drug_row = wild_type_row - 2
        for cell_range in ws.merged_cells:
            if cell_range.min_row != drug_row:
                continue
            cell = ws.cell(cell_range.min_row, cell_range.min_col)
            if cell.value is None:
                continue
            drug_name, *_ = cell.value.split()
            column_headings = {
                col: self._format_heading(ws.cell(header_row, col).value)
                for col in range(cell_range.min_col, cell_range.max_col + 1)
                if ws.cell(header_row, col).value is not None}
            column_headings[1] = 'mutation'
            section = Namespace(drug_name=drug_name, sheet_name=ws.title)
            if not self._read_footer(ws, cell_range, footer_row, section):
                section.not_indicated = True
                yield Namespace(section=section)
                continue
            for row_num in range(wild_type_row, footer_row):
                if ws.cell(row_num, 1).font.strike:
                    continue
                fields = {heading: ws.cell(row_num, col).value
                          for col, heading in column_headings.items()}
                for field in fields:
                    value = fields[field]
                    if value is not None:
                        fields[field] = str(value).strip()
                fields['mutation'] = re.sub(r'del$', 'd', fields['mutation'])
                fields['mutation'] = re.sub(r'ins$', 'i', fields['mutation'])
                entry = Namespace(section=section, **fields)
                if entry.phenotype is not None:
                    yield entry

    def _read_footer(self, ws, cell_range, footer_row, section):
        footer_coordinates = product(
            range(footer_row, ws.max_row + 1),
            range(cell_range.min_col, cell_range.max_col + 1))
        for row, col in footer_coordinates:
            value = str(ws.cell(row, col).value).lower().strip()
            if value == 'not indicated in canada':
                return False
        for footer_reader in self.footer_readers:
            footer_coordinates = product(
                range(footer_row, ws.max_row + 1),
                range(cell_range.min_col, cell_range.max_col + 1))
            footer_reader.read(ws, footer_coordinates, section)
        return True

    @staticmethod
    def _find_row_limits(ws):
        wild_type_row = None
        footer_row = ws.max_row + 1
        for row_num in range(1, ws.max_row + 1):
            mutation = ws.cell(row_num, 1).value
            if wild_type_row is None:
                if mutation and mutation.startswith('WT'):
                    wild_type_row = row_num
            else:
                if mutation is None:
                    footer_row = row_num
                    break
        return wild_type_row, footer_row

    @staticmethod
    def _format_heading(text):
        heading = text.lower()
        heading = re.sub(r'[^a-z0-9]', '_', heading)
        return heading.strip('_')

    def write_errors(self, report):
        if self.missing_wild_types:
            report.write('No wild type found in {}.\n'.format(
                ', '.join(self.missing_wild_types)))
        for footer_reader in self.footer_readers:
            footer_reader.write_errors(report)


class MonitoredPositionsReader:
    def __init__(self):
        self.missing = []

    def read(self, ws, footer_coordinates, section):
        monitored_positions = None
        for row_num, col_num in footer_coordinates:
            cell = ws.cell(row_num, col_num)
            label = cell.value and str(cell.value).lower()
            if label and label.startswith('positions monitored'):
                monitored_positions = ws.cell(row_num + 1, col_num).value
                break
        if monitored_positions is None:
            monitored_positions = []
            self.missing.append(
                '{} in {}'.format(section.drug_name, section.sheet_name))
        else:
            monitored_positions = str(monitored_positions)
            if monitored_positions.lower() == 'none':
                monitored_positions = []
            else:
                monitored_positions = list(map(
                    int,
                    str(monitored_positions).split(',')))
        section.monitored_positions = monitored_positions

    def write_errors(self, report):
        if self.missing:
            report.write('No monitored positions for {}.\n'.format(
                ', '.join(self.missing)))


class FoldRangesReader:
    def __init__(self):
        self.invalid = []
        self.lower_pattern = \
            r'([</=]*) *([\d.]+)(-([\d.]+))?[x ]*FS( range)?, likely susceptibi?le$'
        self.lower_score = PHENOTYPE_SCORES['likely susceptible']
        self.upper_pattern = \
            r'([>/=]*) *([\d.]+)(-([\d.]+))?[x ]*FS( range)?, resistance likely$'
        self.upper_score = PHENOTYPE_SCORES['resistance likely']
        self.middle_score = PHENOTYPE_SCORES['resistance possible']
        self.level_descriptions = {
            score: description
            for description, score in PHENOTYPE_SCORES.items()}

    def read(self, ws, footer_coordinates, section):
        label_positions = {}  # {label: (row, col)}
        for row_num, col_num in footer_coordinates:
            label = ws.cell(row_num, col_num).value
            if not label:
                continue
            label = str(label).lower()
            label = label.replace('virto', 'vitro')
            label = label.replace('susecptibility', 'susceptibility')
            label = label.strip(' -:*')
            label_positions[label] = (row_num, col_num)
        lower_fold, upper_fold = self.find_fold_range(
            'in vitro drug susceptibility',
            label_positions,
            section,
            ws)
        section.lower_fold = lower_fold
        section.upper_fold = upper_fold
        lower_range_fold, upper_range_fold = self.find_fold_range(
                'range fs values',
                label_positions,
                section,
                ws)
        if lower_range_fold is None:
            section.lower_range_fold = section.lower_fold
            section.upper_range_fold = section.upper_fold
        else:
            section.lower_range_fold = lower_range_fold
            section.upper_range_fold = upper_range_fold

    def find_fold_range(self, range_label, label_positions, section, ws):
        row_num, col_num = label_positions.get(range_label,
                                               (None, None))
        if row_num is None:
            return None, None
        lower_text = ws.cell(row_num + 1, col_num).value
        if str(lower_text).lower().strip(' -*') == 'absolute fs values':
            row_num += 1
            lower_text = ws.cell(row_num + 1, col_num).value
        upper_text = ws.cell(row_num + 3, col_num).value
        lower_match = re.match(self.lower_pattern,
                               lower_text,
                               flags=re.IGNORECASE)
        upper_match = re.match(self.upper_pattern,
                               upper_text,
                               flags=re.IGNORECASE)
        if upper_match and lower_match:
            lower_fold = float(lower_match.group(4)
                               or lower_match.group(2))
            upper_fold = float(upper_match.group(2))
        else:
            lower_fold = upper_fold = None
            if lower_match is None:
                self.invalid.append(
                    "Invalid lower fold shift of {!r} for {} in {}.".format(
                        lower_text,
                        section.drug_name,
                        section.sheet_name))
            if upper_match is None:
                self.invalid.append(
                    "Invalid upper fold shift of {!r} for {} in {}.".format(
                        upper_text,
                        section.drug_name,
                        section.sheet_name))
        return lower_fold, upper_fold

    def write_errors(self, report):
        for message in self.invalid:
            print(message, file=report)


class RulesWriter:
    def __init__(self, rules_file, errors_file, references):
        self.rules_file = rules_file
        self.errors_file = errors_file
        self.references = references
        self.unknown_drugs = []
        self.invalid_mutations = []
        self.invalid_mutation_positions = []
        self.combination_changes = []
        self.invalid_fold_shifts = []
        self.phenotype_changes = []
        self.invalid_phenotypes = []
        self.invalid_positions = []
        self.bad_wild_types = []

        # {(sheet, drug_name): [mutation]}
        self.missing_fold_shifts = defaultdict(list)

        self.score_phenotypes = {
            score: phenotype
            for phenotype, score in PHENOTYPE_SCORES.items()}

    def write(self, entries):
        drug_codes = load_drug_codes()
        drug_summaries = {}
        for section, section_entries in groupby(entries, attrgetter('section')):
            try:
                drug_summary = self._find_drug_summary(drug_summaries,
                                                       drug_codes,
                                                       section)
            except KeyError:
                continue
            sheet_region, genotype = section.sheet_name.split('_GT')
            lower_region = sheet_region.lower()
            for known_region in HCV_REGIONS:
                if known_region.lower() == lower_region:
                    region = known_region
                    break
            else:
                region = sheet_region
            genotype = genotype.upper()
            reference = self.references[(genotype, region)]
            reference_name = reference.name
            positions = defaultdict(dict)  # {pos: {score: MutationSet}}
            combinations = {}  # {mutation: score}
            if getattr(section, 'not_indicated', False):
                # noinspection PyTypeChecker
                positions[None]['Not indicated'] = 'TRUE'
            else:
                for entry in section_entries:
                    self._score_mutation(entry,
                                         section,
                                         positions,
                                         combinations,
                                         reference.sequence)
            self._check_combinations(combinations, positions, section)
            self._monitor_positions(section, positions, reference)
            score_formula = build_score_formula(positions)
            drug_summary['genotypes'].append(dict(genotype=genotype,
                                                  region=region,
                                                  reference=reference_name,
                                                  rules=score_formula))
        self._check_not_indicated(drug_summaries)
        drugs = sorted(drug_summaries.values(), key=itemgetter('code'))
        yaml.dump(drugs, self.rules_file, default_flow_style=False, Dumper=SplitDumper)

    def _check_not_indicated(self, drug_summaries):
        all_genotypes = {genotype['genotype']
                         for drug_summary in drug_summaries.values()
                         for genotype in drug_summary['genotypes']}
        positions = {None: {'Not indicated': 'TRUE'}}
        score_formula = build_score_formula(positions)
        for drug_summary in drug_summaries.values():
            drug_genotypes = {genotype['genotype']
                              for genotype in drug_summary['genotypes']}
            missing_genotypes = all_genotypes - drug_genotypes

            drug_region, = {genotype['region']
                            for genotype in drug_summary['genotypes']}
            for genotype_name in missing_genotypes:
                reference = self.references[(genotype_name, drug_region)]
                drug_summary['genotypes'].append(dict(genotype=genotype_name,
                                                      region=drug_region,
                                                      reference=reference.name,
                                                      rules=score_formula))
            drug_summary['genotypes'].sort(key=itemgetter('genotype'))

    def _score_mutation(self,
                        entry,
                        section,
                        positions,
                        combinations,
                        reference):
        mutation = entry.mutation
        phenotype = entry.phenotype.lower()
        try:
            score = PHENOTYPE_SCORES[phenotype]
        except KeyError:
            self.invalid_phenotypes.append(
                '{}, {}, {}: {}'.format(section.sheet_name,
                                        section.drug_name,
                                        entry.mutation,
                                        entry.phenotype))
            return

        self._check_fold_shift(entry, section, phenotype)
        if not score:
            return
        if mutation.startswith('WT'):
            # noinspection PyTypeChecker
            positions[None][score] = 'TRUE'
            return
        if '+' in mutation or ' ' in mutation:
            combinations[mutation] = score
            return
        try:
            new_mutation_set = MutationSet(mutation)
        except ValueError as ex:
            self.invalid_mutations.append('{}: {} ({})'.format(
                section.sheet_name,
                mutation,
                ex))
            return
        max_pos = len(reference)
        if new_mutation_set.pos > max_pos:
            self.invalid_mutation_positions.append(
                f'{section.sheet_name}: {section.drug_name} '
                f'{new_mutation_set} (max {max_pos})')
            return
        expected_wild_type = reference[new_mutation_set.pos-1]
        if expected_wild_type != new_mutation_set.wildtype:
            self.bad_wild_types.append(
                f'{section.sheet_name}: {new_mutation_set} in '
                f'{section.drug_name} expected {expected_wild_type}')
            return
        pos_scores = positions[new_mutation_set.pos]
        old_mutation_set = pos_scores.get(score)
        if old_mutation_set is not None:
            new_mutation_set = MutationSet(
                wildtype=old_mutation_set.wildtype,
                pos=old_mutation_set.pos,
                mutations=(old_mutation_set.mutations |
                           new_mutation_set.mutations))
        pos_scores[score] = new_mutation_set

    def _check_fold_shift(self, entry, section, phenotype):
        if getattr(section, 'lower_fold', None) is None:
            # Range wasn't found, so nothing to check.
            return
        if entry.fold_shift is None:
            mutations = self.missing_fold_shifts[
                (section.sheet_name, section.drug_name)]
            mutations.append(entry.mutation)
            return
        fold_shift_text = str(entry.fold_shift)
        match = re.match(
            r'([<>/=]*) *(([\d.,]+)x?|([\d.,]+)x? *- *>? *([\d.,]+)x)$',
            fold_shift_text,
            flags=re.IGNORECASE)
        if match is None:
            self.invalid_fold_shifts.append(
                f'{section.sheet_name}: {section.drug_name} {entry.mutation} '
                f'{fold_shift_text!r}')
            return
        if match.group(3):
            expected_score = self._calculate_fold_shift_score(
                match.group(3),
                match.group(1),
                section.lower_fold,
                section.upper_fold)
        else:
            expected_score = self._calculate_fold_shift_range_score(
                match.group(5),
                section.lower_range_fold,
                section.upper_range_fold)
        clinical_ras = entry.clinical_ras and str(entry.clinical_ras).lower()
        if clinical_ras == 'yes':
            expected_score = min(8, expected_score + 4)
        expected_phenotype = self.score_phenotypes[expected_score]
        if phenotype != expected_phenotype:
            self.phenotype_changes.append(
                f'{section.sheet_name}: {section.drug_name} {entry.mutation} '
                f'{expected_phenotype} => {phenotype}')

    @staticmethod
    def _calculate_fold_shift_score(fold_shift_text,
                                    comparison,
                                    lower_fold,
                                    upper_fold):
        fold_shift = float(fold_shift_text.replace(',', ''))
        if comparison == '<':
            fold_shift -= 0.1
        elif comparison == '>':
            fold_shift += 0.1
        if fold_shift < lower_fold:
            expected_score = 0
        elif fold_shift > upper_fold:
            expected_score = 8
        else:
            expected_score = 4
        return expected_score

    @staticmethod
    def _calculate_fold_shift_range_score(upper_text,
                                          lower_fold_limit,
                                          upper_fold_limit):
        upper_fold_shift = float(upper_text.replace(',', ''))
        if upper_fold_shift <= lower_fold_limit:
            expected_score = 0
        elif upper_fold_shift > upper_fold_limit:
            expected_score = 8
        else:
            expected_score = 4
        return expected_score

    def _check_combinations(self, combinations, positions, section):
        for combination, combination_score in combinations.items():
            try:
                component_score = calculate_component_score(combination, positions)
            except ValueError as ex:
                self.invalid_mutations.append('{}: {} ({})'.format(
                    section.sheet_name,
                    combination,
                    ex))
                continue

            component_score = min(component_score, 8)
            if component_score != combination_score:
                self.combination_changes.append('{}: {}: {} => {}'.format(
                    section.sheet_name,
                    combination,
                    component_score,
                    combination_score))

    def _find_drug_summary(self, drug_summaries, drug_codes, section):
        drug_name = section.drug_name
        try:
            drug_summary = drug_summaries[drug_name]
        except KeyError:
            drug_code = drug_codes.get(drug_name)
            if drug_code is None:
                self.unknown_drugs.append('{}: {}'.format(
                    section.sheet_name,
                    drug_name))
                raise
            drug_summary = drug_summaries[drug_name] = dict(name=drug_name,
                                                            code=drug_code,
                                                            genotypes=[])
        return drug_summary

    def _monitor_positions(self, section, position_scores, reference):
        monitored_positions = getattr(section, 'monitored_positions', [])
        max_pos = len(reference.sequence)
        for pos in monitored_positions:
            scores = position_scores[pos]
            seen_variants = {mutation.variant
                             for mutation_set in scores.values()
                             for mutation in mutation_set}
            try:
                wild_type = reference.sequence[pos - 1]
            except IndexError:
                self.invalid_positions.append(
                    f'{section.sheet_name}: {section.drug_name} {pos} '
                    f'(max {max_pos})')
                continue
            mutation_set = MutationSet(
                '{}{}!{}{}'.format(wild_type,
                                   pos,
                                   ''.join(seen_variants),
                                   wild_type))
            position_scores[pos]['Effect unknown'] = mutation_set

    def _write_error_group(self, header, errors):
        if not errors:
            return
        if len(errors) == 1:
            self.errors_file.write('{}: {}.\n'.format(header, errors[0]))
        else:
            self.errors_file.write('{}s:\n  {}\n'.format(header,
                                                         '\n  '.join(errors)))

    def write_errors(self):
        missing_fold_shift_reports = []
        for (sheet_name, drug_name), mutations in sorted(
                self.missing_fold_shifts.items()):
            mutation_report = ', '.join(mutations)
            missing_fold_shift_reports.append(
                f'{sheet_name}: {drug_name} {mutation_report}')
        self._write_error_group('Unknown drug', self.unknown_drugs)
        self._write_error_group('Invalid mutation', self.invalid_mutations)
        self._write_error_group('Invalid mutation position',
                                self.invalid_mutation_positions)
        self._write_error_group('Invalid phenotype', self.invalid_phenotypes)
        self._write_error_group('Mismatched wild type', self.bad_wild_types)
        self._write_error_group('Invalid monitored position',
                                self.invalid_positions)
        self._write_error_group('Combination change', self.combination_changes)
        self._write_error_group('Missing fold shift', missing_fold_shift_reports)
        self._write_error_group('Invalid fold shift', self.invalid_fold_shifts)
        self._write_error_group('Phenotype change', self.phenotype_changes)


def main():
    args = parse_args()
    references = load_references()
    worksheets = [ws
                  for ws in args.spreadsheet
                  if ws.title in READY_TABS]
    reader = WorksheetReader(worksheets,
                             MonitoredPositionsReader(),
                             FoldRangesReader())
    writer = RulesWriter(args.rules_yaml, sys.stderr, references)
    writer.write(reader)
    reader.write_errors(sys.stderr)
    writer.write_errors()


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
            if region in HCV_REGIONS:
                reference = Reference(ref_name, sequence)
                references[(genotype, region)] = reference
                if genotype == '6':
                    references[('6E', region)] = reference
    return references


def get_short_drug_name(name):
    return name.split()[0]


def format_score(score):
    try:
        if '"' in score:
            raise ValueError('Text score contained quotes: {}.'.format(score))
        return '"{}"'.format(score)
    except TypeError:
        return score


def calculate_component_score(combination, positions):
    variant_calls = VariantCalls(combination.replace('+', ' '))
    combination_positions = {}
    for mutation_set in variant_calls:
        position_scores = positions[mutation_set.pos]
        if position_scores:
            combination_positions[mutation_set.pos] = position_scores
    if not combination_positions:
        return 0
    score_formula = build_score_formula(combination_positions)
    try:
        rule = HCVR(score_formula)
    except Exception as ex:
        raise ValueError('Bad formula for {}.'.format(combination)) from ex
    component_score = rule(variant_calls)
    return component_score


def build_score_formula(positions):
    wild_type_score = positions.pop(None, {})
    score_terms = sorted((mutation_set, format_score(score))
                         for pos, pos_scores in positions.items()
                         for score, mutation_set in pos_scores.items())
    for score, condition in sorted(wild_type_score.items()):
        score_terms.insert(0, (condition, format_score(score)))
    if not score_terms:
        score_terms.append(('TRUE', '0'))
    score_formula = 'SCORE FROM ( {} )'.format(', '.join(
        '{} => {}'.format(mutation_set, score)
        for mutation_set, score in score_terms))
    return score_formula


if __name__ == '__main__':
    main()
