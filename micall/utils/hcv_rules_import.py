from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter, FileType, Namespace
from collections import namedtuple, defaultdict, Counter
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
        self.lower_pattern = r'<([\d.]+)x? ?FS, likely susceptible$'
        self.lower_score = PHENOTYPE_SCORES['likely susceptible']
        self.upper_pattern = r'>([\d.]+)x? ?FS, resistance likely$'
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
        row_num, col_num = label_positions.get('absolute fs values',
                                               (None, None))
        if row_num is None:
            row_num, col_num = label_positions.get('in vitro drug susceptibility',
                                                   (None, None))

        lower_text = ws.cell(row_num + 1, col_num).value
        upper_text = ws.cell(row_num + 3, col_num).value
        lower_match = re.match(self.lower_pattern,
                               lower_text,
                               flags=re.IGNORECASE)
        upper_match = re.match(self.upper_pattern,
                               upper_text,
                               flags=re.IGNORECASE)
        if upper_match and lower_match:
            section.lower_fold = float(lower_match.group(1))
            section.upper_fold = float(upper_match.group(1))
        else:
            section.lower_fold = section.upper_fold = None
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
        self.invalid_phenotypes = []
        self.invalid_positions = []
        self.bad_wild_types = []

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
        drugs = sorted(drug_summaries.values(), key=itemgetter('code'))
        yaml.dump(drugs, self.rules_file, default_flow_style=False, Dumper=SplitDumper)

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
        self._write_error_group('Unknown drug', self.unknown_drugs)
        self._write_error_group('Invalid mutation', self.invalid_mutations)
        self._write_error_group('Invalid mutation position',
                                self.invalid_mutation_positions)
        self._write_error_group('Invalid phenotype', self.invalid_phenotypes)
        self._write_error_group('Mismatched wild type', self.bad_wild_types)
        self._write_error_group('Invalid monitored position',
                                self.invalid_positions)
        self._write_error_group('Combination change', self.combination_changes)


class FoldShiftChecker:
    def __init__(self):
        self.lower_pattern = r'<(\d+)x? FS, likely susceptible$'
        self.lower_score = PHENOTYPE_SCORES['likely susceptible']
        self.upper_pattern = r'>(\d+)x? FS, resistance likely$'
        self.upper_score = PHENOTYPE_SCORES['resistance likely']
        self.middle_score = PHENOTYPE_SCORES['resistance possible']
        self.level_descriptions = {
            score: description
            for description, score in PHENOTYPE_SCORES.items()}

    def check(self, worksheet, rule_sets, col, section_entries):
        for rule_set in rule_sets:
            if rule_set.first_column <= col <= rule_set.last_column:
                lower_entry = section_entries[0]
                match = re.match(self.lower_pattern, lower_entry)
                lower = int(match.group(1))
                upper_entry = section_entries[2]
                match = re.match(self.upper_pattern, upper_entry)
                if match is None:
                    message = "Invalid upper fold shift of {!r} for {} in {}.".format(
                        upper_entry,
                        rule_set.drug_name,
                        worksheet.title)
                    raise ValueError(message)
                upper = int(match.group(1))
                for mutation, score in rule_set.mutations.items():
                    fold_shift_text = rule_set.fold_shifts.pop(mutation)
                    self.check_text(fold_shift_text,
                                    score,
                                    mutation,
                                    lower,
                                    upper,
                                    worksheet)
                for mutation, fold_shift_text in rule_set.fold_shifts.items():
                    self.check_text(fold_shift_text,
                                    0,
                                    mutation,
                                    lower,
                                    upper,
                                    worksheet)
                rule_set.fold_shifts.clear()
                break

    def check_text(self, fold_shift_text, score, mutation, lower, upper, worksheet):
        match = re.match(r'(\d+)x$', fold_shift_text)
        if match is None:
            message = "Invalid fold shift of {!r} for {} in {}.".format(
                fold_shift_text,
                mutation,
                worksheet.title)
            raise ValueError(message)
        fold_shift = int(match.group(1))
        if fold_shift < lower:
            expected_score = self.lower_score
        elif fold_shift > upper:
            expected_score = self.upper_score
        else:
            expected_score = self.middle_score
        if score != expected_score:
            expected_description = self.level_descriptions[
                expected_score]
            description = self.level_descriptions[score]
            message = ('Expected phenotype {} for {} in {} but '
                       'found {}.').format(expected_description,
                                           mutation,
                                           worksheet.title,
                                           description)
            raise ValueError(message)

    @staticmethod
    def check_missing(worksheet, rule_sets):
        for rule_set in rule_sets:
            if rule_set.fold_shifts:
                message = (
                    'Fold-shift levels not found in {} footer on {}'.format(
                        rule_set.drug_name,
                        worksheet.title))
                raise ValueError(message)


def dump_comments(args):
    worksheets = [ws
                  for ws in args.spreadsheet
                  if ws.title in READY_TABS]
    reader = WorksheetReader(worksheets)

    original_comments = defaultdict(set)  # {lower: {comment}}
    comment_counts = Counter()
    for entry in reader:
        if not entry.comments:
            continue
        comment = re.sub(r'[, ]*[\d.]+%', ' #%', entry.comments.strip(', '))
        comment = re.sub(r' *\([\d, ]+\)', ' (#)', comment)
        comment = re.sub(r'mean of \d+', 'mean of #', comment)
        comment = re.sub(r'[\d.]+x', '#x', comment)
        comment = re.sub(r'(Cliincal|Clincail|Clincal|Clincial|Clincla|'
                         r'clinial|Clinica|Clinicla|Clinucak|Clnical) ',
                         'Clinical ',
                         comment,
                         flags=re.IGNORECASE)
        lower_comment = comment.lower()
        comment_counts[lower_comment] += 1
        original_comments[lower_comment].add(comment)
    for is_clinical_printing in (True, False):
        print('=== Clinical Evidence ==='
              if is_clinical_printing
              else '=== No Clinical Evidence ===')
        for comment, count in sorted(comment_counts.items()):
            is_clinical = ('clinical' in comment and
                           'p/r' not in comment and
                           'no specific clinical' not in comment)
            if is_clinical == is_clinical_printing:
                print('{:-3d} {}'.format(count, min(original_comments[comment])))


def main():
    # TODO: Migrate fold-shift checks to new writer.
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


def read_rule_sets(ws, references, check_phenotypes=False):
    header_rows = []
    monitored_drugs = set()
    fold_shift_checker = FoldShiftChecker()
    rule_sets = None
    previous_rows = []
    for row in ws.rows:
        if rule_sets is None:
            if row[0].value and row[0].value.startswith('WT'):
                rule_sets = find_rule_sets(ws, header_rows, check_phenotypes)
                find_phenotype_scores(row, rule_sets)
            else:
                header_rows.append(row)
        elif row[0].value is not None:
            if not row[0].font.strike:
                find_phenotype_scores(row, rule_sets)
        else:
            # In footer rows.
            if previous_rows:
                previous_row = previous_rows[-1]
                for col, cell in enumerate(previous_row, 1):
                    label = str(cell.value).lower()
                    if label.startswith('positions monitored'):
                        drug_name = monitor_positions(ws,
                                                      rule_sets,
                                                      references,
                                                      col,
                                                      row[col-1].value)
                        monitored_drugs.add(drug_name)
            if len(previous_rows) >= 3 and check_phenotypes:
                label_row = previous_rows[-3]
                for col, cell in enumerate(label_row, 1):
                    label = str(cell.value).lower()
                    if label in ('in vitro drug susceptibility:',
                                 'in virto drug susecptibility:'):
                        section_rows = (previous_rows[-2],
                                        previous_rows[-1],
                                        row)
                        section_entries = [section_row[col-1].value
                                           for section_row in section_rows]
                        fold_shift_checker.check(ws,
                                                 rule_sets,
                                                 col,
                                                 section_entries)
            previous_rows.append(row)
            if len(previous_rows) > 3:
                previous_rows.pop(0)

    if rule_sets is None:
        raise ValueError('No mutation started with WT in {}.'.format(ws.title))
    all_drugs = {rule_set.drug_name for rule_set in rule_sets}
    missing_drugs = all_drugs - monitored_drugs
    if missing_drugs:
        raise ValueError("No 'monitored positions' label for {} in {}.".format(
            ', '.join(sorted(missing_drugs)),
            ws.title))
    if check_phenotypes:
        fold_shift_checker.check_missing(ws, rule_sets)
    return rule_sets


def monitor_positions(worksheet, rule_sets, references, col, positions):
    rule_set = None
    for rule_set in rule_sets:
        if rule_set.first_column <= col <= rule_set.last_column:
            if positions is None:
                raise ValueError(
                    "No list of monitored positions for {} in {}.".format(
                        rule_set.drug_name,
                        worksheet.title))
            if str(positions).lower() == 'none':
                break
            reference = references[(rule_set.genotype.upper(), rule_set.region)]
            try:
                seen_variants = {int(pos): set()
                                 for pos in str(positions).split(',')}
            except ValueError:
                raise ValueError(
                    'Invalid monitored position for {} in {}: {!r}.'.format(
                        rule_set.drug_name,
                        worksheet.title,
                        positions))
            for mutation in rule_set.mutations:
                if mutation is None or '+' in mutation or ' ' in mutation:
                    continue
                match = re.match(r'^([A-Z])(\d+)([A-Zid]+)$', mutation)
                assert match, mutation
                wild_type = match.group(1)
                pos = int(match.group(2))
                expected_wild_type = reference.sequence[pos - 1]
                assert wild_type == expected_wild_type, (reference.name,
                                                         mutation,
                                                         expected_wild_type)
                variants = match.group(3)
                pos_variants = seen_variants.get(pos)
                if pos_variants is not None:
                    pos_variants.update(variants)
            for pos, pos_variants in seen_variants.items():
                try:
                    wild_type = reference.sequence[pos - 1]
                except IndexError:
                    ref_size = len(reference.sequence)
                    raise IndexError(
                        f'Position {pos} is beyond maximum of {ref_size} for '
                        f'{rule_set.drug_name}.')
                mutation = '{}{}!{}{}'.format(wild_type,
                                              pos,
                                              ''.join(sorted(pos_variants)),
                                              wild_type)
                rule_set.mutations[mutation] = 'Effect unknown'
            break
    return rule_set.drug_name


def find_phenotype_scores(row, rule_sets):
    mutation = row[0].value
    for rule_set in rule_sets:
        condition = None if mutation.startswith('WT') else mutation
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
                rule_set.mutations[condition] = score
        if rule_set.fold_shifts is not None:
            fold_shift_cell = row[rule_set.phenotype_column-3]
            fold_shift = fold_shift_cell.value
            if fold_shift is not None and fold_shift not in ('1x', '1'):
                rule_set.fold_shifts[condition] = fold_shift


def create_rule_set(sheet_name,
                    drug_name,
                    first_column,
                    phenotype_column,
                    last_column,
                    check_phenotypes=False):
    sheet_region, genotype = sheet_name.split('_GT')
    lower_region = sheet_region.lower()
    for known_region in HCV_REGIONS:
        if known_region.lower() == lower_region:
            region = known_region
            break
    else:
        region = sheet_region
    short_name = get_short_drug_name(drug_name)
    fold_shifts = {} if check_phenotypes else None
    return RuleSet(region,
                   genotype,
                   short_name,
                   first_column,
                   phenotype_column,
                   last_column,
                   {},
                   fold_shifts)


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


def find_rule_sets(ws, header_rows, check_phenotypes):
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
                                                             drug_range.max_col,
                                                             check_phenotypes))
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
            drug_code = drug_codes.get(drug_name)
            if drug_code is None:
                message = 'Unknown drug: {} found in {} {}.'.format(
                    drug_name,
                    rule_set.region,
                    rule_set.genotype)
                raise ValueError(message)
            drug_summary = drug_summaries[drug_name] = dict(name=drug_name,
                                                            code=drug_code,
                                                            genotypes=[])
        positions = defaultdict(dict)  # {pos: {score: MutationSet}}
        combinations = {}  # {mutation: score}
        for mutation, score in rule_set.mutations.items():
            if mutation is None:
                # noinspection PyTypeChecker
                positions[None][score] = 'TRUE'
                continue
            if '+' in mutation or ' ' in mutation:
                combinations[mutation] = score
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
        combination_changes = []
        for combination, combination_score in combinations.items():
            try:
                component_score = calculate_component_score(combination, positions)
            except ValueError:
                message = 'Components score could not be calculated for {} {} {}.'.format(
                    rule_set.region,
                    rule_set.genotype,
                    combination)
                raise ValueError(message)

            component_score = min(component_score, 8)
            if component_score != combination_score:
                combination_changes.append((combination, combination_score, component_score))
        if combination_changes:
            print('Combination changes for',
                  rule_set.region,
                  rule_set.genotype,
                  rule_set.drug_name)
        for combination, combination_score, component_score in combination_changes:
            print(' ', combination, combination_score, component_score)
        genotype = rule_set.genotype.upper()
        reference = references[(genotype, rule_set.region)]
        reference_name = reference.name
        score_formula = build_score_formula(positions)
        drug_summary['genotypes'].append(dict(genotype=genotype,
                                              region=rule_set.region,
                                              reference=reference_name,
                                              rules=score_formula))
    drugs = sorted(drug_summaries.values(), key=itemgetter('code'))
    yaml.dump(drugs, rules_file, default_flow_style=False, Dumper=SplitDumper)


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
        score_terms.insert(0, (condition, score))
    if not score_terms:
        score_terms.append(('TRUE', '0'))
    score_formula = 'SCORE FROM ( {} )'.format(', '.join(
        '{} => {}'.format(mutation_set, score)
        for mutation_set, score in score_terms))
    return score_formula


if __name__ == '__main__':
    main()
