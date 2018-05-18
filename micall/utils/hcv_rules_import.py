import re
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter, FileType, Namespace
from collections import namedtuple, defaultdict, Counter
from copy import copy
from functools import partial
from itertools import product
from operator import itemgetter
import os

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
STRIKE = '**strike**'
DATA_CHANGES = {
    'NS3_GT1a': {
        # 'A1': ('OLD', 'NEW'),
        'A131': ('S343N', 'T343N'),  # Wild type mismatch
        'A222': ('R155K_D168A', 'R155K+D168A'),
        'A149': ('Q9Any+A156T', 'Q9!Q+A156T'),
        'A259': ('E357K+D168Any', 'E357K+D168!D')
    },
    'NS3_GT1b': {
        'A40': ('V71I', 'I71I')
    },
    'NS3_GT2': {
        'Y6': ('Aunaprevir in GT2', 'Asunaprevir in GT2'),
        'E120': (None, 'positions monitored'),
        'E121': (None, 'None'),  # Replace blank with the word None so I know it's not a mistake.
        'J120': (None, 'positions monitored'),
        'J121': (None, 'None'),
        'AA120': (None, 'positions monitored'),
        'AA121': (None, 'None'),
    },
    'NS3_GT3': {
        # 'J26': # TODO: How to handle resistant WT with other mutations? Ignored for now.
        'I131': (None, 'None'),
        'AB131': (None, 'None'),
        'AI130': (None, 'positions monitored'),
        'AI131': (None, 'None'),
        'A76': ('Q186L', 'D186L'),
        'A77': ('Q186R', 'D186R'),
        'A98': ('Y56H+D168V', 'Y56H+Q168V'),  # TODO: Any chance there's confusion between D186 and Q168?
        'A85': ('Q41R+V551', 'Q41R+V551!V'),
        'A99': ('Q80R+S166T', 'Q80R+A166T'),
    },
    'NS3_GT4': {
        'T7': ('Gelcaprevir in GT4', 'Glecaprevir in GT4'),
        'U106': (156168, '156, 168'),  # Put a space after comma so Excel knows it's text.
        'N111': (156168, '156, 168'),
        'AG106': (None, 'None'),
    },
    'NS3_GT5': {
        'G6': ('Paritapavir in GT5', 'Paritaprevir in GT5'),
        'D85': (None, 'positions monitored'),
        'D86': (None, 'None'),
        'I85': (None, 'positions monitored'),
        'I86': (None, 'None'),
        'Y85': (None, 'positions monitored'),
        'Y86': (None, 'None'),
        'S86': (156168, '156, 168'),
    },
    'NS3_GT6': {
        'J106': (None, 'positions monitored'),
        'J107': (None, 'None'),
        'Z106': (None, 'positions monitored'),
        'Z107': (None, 'None'),
        'AE106': (None, 'positions monitored'),
        'AE107': (None, 'None'),
        'N113': (41156168, '41, 156, 168'),
        'T108': (80156168, '80, 156, 168'),
    },
    'NS5A_GT1a': {
        'O30': (17, 'likely susceptible'),
        'O31': (17, 'likely susceptible'),
        'O32': (17, 'likely susceptible'),
        'A45': ('Q30del', 'Q30d'),
        'A119': ('T213A', 'A213A'),
        'A126': ('A377T', 'T377T'),
        'A145': ('M28A+M28T', STRIKE),  # TODO: What does it mean to have two mutations at same position?
        'A230': ('L24M/T+L31IV+Y93H', 'K24MT+L31IV+Y93H'),
        'T19': ('Velpastasvir in GT1a', 'Velpatasvir in GT1a'),
        'AX19': ('Odalasvir in GT1a', None),  # TODO: delete it instead of hiding it?
    },
    'NS5A_GT1b': {
        'Q15': ('Velpastasvir in GT1b', 'Velpatasvir in GT1b'),
        'A26': ('P29del', 'P29d'),
        'A50': ('P32del', 'P32d'),
        'A90': ('L28M+P32del', 'L28M+P32d'),
        'A107': ('L31F+P32del', 'L31F+P32d'),
        'A111': ('L31I/V+Y93H', 'L31IV+Y93H'),
        'A121': ('L31V+P32del', 'L31V+P32d'),
        'A185': ('L31I/M+Q54H+Q62E+Y93H', 'L31IM+Q54H+Q62E+Y93H'),
    },
    'NS5A_GT2': {
        'M9': ('Velpastasvir in GT2', 'Velpatasvir in GT2'),
        'A11': (' WT (2a_JFH)_L31', 'WT (2a_JFH)_L31'),
        'A26': ('L28F ', 'F28F'),  # TODO: Is position 28 for different subtypes?
        'A27': ('L28L', 'F28L'),
        'A42': ('M31I (GT2a_J6)', STRIKE),  # TODO: Are we using subtypes now?
        'A43': ('M31L (GT2a_J6)', STRIKE),
        'A47': ('L31M (GT2b_MD2b)', STRIKE),
        'A51': ('P58A (GT2b)', STRIKE),
        'A61': ('C92S (GT2a)', STRIKE),
        'A62': ('C92S (GT2b)', STRIKE),
        'A63': ('C92T (GT2a)', STRIKE),
        'A64': ('C92T (GT2b)', STRIKE),
        'A68': ('Y93F (GT2a)', STRIKE),
        'A69': ('Y93F (GT2b)', STRIKE),
        'A98': ('M31V+Y93H', 'L31V+Y93H'),
        'A102': ('F28+K30R+L31M', 'F28!F+K30R+L31M'),
        'T133': (None, 'None'),
        'AE139': (None, 'None'),
        'I139': (None, 'positions monitored'),
        'I140': (None, 'None'),
    },
    'NS5A_GT3': {
        'O14': ('Velpastasvir in GT3', 'Velpatasvir in GT3'),
        'AA14': ('Pibrentasivr in GT3', 'Pibrentasvir in GT3'),
        'A53': ('Q41K', 'K41K'),
        'A95': ('M28V_Y93H', 'M28V+Y93H'),
        'A97': ('A30E_S62T', 'A30E+S62T'),
        'A125': ('S62L_Y93H', 'S62L+Y93H'),
        'A127': ('S62T_Y93H', 'S62T+Y93H'),
        'A128': ('E92A_Y93H', 'E92A+Y93H'),
        'A133': ('A30K_P58L_Y93H', 'A30K+P58L+Y93H'),
        'W160': (None, 'None'),
        'K166': (None, 'positions monitored'),
        'K167': (None, 'None'),
    },
    'NS5A_GT4': {
        'O10': ('Velpastasvir in GT4', 'Velpatasvir in GT4'),
        'A20': ('L28V (GT4a)', STRIKE),
        'A21': ('L28V (GT4d)', STRIKE),
        'A72': ('L28V+T58S', 'L28V+P58S'),
        'A73': ('L28V+T58S', 'L28V+P58S'),
        'D117': (None, 'None'),
        'J124': (None, 'None'),
    },
    'NS5A_GT5': {
        'L7': ('Velpastasvir in GT5', 'Velpatasvir in GT5'),
        'C70': (None, 'None'),
        'S70': (None, 'None'),
        'AD70': (None, 'None'),
        'H77': (None, 'None'),
    },
    'NS5A_GT6': {
        'N9': ('Velpastasvir in GT6', 'Velpatasvir in GT6'),
        'D82': (None, 'None'),
        'U82': (None, 'None'),
        'AF82': (None, 'None'),
        'J89': (None, 'None'),
    }}
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
            section = Namespace(drug_name=drug_name,
                                sheet_name=ws.title)
            for footer_reader in self.footer_readers:
                footer_coordinates = product(
                    range(footer_row, ws.max_row + 1),
                    range(cell_range.min_col, cell_range.max_col + 1))
                footer_reader.read(ws, footer_coordinates, section)
            for row_num in range(wild_type_row, footer_row):
                if ws.cell(row_num, 1).font.strike:
                    continue
                fields = {heading: ws.cell(row_num, col).value
                          for col, heading in column_headings.items()}
                entry = Namespace(section=section, **fields)
                yield entry

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
        elif monitored_positions.lower() == 'none':
            monitored_positions = []
        else:
            monitored_positions = list(map(
                int,
                str(monitored_positions).split(', ')))
        section.monitored_positions = monitored_positions

    def write_errors(self, report):
        if self.missing:
            report.write('No monitored positions for {}.\n'.format(
                ', '.join(self.missing)))


class FoldRangesReader:
    def __init__(self):
        self.invalid = []
        self.lower_pattern = r'<(\d+)x? FS, likely susceptible$'
        self.lower_score = PHENOTYPE_SCORES['likely susceptible']
        self.upper_pattern = r'>(\d+)x? FS, resistance likely$'
        self.upper_score = PHENOTYPE_SCORES['resistance likely']
        self.middle_score = PHENOTYPE_SCORES['resistance possible']
        self.level_descriptions = {
            score: description
            for description, score in PHENOTYPE_SCORES.items()}

    def read(self, ws, footer_coordinates, section):
        for row_num, col_num in footer_coordinates:
            label = ws.cell(row_num, col_num).value
            if label and label.lower() in ('in vitro drug susceptibility:',
                                           'in virto drug susecptibility:'):
                lower_text = ws.cell(row_num + 1, col_num).value
                upper_text = ws.cell(row_num + 3, col_num).value
                lower_match = re.match(self.lower_pattern, lower_text)
                upper_match = re.match(self.upper_pattern, upper_text)
                if upper_match and lower_match:
                    section.lower_fold = int(lower_match.group(1))
                    section.upper_fold = int(upper_match.group(1))
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
                break

    def write_errors(self, report):
        for message in self.invalid:
            print(message, file=report)


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
        rule_sets = read_rule_sets(ws, references, check_phenotypes=True)
        all_rule_sets.extend(rule_sets)
    write_rules(all_rule_sets, references, args.rules_yaml)


def make_data_changes(wb):
    for sheet_name, changes in DATA_CHANGES.items():
        ws = wb[sheet_name]
        for coordinate, (old, new) in changes.items():
            cell = ws[coordinate]
            assert cell.value == old, (sheet_name, coordinate, cell.value)
            if new == STRIKE:
                new_font = copy(cell.font)
                new_font.strike = True
                cell.font = new_font
            else:
                ws[coordinate] = new


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
