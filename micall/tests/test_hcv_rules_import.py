from argparse import Namespace
from copy import copy
from unittest import TestCase, skip

from io import StringIO

from openpyxl import Workbook

from micall.utils.hcv_rules_import import create_rule_set, RuleSet, \
    write_rules, read_rule_sets, load_references, WorksheetReader, \
    MonitoredPositionsReader, FoldRangesReader, RulesWriter

REFERENCES = load_references()


class CreateRuleSetTest(TestCase):
    def test(self):
        sheet_name = 'R1_GT2x'
        drug_name = 'Paritaprevir in GT1a'
        first_column = 21
        phenotype_column = 23
        last_column = 25
        expected_region = 'R1'
        expected_genotype = '2x'
        expected_drug_name = 'Paritaprevir'

        rule_set = create_rule_set(sheet_name,
                                   drug_name,
                                   first_column,
                                   phenotype_column,
                                   last_column)

        self.assertEqual(first_column, rule_set.first_column)
        self.assertEqual(phenotype_column, rule_set.phenotype_column)
        self.assertEqual(last_column, rule_set.last_column)
        self.assertEqual(expected_region, rule_set.region)
        self.assertEqual(expected_genotype, rule_set.genotype)
        self.assertEqual(expected_drug_name, rule_set.drug_name)


def create_worksheet(title, row_data):
    """ Create a worksheet from a list of lists of cell values.

    Use these special cell values:
    '' - empty cell
    '<' - merge cell with left neighbour
    :param str title: the title for the worksheet
    :param list row_data: a list of lists of cell values
    """
    wb = Workbook()
    ws = wb.create_sheet(title)
    for i, row in enumerate(row_data, 1):
        row = [None if cell == '' else cell for cell in row]
        ws.append(row)
        start_column = None
        for j, cell in enumerate(row + [None], 1):
            if cell == '<':
                if start_column is None:
                    start_column = j - 1
            elif start_column is not None:
                end_column = j - 1
                ws.merge_cells(start_row=i,
                               end_row=i,
                               start_column=start_column,
                               end_column=end_column)
                start_column = None
    return ws


class WorksheetReaderTest(TestCase):
    def setUp(self):
        self.expected_errors = ''

    def assertReads(self, expected_entries, worksheets, *footer_readers):
        errors = StringIO()
        reader = WorksheetReader(worksheets, *footer_readers)

        entries = list(reader)
        reader.write_errors(errors)

        self.assertEqual(self.expected_errors, errors.getvalue())
        self.assertEqual(expected_entries, entries)

    def test_one_sheet(self):
        row_data = [
            ['',     'Ignored header row'],
            ['',     'Example1 drug', '<',                  '<', 'Example2', '<', '<'],
            ['',     'a',             'Phenotype',          'b', 'a',        'b', 'Phenotype'],
            ['WT',   '',              'likely susceptible', '',  '',         '',  'likely susceptible'],
            ['V36A', '',              'resistance likely',  '',  '',         '',  'likely susceptible'],
            ['T40A', '',              'likely susceptible', '',  '',         '',  'resistance possible']]
        worksheets = [create_worksheet('NS3_GT1a', row_data)]
        expected_section1 = Namespace(drug_name='Example1',
                                      sheet_name='NS3_GT1a')
        expected_section2 = Namespace(drug_name='Example2',
                                      sheet_name='NS3_GT1a')
        expected_entries = [Namespace(mutation='WT',
                                      section=expected_section1,
                                      phenotype='likely susceptible',
                                      a=None,
                                      b=None),
                            Namespace(mutation='V36A',
                                      section=expected_section1,
                                      phenotype='resistance likely',
                                      a=None,
                                      b=None),
                            Namespace(mutation='T40A',
                                      section=expected_section1,
                                      phenotype='likely susceptible',
                                      a=None,
                                      b=None),
                            Namespace(mutation='WT',
                                      section=expected_section2,
                                      phenotype='likely susceptible',
                                      a=None,
                                      b=None),
                            Namespace(mutation='V36A',
                                      section=expected_section2,
                                      phenotype='likely susceptible',
                                      a=None,
                                      b=None),
                            Namespace(mutation='T40A',
                                      section=expected_section2,
                                      phenotype='resistance possible',
                                      a=None,
                                      b=None)]

        self.assertReads(expected_entries, worksheets)

    def test_two_sheets(self):
        row_data1 = [
            ['',     'Ignored header row'],
            ['',     'Example1 drug', '<',                  '<'],
            ['',     'a',             'Phenotype',          'b'],
            ['WT',   '',              'likely susceptible', ''],
            ['V36A', '',              'resistance likely',  ''],
            ['T40A', '',              'likely susceptible', '']]
        row_data2 = [
            ['',     'Example2', '<', '<'],
            ['',     'a',        'b', 'Phenotype'],
            ['WT',   '',         '',  'likely susceptible'],
            ['V36A', '',         '',  'likely susceptible'],
            ['T40A', '',         '',  'resistance possible']]
        worksheets = [create_worksheet('NS3_GT1a', row_data1),
                      create_worksheet('NS3_GT1b', row_data2)]
        expected_section1 = Namespace(drug_name='Example1',
                                      sheet_name='NS3_GT1a')
        expected_section2 = Namespace(drug_name='Example2',
                                      sheet_name='NS3_GT1b')
        expected_entries = [Namespace(mutation='WT',
                                      section=expected_section1,
                                      phenotype='likely susceptible',
                                      a=None,
                                      b=None),
                            Namespace(mutation='V36A',
                                      section=expected_section1,
                                      phenotype='resistance likely',
                                      a=None,
                                      b=None),
                            Namespace(mutation='T40A',
                                      section=expected_section1,
                                      phenotype='likely susceptible',
                                      a=None,
                                      b=None),
                            Namespace(mutation='WT',
                                      section=expected_section2,
                                      phenotype='likely susceptible',
                                      a=None,
                                      b=None),
                            Namespace(mutation='V36A',
                                      section=expected_section2,
                                      phenotype='likely susceptible',
                                      a=None,
                                      b=None),
                            Namespace(mutation='T40A',
                                      section=expected_section2,
                                      phenotype='resistance possible',
                                      a=None,
                                      b=None)]

        self.assertReads(expected_entries, worksheets)

    def test_missing_wild_type(self):
        row_data1 = [
            ['',     'Ignored header row'],
            ['',     'Example1 drug', '<',                  '<'],
            ['',     'a',             'Phenotype',          'b'],
            ['V36A', '',              'resistance likely',  ''],
            ['T40A', '',              'likely susceptible', '']]
        row_data2 = [
            ['',     'Example2', '<', '<'],
            ['',     'a',        'b', 'Phenotype'],
            ['WT',   '',         '',  'likely susceptible'],
            ['V36A', '',         '',  'likely susceptible'],
            ['T40A', '',         '',  'resistance possible']]
        worksheets = [create_worksheet('NS3_GT1a', row_data1),
                      create_worksheet('NS3_GT1b', row_data2)]
        expected_section2 = Namespace(drug_name='Example2',
                                      sheet_name='NS3_GT1b')
        expected_entries = [Namespace(mutation='WT',
                                      section=expected_section2,
                                      phenotype='likely susceptible',
                                      a=None,
                                      b=None),
                            Namespace(mutation='V36A',
                                      section=expected_section2,
                                      phenotype='likely susceptible',
                                      a=None,
                                      b=None),
                            Namespace(mutation='T40A',
                                      section=expected_section2,
                                      phenotype='resistance possible',
                                      a=None,
                                      b=None)]
        self.expected_errors = """\
No wild type found in NS3_GT1a.
"""

        self.assertReads(expected_entries, worksheets)

    def test_blank(self):
        row_data = [
            ['',     'Example', '<'],
            ['',     'Phenotype'],
            ['WT',   'likely susceptible'],
            ['V36A', ''],
            ['T40A', 'resistance possible']]
        worksheets = [create_worksheet('NS3_GT1a', row_data)]
        expected_section = Namespace(drug_name='Example',
                                     sheet_name='NS3_GT1a')
        expected_entries = [Namespace(mutation='WT',
                                      section=expected_section,
                                      phenotype='likely susceptible'),
                            Namespace(mutation='T40A',
                                      section=expected_section,
                                      phenotype='resistance possible')]

        self.assertReads(expected_entries, worksheets)

    def test_mutations_stripped(self):
        row_data = [
            ['',       'Example', '<'],
            ['',       'Phenotype'],
            ['WT',     'likely susceptible'],
            [' V36A ', 'likely susceptible']]
        worksheets = [create_worksheet('NS3_GT1a', row_data)]
        expected_section = Namespace(drug_name='Example',
                                     sheet_name='NS3_GT1a')
        expected_entries = [Namespace(mutation='WT',
                                      section=expected_section,
                                      phenotype='likely susceptible'),
                            Namespace(mutation='V36A',
                                      section=expected_section,
                                      phenotype='likely susceptible')]
        self.assertReads(expected_entries, worksheets)

    def test_mutation_deletion(self):
        row_data = [
            ['',       'Example', '<'],
            ['',       'Phenotype'],
            ['WT',     'likely susceptible'],
            ['V36del', 'likely susceptible']]
        worksheets = [create_worksheet('NS3_GT1a', row_data)]
        expected_section = Namespace(drug_name='Example',
                                     sheet_name='NS3_GT1a')
        expected_entries = [Namespace(mutation='WT',
                                      section=expected_section,
                                      phenotype='likely susceptible'),
                            Namespace(mutation='V36d',
                                      section=expected_section,
                                      phenotype='likely susceptible')]
        self.assertReads(expected_entries, worksheets)

    def test_other_merged_section(self):
        row_data = [
            ['',     'Example', '<', '<'],
            ['',     'a',        'b', 'Phenotype'],
            ['WT',   '',         '',  'likely susceptible'],
            ['V36A', '',         '',  'likely susceptible'],
            ['T40A', '',         '',  'resistance possible'],
            ['',     'Other merged section', '<']]
        worksheets = [create_worksheet('NS3_GT1b', row_data)]
        expected_section = Namespace(drug_name='Example',
                                     sheet_name='NS3_GT1b')
        expected_entries = [Namespace(mutation='WT',
                                      section=expected_section,
                                      phenotype='likely susceptible',
                                      a=None,
                                      b=None),
                            Namespace(mutation='V36A',
                                      section=expected_section,
                                      phenotype='likely susceptible',
                                      a=None,
                                      b=None),
                            Namespace(mutation='T40A',
                                      section=expected_section,
                                      phenotype='resistance possible',
                                      a=None,
                                      b=None)]

        self.assertReads(expected_entries, worksheets)

    def test_monitored_positions(self):
        row_data = [
            ['',     'Example', '<'],
            ['',     'Phenotype'],
            ['WT',   'likely susceptible'],
            ['V36A', 'likely susceptible'],
            ['T40A', 'resistance possible'],
            ['',     'Positions monitored:'],
            ['',     '36, 99']]
        worksheets = [create_worksheet('NS3_GT1b', row_data)]
        expected_section = Namespace(drug_name='Example',
                                     sheet_name='NS3_GT1b',
                                     monitored_positions=[36, 99])
        expected_entries = [Namespace(mutation='WT',
                                      section=expected_section,
                                      phenotype='likely susceptible'),
                            Namespace(mutation='V36A',
                                      section=expected_section,
                                      phenotype='likely susceptible'),
                            Namespace(mutation='T40A',
                                      section=expected_section,
                                      phenotype='resistance possible')]

        self.assertReads(expected_entries, worksheets, MonitoredPositionsReader())

    def test_monitored_position_int(self):
        row_data = [
            ['',     'Example', '<'],
            ['',     'Phenotype'],
            ['WT',   'likely susceptible'],
            ['V36A', 'likely susceptible'],
            ['T40A', 'resistance possible'],
            ['',     'Positions monitored:'],
            ['',     36]]
        worksheets = [create_worksheet('NS3_GT1b', row_data)]
        expected_section = Namespace(drug_name='Example',
                                     sheet_name='NS3_GT1b',
                                     monitored_positions=[36])
        expected_entries = [Namespace(mutation='WT',
                                      section=expected_section,
                                      phenotype='likely susceptible'),
                            Namespace(mutation='V36A',
                                      section=expected_section,
                                      phenotype='likely susceptible'),
                            Namespace(mutation='T40A',
                                      section=expected_section,
                                      phenotype='resistance possible')]

        self.assertReads(expected_entries, worksheets, MonitoredPositionsReader())

    def test_no_monitored_positions(self):
        row_data = [
            ['',     'Example', '<'],
            ['',     'Phenotype'],
            ['WT',   'likely susceptible'],
            ['V36A', 'likely susceptible'],
            ['T40A', 'resistance possible'],
            ['',     'Nothing monitored.']]
        worksheets = [create_worksheet('NS3_GT1b', row_data)]
        self.expected_errors = """\
No monitored positions for Example in NS3_GT1b.
"""
        expected_section = Namespace(drug_name='Example',
                                     sheet_name='NS3_GT1b',
                                     monitored_positions=[])
        expected_entries = [Namespace(mutation='WT',
                                      section=expected_section,
                                      phenotype='likely susceptible'),
                            Namespace(mutation='V36A',
                                      section=expected_section,
                                      phenotype='likely susceptible'),
                            Namespace(mutation='T40A',
                                      section=expected_section,
                                      phenotype='resistance possible')]

        self.assertReads(expected_entries, worksheets, MonitoredPositionsReader())

    def test_not_indicated(self):
        row_data = [
            ['',     'Example', '<'],
            ['',     'Phenotype'],
            ['WT',   'likely susceptible'],
            ['V36A', 'likely susceptible'],
            ['T40A', 'resistance possible'],
            ['',     'Nothing monitored.'],
            ['',     'Not indicated in Canada']]
        worksheets = [create_worksheet('NS3_GT1b', row_data)]
        expected_entries = []

        self.assertReads(expected_entries, worksheets, MonitoredPositionsReader())

    def test_fold_shift_ranges(self):
        row_data = [
            ['',     'Example', '<', '<'],
            ['',     'Phenotype'],
            ['WT',   'likely susceptible'],
            ['V36A', 'likely susceptible'],
            ['T40A', 'resistance possible'],
            ['',     'In vitro drug susceptibility:'],
            ['',     '<20x FS, likely susceptible'],
            ['',     '20-100x FS, resistance possible'],
            ['',     '>100x FS, resistance likely'],
            ['',     '',         'Positions monitored:'],
            ['',     '',         '36, 99']]
        worksheets = [create_worksheet('NS3_GT1b', row_data)]
        expected_section = Namespace(drug_name='Example',
                                     sheet_name='NS3_GT1b',
                                     monitored_positions=[36, 99],
                                     lower_fold=20,
                                     upper_fold=100)
        expected_entries = [Namespace(mutation='WT',
                                      section=expected_section,
                                      phenotype='likely susceptible'),
                            Namespace(mutation='V36A',
                                      section=expected_section,
                                      phenotype='likely susceptible'),
                            Namespace(mutation='T40A',
                                      section=expected_section,
                                      phenotype='resistance possible')]

        self.assertReads(expected_entries,
                         worksheets,
                         MonitoredPositionsReader(),
                         FoldRangesReader())

    def test_fold_shift_typos(self):
        row_data = [
            ['',     'Example', '<', '<'],
            ['',     'Phenotype'],
            ['WT',   'likely susceptible'],
            ['V36A', 'likely susceptible'],
            ['T40A', 'resistance possible'],
            ['',     'In virto Drug Susecptibility:'],
            ['',     '<20x FS, likely susceptible'],
            ['',     '20-100x FS, resistance possible'],
            ['',     '>100x FS, resistance likely']]
        worksheets = [create_worksheet('NS3_GT1b', row_data)]
        expected_section = Namespace(drug_name='Example',
                                     sheet_name='NS3_GT1b',
                                     lower_fold=20,
                                     upper_fold=100)
        expected_entries = [Namespace(mutation='WT',
                                      section=expected_section,
                                      phenotype='likely susceptible'),
                            Namespace(mutation='V36A',
                                      section=expected_section,
                                      phenotype='likely susceptible'),
                            Namespace(mutation='T40A',
                                      section=expected_section,
                                      phenotype='resistance possible')]

        self.assertReads(expected_entries, worksheets, FoldRangesReader())

    def test_fold_shift_absolute(self):
        row_data = [
            ['',     'Example', '<', '<'],
            ['',     'Phenotype'],
            ['WT',   'likely susceptible'],
            ['V36A', 'likely susceptible'],
            ['T40A', 'resistance possible'],
            ['',     'In vitro Drug Susceptibility:'],
            ['',     'Absolute FS values -'],
            ['',     '<20x FS, likely susceptible'],
            ['',     '20-100x FS, resistance possible'],
            ['',     '>100x FS, resistance likely']]
        worksheets = [create_worksheet('NS3_GT1b', row_data)]
        expected_section = Namespace(drug_name='Example',
                                     sheet_name='NS3_GT1b',
                                     lower_fold=20,
                                     upper_fold=100)
        expected_entries = [Namespace(mutation='WT',
                                      section=expected_section,
                                      phenotype='likely susceptible'),
                            Namespace(mutation='V36A',
                                      section=expected_section,
                                      phenotype='likely susceptible'),
                            Namespace(mutation='T40A',
                                      section=expected_section,
                                      phenotype='resistance possible')]

        self.assertReads(expected_entries, worksheets, FoldRangesReader())

    def test_invalid_fold_shift(self):
        row_data = [
            ['',     'Example', '<', '<'],
            ['',     'Phenotype'],
            ['WT',   'likely susceptible'],
            ['V36A', 'likely susceptible'],
            ['T40A', 'resistance possible'],
            ['',     'In vitro drug susceptibility:'],
            ['',     'less than 20x FS, moist'],
            ['',     '20-100x FS, resistance possible'],
            ['',     '>100x FS, resistance likely']]
        worksheets = [create_worksheet('NS3_GT1b', row_data)]
        self.expected_errors = """\
Invalid lower fold shift of 'less than 20x FS, moist' for Example in NS3_GT1b.
"""
        expected_section = Namespace(drug_name='Example',
                                     sheet_name='NS3_GT1b',
                                     lower_fold=None,
                                     upper_fold=None)
        expected_entries = [Namespace(mutation='WT',
                                      section=expected_section,
                                      phenotype='likely susceptible'),
                            Namespace(mutation='V36A',
                                      section=expected_section,
                                      phenotype='likely susceptible'),
                            Namespace(mutation='T40A',
                                      section=expected_section,
                                      phenotype='resistance possible')]

        self.assertReads(expected_entries, worksheets, FoldRangesReader())

    def test_heading_underscores(self):
        row_data = [
            ['',     'Example',    '<'],
            ['',     'Fold-Shift', 'Phenotype'],
            ['WT',   '1x',         'likely susceptible'],
            ['V36A', '2x',         'likely susceptible'],
            ['T40A', '50x',        'resistance possible']]
        worksheets = [create_worksheet('NS3_GT1b', row_data)]
        expected_section = Namespace(drug_name='Example',
                                     sheet_name='NS3_GT1b')
        expected_entries = [Namespace(mutation='WT',
                                      section=expected_section,
                                      phenotype='likely susceptible',
                                      fold_shift='1x'),
                            Namespace(mutation='V36A',
                                      section=expected_section,
                                      phenotype='likely susceptible',
                                      fold_shift='2x'),
                            Namespace(mutation='T40A',
                                      section=expected_section,
                                      phenotype='resistance possible',
                                      fold_shift='50x')]

        self.assertReads(expected_entries, worksheets)

    def test_heading_stripped(self):
        row_data = [
            ['',     'Example',             '<'],
            ['',     'Phenotype',           'Comment*'],
            ['WT',   'likely susceptible'],
            ['V36A', 'likely susceptible'],
            ['T40A', 'resistance possible', 'unreliable']]
        worksheets = [create_worksheet('NS3_GT1b', row_data)]
        expected_section = Namespace(drug_name='Example',
                                     sheet_name='NS3_GT1b')
        expected_entries = [Namespace(mutation='WT',
                                      section=expected_section,
                                      phenotype='likely susceptible',
                                      comment=None),
                            Namespace(mutation='V36A',
                                      section=expected_section,
                                      phenotype='likely susceptible',
                                      comment=None),
                            Namespace(mutation='T40A',
                                      section=expected_section,
                                      phenotype='resistance possible',
                                      comment='unreliable')]
        reader = WorksheetReader(worksheets)

        entries = list(reader)

        self.assertEqual(expected_entries, entries)

    def test_strike_through(self):
        row_data = [
            ['',     'Example',             '<'],
            ['',     'Phenotype',           'Comment'],
            ['WT',   'likely susceptible'],
            ['V36A', 'likely susceptible', 'strike through to ignore'],
            ['T40A', 'resistance possible']]
        worksheets = [create_worksheet('NS3_GT1b', row_data)]
        cell = worksheets[0].cell(4, 1)
        font = copy(cell.font)
        font.strike = True
        cell.font = font
        expected_section = Namespace(drug_name='Example',
                                     sheet_name='NS3_GT1b')
        expected_entries = [Namespace(mutation='WT',
                                      section=expected_section,
                                      phenotype='likely susceptible',
                                      comment=None),
                            Namespace(mutation='T40A',
                                      section=expected_section,
                                      phenotype='resistance possible',
                                      comment=None)]

        self.assertReads(expected_entries, worksheets)


class ReadRuleSetTest(TestCase):
    def test_numeric_scores(self):
        row_data = [
            ['',     'Ignored header row'],
            ['',     'Example1 drug', '<',                  '<', 'Example2', '<', '<'],
            ['',     'a',             'Phenotype',          'b', 'a',        'b', 'Phenotype'],
            ['WT',   '',              'likely susceptible', '',  '',         '',  'likely susceptible'],
            ['V36A', '',              'resistance likely',  '',  '',         '',  'likely susceptible'],
            ['T40A', '',              'likely susceptible', '',  '',         '',  'resistance possible'],
            ['',     'Ignored footer row'],
            ['',     'Positions monitored...'],
            ['',     'None',          '',                   '',  '',         'Positions monitored:'],
            ['',     '',              '',                   '',  '',         'None']]
        expected_rule_sets = [RuleSet('NS3',
                                      '1a',
                                      'Example1',
                                      2,
                                      3,
                                      4,
                                      {'V36A': 8}),
                              RuleSet('NS3',
                                      '1a',
                                      'Example2',
                                      5,
                                      7,
                                      7,
                                      {'T40A': 4})]

        ws = create_worksheet('NS3_GT1a', row_data)

        rule_sets = read_rule_sets(ws, REFERENCES)

        self.assertEqual(expected_rule_sets, rule_sets)

    def test_check_phenotypes(self):
        row_data = [
            ['',     'Ignored header row'],
            ['',     'Example1 drug', '<', '<',               '<'],
            ['',     'Fold-shift', 'a', 'Phenotype',          'b'],
            ['WT',   '1x',         '',  'likely susceptible', ''],
            ['V36A', '200x',       '',  'resistance likely',  ''],
            ['T40A', '2x',         '',  'likely susceptible', ''],
            ['',     'Ignored footer row'],
            ['',     'In vitro drug susceptibility:'],
            ['',     '<20x FS, likely susceptible'],
            ['',     '20-100x FS, resistance possible'],
            ['',     '>100x FS, resistance likely'],
            ['',     'Positions monitored...'],
            ['',     'None']]
        expected_rule_sets = [RuleSet('NS3',
                                      '1a',
                                      'Example1',
                                      2,
                                      4,
                                      5,
                                      {'V36A': 8},
                                      {})]

        ws = create_worksheet('NS3_GT1a', row_data)

        rule_sets = read_rule_sets(ws, REFERENCES, check_phenotypes=True)

        self.assertEqual(expected_rule_sets, rule_sets)

    def test_phenotype_typos(self):
        row_data = [
            ['',     'Ignored header row'],
            ['',     'Example1 drug', '<', '<',               '<'],
            ['',     'Fold-shift', 'a', 'Phenotype',          'b'],
            ['WT',   '1x',         '',  'likely susceptible', ''],
            ['V36A', '200x',       '',  'resistance likely',  ''],
            ['T40A', '2x',         '',  'likely susceptible', ''],
            ['',     'Ignored footer row'],
            ['',     'In virto Drug Susecptibility:'],
            ['',     '<20x FS, likely susceptible'],
            ['',     '20-100x FS, resistance possible'],
            ['',     '>100x FS, resistance likely'],
            ['',     'Positions monitored...'],
            ['',     'None']]
        expected_rule_sets = [RuleSet('NS3',
                                      '1a',
                                      'Example1',
                                      2,
                                      4,
                                      5,
                                      {'V36A': 8},
                                      {})]

        ws = create_worksheet('NS3_GT1a', row_data)

        rule_sets = read_rule_sets(ws, REFERENCES, check_phenotypes=True)

        self.assertEqual(expected_rule_sets, rule_sets)

    def test_bad_phenotype(self):
        row_data = [
            ['',     'Ignored header row'],
            ['',     'Example1 drug', '<', '<',                '<'],
            ['',     'Fold-shift', 'a', 'Phenotype',           'b'],
            ['WT',   '1x',         '',  'likely susceptible',  ''],
            ['V36A', '200x',       '',  'resistance possible', ''],
            ['T40A', '2x',         '',  'likely susceptible',  ''],
            ['',     'Ignored footer row'],
            ['',     'In vitro drug susceptibility:'],
            ['',     '<20x FS, likely susceptible'],
            ['',     '20-100x FS, resistance possible'],
            ['',     '>100x FS, resistance likely'],
            ['',     'Positions monitored...'],
            ['',     'None']]

        ws = create_worksheet('NS3_GT1a', row_data)

        with self.assertRaisesRegex(ValueError,
                                    r'Expected phenotype resistance likely for '
                                    r'V36A in NS3_GT1a but found resistance '
                                    r'possible'):
            read_rule_sets(ws, REFERENCES, check_phenotypes=True)

    def test_invalid_fold_shift(self):
        row_data = [
            ['',     'Ignored header row'],
            ['',     'Example1 drug', '<', '<',                 '<'],
            ['',     'Fold-shift',  'a', 'Phenotype',           'b'],
            ['WT',   '1x',          '',  'likely susceptible',  ''],
            ['V36A', 'three times', '',  'resistance possible', ''],
            ['T40A', '2x',          '',  'likely susceptible',  ''],
            ['',     'Ignored footer row'],
            ['',     'In vitro drug susceptibility:'],
            ['',     '<20x FS, likely susceptible'],
            ['',     '20-100x FS, resistance possible'],
            ['',     '>100x FS, resistance likely'],
            ['',     'Positions monitored...'],
            ['',     'None']]

        ws = create_worksheet('NS3_GT1a', row_data)

        with self.assertRaisesRegex(ValueError,
                                    r"Invalid fold shift of 'three times' for "
                                    r"V36A in NS3_GT1a\."):
            read_rule_sets(ws, REFERENCES, check_phenotypes=True)

    def test_invalid_upper_fold_shift(self):
        row_data = [
            ['',     'Ignored header row'],
            ['',     'Example1 drug', '<', '<',                 '<'],
            ['',     'Fold-shift',  'a', 'Phenotype',           'b'],
            ['WT',   '1x',          '',  'likely susceptible',  ''],
            ['V36A', '200x',        '',  'resistance possible', ''],
            ['T40A', '2x',          '',  'likely susceptible',  ''],
            ['',     'Ignored footer row'],
            ['',     'In vitro drug susceptibility:'],
            ['',     '<20x FS, likely susceptible'],
            ['',     '20-100x FS, resistance possible'],
            ['',     'over 100, bad'],
            ['',     'Positions monitored...'],
            ['',     'None']]

        ws = create_worksheet('NS3_GT1a', row_data)

        with self.assertRaisesRegex(
                ValueError,
                r"Invalid upper fold shift of 'over 100, bad' for Example1 "
                r"in NS3_GT1a\."):
            read_rule_sets(ws, REFERENCES, check_phenotypes=True)

    def test_positions_monitored(self):
        row_data = [
            ['',     'Ignored header row'],
            ['',     'Example1 drug', '<',                  '<', 'Example2', '<', '<'],
            ['',     'a',             'Phenotype',          'b', 'a',        'b', 'Phenotype'],
            ['WT',   '',              'likely susceptible', '',  '',         '',  'likely susceptible'],
            ['V36A', '',              'resistance likely',  '',  '',         '',  'likely susceptible'],
            ['T40A', '',              'resistance likely', '',  '',         '',  'resistance possible'],
            ['',     'Ignored footer row'],
            ['',     'Positions monitored...'],
            ['',     '36,54, 62',     '',                   '',  '',         'Positions monitored:'],
            ['',     '',              '',                   '',  '',         'None']]
        expected_rule_sets = [RuleSet('NS3',
                                      '1a',
                                      'Example1',
                                      2,
                                      3,
                                      4,
                                      {'V36A': 8,
                                       'T40A': 8,
                                       'V36!AV': 'Effect unknown',
                                       'T54!T': 'Effect unknown',
                                       'R62!R': 'Effect unknown'}),
                              RuleSet('NS3',
                                      '1a',
                                      'Example2',
                                      5,
                                      7,
                                      7,
                                      {'T40A': 4})]

        ws = create_worksheet('NS3_GT1a', row_data)

        rule_sets = read_rule_sets(ws, REFERENCES)

        self.assertEqual(expected_rule_sets, rule_sets)

    def test_ns5a_positions_monitored(self):
        row_data = [
            ['',     'Ignored header row'],
            ['',     'Example1 drug', '<',                  '<', 'Example2', '<', '<'],
            ['',     'a',             'Phenotype',          'b', 'a',        'b', 'Phenotype'],
            ['WT',   '',              'likely susceptible', '',  '',         '',  'likely susceptible'],
            ['F36A', '',              'resistance likely',  '',  '',         '',  'likely susceptible'],
            ['Q40A', '',              'resistance likely', '',  '',         '',  'resistance possible'],
            ['',     'Ignored footer row'],
            ['',     'Positions monitored...'],
            ['',     '36,54, 62',     '',                   '',  '',         'Positions monitored:'],
            ['',     '',              '',                   '',  '',         'None']]
        expected_rule_sets = [RuleSet('NS5a',
                                      '1a',
                                      'Example1',
                                      2,
                                      3,
                                      4,
                                      {'F36A': 8,
                                       'Q40A': 8,
                                       'F36!AF': 'Effect unknown',
                                       'H54!H': 'Effect unknown',
                                       'E62!E': 'Effect unknown'}),
                              RuleSet('NS5a',
                                      '1a',
                                      'Example2',
                                      5,
                                      7,
                                      7,
                                      {'Q40A': 4})]

        ws = create_worksheet('NS5A_GT1a', row_data)

        rule_sets = read_rule_sets(ws, REFERENCES)

        self.assertEqual(expected_rule_sets, rule_sets)

    def test_combination(self):
        row_data = [
            ['',     'Ignored header row'],
            ['',     'Example1 drug', '<',                  '<', 'Example2', '<', '<'],
            ['',     'a',             'Phenotype',          'b', 'a',        'b', 'Phenotype'],
            ['WT',   '',              'likely susceptible', '',  '',         '',  'likely susceptible'],
            ['V36A', '',              'resistance possible', '',  '',        '',  'likely susceptible'],
            ['T40A', '',              'resistance possible', '',  '',        '',  'likely susceptible'],
            ['V36A+T40A', '',         'resistance likely',  '',  '',         '',  'likely susceptible'],
            ['',     'Positions monitored...', '',          '',  '',         'Positions monitored:'],
            ['',     '36',          '',                   '',  '',         'None']]
        expected_rule_sets = [RuleSet('NS3',
                                      '1a',
                                      'Example1',
                                      2,
                                      3,
                                      4,
                                      {'V36A': 4,
                                       'T40A': 4,
                                       'V36A+T40A': 8,
                                       'V36!AV': 'Effect unknown'}),
                              RuleSet('NS3',
                                      '1a',
                                      'Example2',
                                      5,
                                      7,
                                      7,
                                      {})]

        ws = create_worksheet('NS3_GT1a', row_data)

        rule_sets = read_rule_sets(ws, REFERENCES)

        self.assertEqual(expected_rule_sets, rule_sets)

    def test_no_positions_monitored_label(self):
        row_data = [
            ['',     'Ignored header row'],
            ['',     'Example1 drug', '<',                  '<', 'Example2', '<', '<'],
            ['',     'a',             'Phenotype',          'b', 'a',        'b', 'Phenotype'],
            ['WT',   '',              'likely susceptible', '',  '',         '',  'likely susceptible'],
            ['V36A', '',              'resistance likely',  '',  '',         '',  'likely susceptible'],
            ['T40A', '',              'likely susceptible', '',  '',         '',  'resistance possible'],
            ['',     'Ignored footer row'],
            ['',     'Nothing monitored...'],
            ['',     '36,54',         '',                   '',  '',         'Positions monitored:'],
            ['',     '',              '',                   '',  '',         'None']]

        ws = create_worksheet('NS3_GT1a', row_data)

        with self.assertRaisesRegex(
                ValueError,
                r"No 'monitored positions' label for Example1 in NS3_GT1a\."):
            read_rule_sets(ws, REFERENCES)

    def test_no_positions_monitored_list(self):
        row_data = [
            ['',     'Ignored header row'],
            ['',     'Example1 drug', '<',                  '<', 'Example2', '<', '<'],
            ['',     'a',             'Phenotype',          'b', 'a',        'b', 'Phenotype'],
            ['WT',   '',              'likely susceptible', '',  '',         '',  'likely susceptible'],
            ['V36A', '',              'resistance likely',  '',  '',         '',  'likely susceptible'],
            ['T40A', '',              'likely susceptible', '',  '',         '',  'resistance possible'],
            ['',     'Ignored footer row'],
            ['',     'Positions monitored...'],
            ['',     '',              '',                   '',  '',         'Positions monitored:'],
            ['',     '',              '',                   '',  '',         'None']]

        ws = create_worksheet('NS3_GT1a', row_data)

        with self.assertRaisesRegex(
                ValueError,
                r"No list of monitored positions for Example1 in NS3_GT1a\."):
            read_rule_sets(ws, REFERENCES)

    def test_bad_position(self):
        row_data = [
            ['',     'Ignored header row'],
            ['',     'Example1 drug', '<',                  '<', 'Example2', '<', '<'],
            ['',     'a',             'Phenotype',          'b', 'a',        'b', 'Phenotype'],
            ['WT',   '',              'likely susceptible', '',  '',         '',  'likely susceptible'],
            ['V36A', '',              'resistance likely',  '',  '',         '',  'likely susceptible'],
            ['T40A', '',              'likely susceptible', '',  '',         '',  'resistance possible'],
            ['',     'Ignored footer row'],
            ['',     'Positions monitored...'],
            ['',     '23,badint,29',              '',                   '',  '',         'Positions monitored:'],
            ['',     '',              '',                   '',  '',         'None']]

        ws = create_worksheet('NS3_GT1a', row_data)

        with self.assertRaisesRegex(
                ValueError,
                r"Invalid monitored position for Example1 in NS3_GT1a: '23,badint,29'\."):
            read_rule_sets(ws, REFERENCES)

    def test_no_wildtype(self):
        row_data = [
            ['',     'Ignored header row'],
            ['',     'Example1 drug', '<',                  '<', 'Example2', '<', '<'],
            ['',     'a',             'Phenotype',          'b', 'a',        'b', 'Phenotype'],
            ['XWT',  '',              'likely susceptible', '',  '',         '',  'likely susceptible'],
            ['V36A', '',              'resistance likely',  '',  '',         '',  'likely susceptible'],
            ['T40A', '',              'likely susceptible', '',  '',         '',  'resistance possible'],
            ['',     'Ignored footer row'],
            ['',     'Positions monitored...'],
            ['',     'None',          '',                   '',  '',         'Positions monitored:'],
            ['',     '',              '',                   '',  '',         'None']]

        ws = create_worksheet('NS3_GT1a', row_data)

        with self.assertRaisesRegex(
                ValueError,
                r"No mutation started with WT in NS3_GT1a\."):
            read_rule_sets(ws, REFERENCES)

    def test_wild_type_resistance(self):
        row_data = [
            ['',     'Ignored header row'],
            ['',     'Example1 drug', '<',                  '<', 'Example2', '<', '<'],
            ['',     'a',             'Phenotype',          'b', 'a',        'b', 'Phenotype'],
            ['WT',   '',              'resistance likely',  '',  '',         '',  'likely susceptible'],
            ['V36A', '',              '',                   '',  '',         '',  'likely susceptible'],
            ['T40A', '',              '',                   '',  '',         '',  'resistance possible'],
            ['',     'Ignored footer row'],
            ['',     'Positions monitored...'],
            ['',     'None',          '',                   '',  '',         'Positions monitored:'],
            ['',     '',              '',                   '',  '',         'None']]
        expected_rule_sets = [RuleSet('NS3',
                                      '1a',
                                      'Example1',
                                      2,
                                      3,
                                      4,
                                      {None: 8}),
                              RuleSet('NS3',
                                      '1a',
                                      'Example2',
                                      5,
                                      7,
                                      7,
                                      {'T40A': 4})]

        ws = create_worksheet('NS3_GT1a', row_data)

        rule_sets = read_rule_sets(ws, REFERENCES)

        self.assertEqual(expected_rule_sets, rule_sets)

    def test_wild_type_resistance_with_monitoring(self):
        row_data = [
            ['',     'Ignored header row'],
            ['',     'Example1 drug', '<',                  '<', 'Example2', '<', '<'],
            ['',     'a',             'Phenotype',          'b', 'a',        'b', 'Phenotype'],
            ['WT',   '',              'resistance likely',  '',  '',         '',  'likely susceptible'],
            ['V36A', '',              '',                   '',  '',         '',  'likely susceptible'],
            ['T40A', '',              '',                   '',  '',         '',  'resistance possible'],
            ['',     'Ignored footer row'],
            ['',     'Positions monitored...'],
            ['',     '36',          '',                   '',  '',         'Positions monitored:'],
            ['',     '',              '',                   '',  '',         'None']]
        expected_rule_sets = [RuleSet('NS3',
                                      '1a',
                                      'Example1',
                                      2,
                                      3,
                                      4,
                                      {None: 8,
                                       'V36!V': 'Effect unknown'}),
                              RuleSet('NS3',
                                      '1a',
                                      'Example2',
                                      5,
                                      7,
                                      7,
                                      {'T40A': 4})]

        ws = create_worksheet('NS3_GT1a', row_data)

        rule_sets = read_rule_sets(ws, REFERENCES)

        self.assertEqual(expected_rule_sets, rule_sets)

    def test_unknown_phenotype(self):
        row_data = [
            ['',     'Ignored header row'],
            ['',     'Example1 drug', '<',                  '<', 'Example2', '<', '<'],
            ['',     'a',             'Phenotype',          'b', 'a',        'b', 'Phenotype'],
            ['WT',   '',              'likely susceptible', '',  '',         '',  'likely susceptible'],
            ['V36A', '',              'trouble likely',  '',  '',         '',  'likely susceptible'],
            ['T40A', '',              'likely susceptible', '',  '',         '',  'resistance possible'],
            ['',     'Ignored footer row'],
            ['',     'Positions monitored...'],
            ['',     'None',          '',                   '',  '',         'Positions monitored:'],
            ['',     '',              '',                   '',  '',         'None']]

        ws = create_worksheet('NS3_GT1a', row_data)

        with self.assertRaisesRegex(
                ValueError,
                r"Unknown phenotype for Example1 V36A: 'trouble likely'\."):
            read_rule_sets(ws, REFERENCES)

    def test_missing_phenotype(self):
        row_data = [
            ['',     'Ignored header row'],
            ['',     'Example1 drug', '<',                  '<', 'Example2', '<', '<'],
            ['',     'a',             'x',                  'b', 'a',        'b', 'Phenotype'],
            ['WT',   '',              'likely susceptible', '',  '',         '',  'likely susceptible'],
            ['V36A', '',              'resistance likely',  '',  '',         '',  'likely susceptible'],
            ['T40A', '',              'likely susceptible', '',  '',         '',  'resistance possible'],
            ['',     'Ignored footer row'],
            ['',     'Positions monitored...'],
            ['',     'None',          '',                   '',  '',         'Positions monitored:'],
            ['',     '',              '',                   '',  '',         'None']]

        ws = create_worksheet('NS3_GT1a', row_data)

        with self.assertRaisesRegex(
                ValueError,
                r'No Phenotype column between columns 2 and 4 of NS3_GT1a\.'):
            read_rule_sets(ws, REFERENCES)


class RulesWriterTest(TestCase):
    def setUp(self):
        self.entries = []
        self.expected_rules = self.expected_errors = ''

    def write(self):
        rules = StringIO()
        errors = StringIO()
        writer = RulesWriter(rules, errors, REFERENCES)

        writer.write(self.entries)
        writer.write_errors()

        self.assertEqual(self.expected_rules, rules.getvalue())
        self.assertEqual(self.expected_errors, errors.getvalue())

    def assertWrites(self, expected_rules, entries):
        rules = StringIO()
        errors = StringIO()
        writer = RulesWriter(rules, errors, REFERENCES)

        writer.write(entries)
        writer.write_errors()

        self.assertEqual(self.expected_errors, errors.getvalue())
        self.assertEqual(expected_rules, rules.getvalue())

    def test_empty(self):
        """ Empty spreadsheet generates an empty list of rules. """
        entries = []
        expected_rules = """\
[]
"""

        self.assertWrites(expected_rules, entries)

    def test_simple(self):
        section = Namespace(drug_name='Paritaprevir', sheet_name='NS3_GT1a')
        entries = [Namespace(mutation='WT',
                             section=section,
                             phenotype='likely susceptible',
                             comment=None),
                   Namespace(mutation='T40A',
                             section=section,
                             phenotype='resistance possible',
                             comment=None)]
        expected_rules = """\
- code: PTV
  genotypes:
  - genotype: 1A
    reference: HCV1A-H77-NS3
    region: NS3
    rules: SCORE FROM ( T40A => 4 )
  name: Paritaprevir
"""

        self.assertWrites(expected_rules, entries)

    def test_no_resistance(self):
        section = Namespace(drug_name='Paritaprevir', sheet_name='NS3_GT1a')
        entries = [Namespace(mutation='WT',
                             section=section,
                             phenotype='likely susceptible',
                             comment=None)]
        expected_rules = """\
- code: PTV
  genotypes:
  - genotype: 1A
    reference: HCV1A-H77-NS3
    region: NS3
    rules: SCORE FROM ( TRUE => 0 )
  name: Paritaprevir
"""

        self.assertWrites(expected_rules, entries)

    def test_exclude_zeroes(self):
        section = Namespace(drug_name='Paritaprevir', sheet_name='NS3_GT1a')
        entries = [Namespace(mutation='WT',
                             section=section,
                             phenotype='likely susceptible',
                             comment=None),
                   Namespace(mutation='R20A',
                             section=section,
                             phenotype='likely susceptible',
                             comment=None),
                   Namespace(mutation='T40A',
                             section=section,
                             phenotype='resistance likely',
                             comment=None)]
        expected_rules = """\
- code: PTV
  genotypes:
  - genotype: 1A
    reference: HCV1A-H77-NS3
    region: NS3
    rules: SCORE FROM ( T40A => 8 )
  name: Paritaprevir
"""

        self.assertWrites(expected_rules, entries)

    def test_phenotype_case_insensitive(self):
        section = Namespace(drug_name='Paritaprevir', sheet_name='NS3_GT1a')
        entries = [Namespace(mutation='WT',
                             section=section,
                             phenotype='likely susceptible',
                             comment=None),
                   Namespace(mutation='T40A',
                             section=section,
                             phenotype='resiSTANce possible',
                             comment=None)]
        expected_rules = """\
- code: PTV
  genotypes:
  - genotype: 1A
    reference: HCV1A-H77-NS3
    region: NS3
    rules: SCORE FROM ( T40A => 4 )
  name: Paritaprevir
"""

        self.assertWrites(expected_rules, entries)

    def test_phenotype_invalid(self):
        section = Namespace(drug_name='Paritaprevir', sheet_name='NS3_GT1a')
        entries = [Namespace(mutation='WT',
                             section=section,
                             phenotype='likely susceptible',
                             comment=None),
                   Namespace(mutation='R20S',
                             section=section,
                             phenotype='bogus resistance',
                             comment=None),
                   Namespace(mutation='T40A',
                             section=section,
                             phenotype='resistance possible',
                             comment=None)]
        expected_rules = """\
- code: PTV
  genotypes:
  - genotype: 1A
    reference: HCV1A-H77-NS3
    region: NS3
    rules: SCORE FROM ( T40A => 4 )
  name: Paritaprevir
"""
        self.expected_errors = """\
Invalid phenotype: NS3_GT1a, Paritaprevir, R20S: bogus resistance.
"""

        self.assertWrites(expected_rules, entries)

    def test_unknown_drug(self):
        section1 = Namespace(drug_name='Paulrevir', sheet_name='NS3_GT1a')
        section2 = Namespace(drug_name='Paritaprevir', sheet_name='NS3_GT1a')
        entries = [Namespace(mutation='WT',
                             section=section1,
                             phenotype='likely susceptible',
                             comment=None),
                   Namespace(mutation='T41R',
                             section=section1,
                             phenotype='resistance likely',
                             comment=None),
                   Namespace(mutation='WT',
                             section=section2,
                             phenotype='likely susceptible',
                             comment=None),
                   Namespace(mutation='T40A',
                             section=section2,
                             phenotype='resistance possible',
                             comment=None)]
        expected_rules = """\
- code: PTV
  genotypes:
  - genotype: 1A
    reference: HCV1A-H77-NS3
    region: NS3
    rules: SCORE FROM ( T40A => 4 )
  name: Paritaprevir
"""
        self.expected_errors = """\
Unknown drug: NS3_GT1a: Paulrevir.
"""

        self.assertWrites(expected_rules, entries)

    def test_unknown_drugs(self):
        section1 = Namespace(drug_name='Paulrevir', sheet_name='NS3_GT1a')
        section2 = Namespace(drug_name='Paritaprevir', sheet_name='NS3_GT1a')
        section3 = Namespace(drug_name='Saynomorevir', sheet_name='NS3_GT1a')
        entries = [Namespace(mutation='WT',
                             section=section1,
                             phenotype='likely susceptible',
                             comment=None),
                   Namespace(mutation='T41R',
                             section=section1,
                             phenotype='resistance likely',
                             comment=None),
                   Namespace(mutation='WT',
                             section=section2,
                             phenotype='likely susceptible',
                             comment=None),
                   Namespace(mutation='T40A',
                             section=section2,
                             phenotype='resistance possible',
                             comment=None),
                   Namespace(mutation='WT',
                             section=section3,
                             phenotype='likely susceptible',
                             comment=None),
                   Namespace(mutation='T20A',
                             section=section3,
                             phenotype='resistance possible',
                             comment=None)]
        expected_rules = """\
- code: PTV
  genotypes:
  - genotype: 1A
    reference: HCV1A-H77-NS3
    region: NS3
    rules: SCORE FROM ( T40A => 4 )
  name: Paritaprevir
"""
        self.expected_errors = """\
Unknown drugs:
  NS3_GT1a: Paulrevir
  NS3_GT1a: Saynomorevir
"""

        self.assertWrites(expected_rules, entries)

    def test_wild_type_resistant(self):
        section = Namespace(drug_name='Paritaprevir', sheet_name='NS3_GT1a')
        entries = [Namespace(mutation='WT',
                             section=section,
                             phenotype='resistance possible',
                             comment=None),
                   Namespace(mutation='T40A',
                             section=section,
                             phenotype='resistance possible',
                             comment=None)]
        expected_rules = """\
- code: PTV
  genotypes:
  - genotype: 1A
    reference: HCV1A-H77-NS3
    region: NS3
    rules: SCORE FROM ( TRUE => 4, T40A => 4 )
  name: Paritaprevir
"""

        self.assertWrites(expected_rules, entries)

    def test_wild_type_mismatch(self):
        section = Namespace(drug_name='Paritaprevir', sheet_name='NS3_GT3')
        entries = [Namespace(mutation='WT',
                             section=section,
                             phenotype='likely susceptible'),
                   Namespace(mutation='A156G',
                             section=section,
                             phenotype='resistance likely'),
                   Namespace(mutation='Q186L',
                             section=section,
                             phenotype='resistance possible')]
        expected_rules = """\
- code: PTV
  genotypes:
  - genotype: '3'
    reference: HCV3-S52-NS3
    region: NS3
    rules: SCORE FROM ( A156G => 8 )
  name: Paritaprevir
"""
        self.expected_errors = """\
Mismatched wild type: NS3_GT3: Q186L in Paritaprevir expected D.
"""

        self.assertWrites(expected_rules, entries)

    def test_invalid_mutation(self):
        section = Namespace(drug_name='Paritaprevir', sheet_name='NS3_GT1a')
        entries = [Namespace(mutation='WT',
                             section=section,
                             phenotype='likely susceptible',
                             comment=None),
                   Namespace(mutation='T30A?',
                             section=section,
                             phenotype='resistance possible',
                             comment=None),
                   Namespace(mutation='T40A',
                             section=section,
                             phenotype='resistance possible',
                             comment=None)]
        expected_rules = """\
- code: PTV
  genotypes:
  - genotype: 1A
    reference: HCV1A-H77-NS3
    region: NS3
    rules: SCORE FROM ( T40A => 4 )
  name: Paritaprevir
"""
        self.expected_errors = """\
Invalid mutation: NS3_GT1a: T30A? (MutationSet text expects wild type \
(optional), position, and one or more variants.).
"""

        self.assertWrites(expected_rules, entries)

    def test_invalid_mutations(self):
        section = Namespace(drug_name='Paritaprevir', sheet_name='NS3_GT1a')
        entries = [Namespace(mutation='WT',
                             section=section,
                             phenotype='likely susceptible',
                             comment=None),
                   Namespace(mutation='T30A?',
                             section=section,
                             phenotype='resistance possible',
                             comment=None),
                   Namespace(mutation='T30W?',
                             section=section,
                             phenotype='resistance possible',
                             comment=None),
                   Namespace(mutation='T40A',
                             section=section,
                             phenotype='resistance possible',
                             comment=None)]
        expected_rules = """\
- code: PTV
  genotypes:
  - genotype: 1A
    reference: HCV1A-H77-NS3
    region: NS3
    rules: SCORE FROM ( T40A => 4 )
  name: Paritaprevir
"""
        self.expected_errors = """\
Invalid mutations:
  NS3_GT1a: T30A? (MutationSet text expects wild type \
(optional), position, and one or more variants.)
  NS3_GT1a: T30W? (MutationSet text expects wild type \
(optional), position, and one or more variants.)
"""

        self.assertWrites(expected_rules, entries)

    def test_combination_matches(self):
        section = Namespace(drug_name='Paritaprevir', sheet_name='NS3_GT1a')
        entries = [Namespace(mutation='WT',
                             section=section,
                             phenotype='likely susceptible',
                             comment=None),
                   Namespace(mutation='E30A',
                             section=section,
                             phenotype='resistance possible',
                             comment=None),
                   Namespace(mutation='T40W',
                             section=section,
                             phenotype='resistance possible',
                             comment=None),
                   Namespace(mutation='E30A+T40W',
                             section=section,
                             phenotype='resistance likely',
                             comment=None)]
        expected_rules = """\
- code: PTV
  genotypes:
  - genotype: 1A
    reference: HCV1A-H77-NS3
    region: NS3
    rules: SCORE FROM ( E30A => 4, T40W => 4 )
  name: Paritaprevir
"""

        self.assertWrites(expected_rules, entries)

    def test_combination_change(self):
        """ Combination score differs from parts, just report it for now. """
        section = Namespace(drug_name='Paritaprevir', sheet_name='NS3_GT1a')
        entries = [Namespace(mutation='WT',
                             section=section,
                             phenotype='likely susceptible'),
                   Namespace(mutation='E30A',
                             section=section,
                             phenotype='resistance possible'),
                   Namespace(mutation='E30A+R40W',
                             section=section,
                             phenotype='resistance likely')]
        expected_rules = """\
- code: PTV
  genotypes:
  - genotype: 1A
    reference: HCV1A-H77-NS3
    region: NS3
    rules: SCORE FROM ( E30A => 4 )
  name: Paritaprevir
"""
        self.expected_errors = """\
Combination change: NS3_GT1a: E30A+R40W: 4 => 8.
"""

        self.assertWrites(expected_rules, entries)

    def test_invalid_combination(self):
        section = Namespace(drug_name='Paritaprevir', sheet_name='NS3_GT1a')
        entries = [Namespace(mutation='WT',
                             section=section,
                             phenotype='likely susceptible',
                             comment=None),
                   Namespace(mutation='E30A',
                             section=section,
                             phenotype='resistance possible',
                             comment=None),
                   Namespace(mutation='E30A+R40W*',
                             section=section,
                             phenotype='resistance likely',
                             comment=None)]
        expected_rules = """\
- code: PTV
  genotypes:
  - genotype: 1A
    reference: HCV1A-H77-NS3
    region: NS3
    rules: SCORE FROM ( E30A => 4 )
  name: Paritaprevir
"""
        self.expected_errors = """\
Invalid mutation: NS3_GT1a: E30A+R40W* (MutationSet text expects wild type \
(optional), position, and one or more variants.).
"""

        self.assertWrites(expected_rules, entries)

    def test_monitored_positions(self):
        section = Namespace(drug_name='Paritaprevir',
                            sheet_name='NS3_GT1a',
                            monitored_positions=[80])
        entries = [Namespace(mutation='WT',
                             section=section,
                             phenotype='likely susceptible'),
                   Namespace(mutation='Q80K',
                             section=section,
                             phenotype='resistance possible'),
                   Namespace(mutation='Q80L',
                             section=section,
                             phenotype='resistance likely')]
        expected_rules = """\
- code: PTV
  genotypes:
  - genotype: 1A
    reference: HCV1A-H77-NS3
    region: NS3
    rules: SCORE FROM ( Q80K => 4, Q80L => 8, Q80!KLQ => "Effect unknown" )
  name: Paritaprevir
"""

        self.assertWrites(expected_rules, entries)

    def test_monitored_position_invalid(self):
        section = Namespace(drug_name='Paritaprevir',
                            sheet_name='NS3_GT1a',
                            monitored_positions=[800])
        entries = [Namespace(mutation='WT',
                             section=section,
                             phenotype='likely susceptible'),
                   Namespace(mutation='Q80K',
                             section=section,
                             phenotype='resistance possible'),
                   Namespace(mutation='Q80L',
                             section=section,
                             phenotype='resistance likely')]
        expected_rules = """\
- code: PTV
  genotypes:
  - genotype: 1A
    reference: HCV1A-H77-NS3
    region: NS3
    rules: SCORE FROM ( Q80K => 4, Q80L => 8 )
  name: Paritaprevir
"""
        self.expected_errors = """\
Invalid monitored position: NS3_GT1a: Paritaprevir 800 (max 631).
"""

        self.assertWrites(expected_rules, entries)

    def test_mutation_invalid_position(self):
        section = Namespace(drug_name='Paritaprevir',
                            sheet_name='NS3_GT1a')
        entries = [Namespace(mutation='WT',
                             section=section,
                             phenotype='likely susceptible'),
                   Namespace(mutation='Q80K',
                             section=section,
                             phenotype='resistance possible'),
                   Namespace(mutation='Q800L',
                             section=section,
                             phenotype='resistance likely')]
        expected_rules = """\
- code: PTV
  genotypes:
  - genotype: 1A
    reference: HCV1A-H77-NS3
    region: NS3
    rules: SCORE FROM ( Q80K => 4 )
  name: Paritaprevir
"""
        self.expected_errors = """\
Invalid mutation position: NS3_GT1a: Paritaprevir Q800L (max 631).
"""

        self.assertWrites(expected_rules, entries)


class WriteRulesTest(TestCase):
    def test_single(self):
        rule_sets = [RuleSet('NS3', '1a', 'Paritaprevir', 1, 2, 3, {'R10V': 4})]
        expected_rules_text = """\
- code: PTV
  genotypes:
  - genotype: 1A
    reference: HCV1A-H77-NS3
    region: NS3
    rules: SCORE FROM ( R10V => 4 )
  name: Paritaprevir
"""
        rules_file = StringIO()

        write_rules(rule_sets, REFERENCES, rules_file)

        self.assertEqual(expected_rules_text, rules_file.getvalue())

    def test_resistant_wild_type(self):
        rule_sets = [RuleSet('NS3', '1a', 'Paritaprevir', 1, 2, 3, {None: 4})]
        expected_rules_text = """\
- code: PTV
  genotypes:
  - genotype: 1A
    reference: HCV1A-H77-NS3
    region: NS3
    rules: SCORE FROM ( TRUE => 4 )
  name: Paritaprevir
"""
        rules_file = StringIO()

        write_rules(rule_sets, REFERENCES, rules_file)

        self.assertEqual(expected_rules_text, rules_file.getvalue())

    def test_resistant_wild_type_plus_other(self):
        rule_sets = [RuleSet('NS3', '1a', 'Paritaprevir', 1, 2, 3, {None: 4,
                                                                    'R10V': 4})]
        expected_rules_text = """\
- code: PTV
  genotypes:
  - genotype: 1A
    reference: HCV1A-H77-NS3
    region: NS3
    rules: SCORE FROM ( TRUE => 4, R10V => 4 )
  name: Paritaprevir
"""
        rules_file = StringIO()

        write_rules(rule_sets, REFERENCES, rules_file)

        self.assertEqual(expected_rules_text, rules_file.getvalue())

    def test_empty_rule_set(self):
        rule_sets = [RuleSet('NS3', '1a', 'Paritaprevir', 1, 2, 3, {})]
        expected_rules_text = """\
- code: PTV
  genotypes:
  - genotype: 1A
    reference: HCV1A-H77-NS3
    region: NS3
    rules: SCORE FROM ( TRUE => 0 )
  name: Paritaprevir
"""
        rules_file = StringIO()

        write_rules(rule_sets, REFERENCES, rules_file)

        self.assertEqual(expected_rules_text, rules_file.getvalue())

    def test_two_genotypes(self):
        rule_sets = [RuleSet('NS3', '1a', 'Paritaprevir', 1, 2, 3, {'R10V': 4}),
                     RuleSet('NS3', '1b', 'Paritaprevir', 1, 2, 3, {'A20L': 8})]
        expected_rules_text = """\
- code: PTV
  genotypes:
  - genotype: 1A
    reference: HCV1A-H77-NS3
    region: NS3
    rules: SCORE FROM ( R10V => 4 )
  - genotype: 1B
    reference: HCV1B-Con1-NS3
    region: NS3
    rules: SCORE FROM ( A20L => 8 )
  name: Paritaprevir
"""
        rules_file = StringIO()

        write_rules(rule_sets, REFERENCES, rules_file)

        self.assertEqual(expected_rules_text, rules_file.getvalue())

    def test_two_drugs(self):
        rule_sets = [RuleSet('NS3', '1a', 'Paritaprevir', 1, 2, 3, {'R10V': 4}),
                     RuleSet('NS3', '1b', 'Boceprevir', 1, 2, 3, {'A20L': 8})]
        expected_rules_text = """\
- code: BPV
  genotypes:
  - genotype: 1B
    reference: HCV1B-Con1-NS3
    region: NS3
    rules: SCORE FROM ( A20L => 8 )
  name: Boceprevir
- code: PTV
  genotypes:
  - genotype: 1A
    reference: HCV1A-H77-NS3
    region: NS3
    rules: SCORE FROM ( R10V => 4 )
  name: Paritaprevir
"""
        rules_file = StringIO()

        write_rules(rule_sets, REFERENCES, rules_file)

        self.assertEqual(expected_rules_text, rules_file.getvalue())

    def test_two_mutations(self):
        rule_sets = [RuleSet('NS3', '1a', 'Paritaprevir', 1, 2, 3, {'R10V': 4,
                                                                    'A20L': 8})]
        expected_rules_text = """\
- code: PTV
  genotypes:
  - genotype: 1A
    reference: HCV1A-H77-NS3
    region: NS3
    rules: SCORE FROM ( R10V => 4, A20L => 8 )
  name: Paritaprevir
"""
        rules_file = StringIO()

        write_rules(rule_sets, REFERENCES, rules_file)

        self.assertEqual(expected_rules_text, rules_file.getvalue())

    def test_two_variants(self):
        rule_sets = [RuleSet('NS3', '1a', 'Paritaprevir', 1, 2, 3, {'R10A': 4,
                                                                    'R10V': 4})]
        expected_rules_text = """\
- code: PTV
  genotypes:
  - genotype: 1A
    reference: HCV1A-H77-NS3
    region: NS3
    rules: SCORE FROM ( R10AV => 4 )
  name: Paritaprevir
"""
        rules_file = StringIO()

        write_rules(rule_sets, REFERENCES, rules_file)

        self.assertEqual(expected_rules_text, rules_file.getvalue())

    def test_combination_unchanged(self):
        rule_sets = [RuleSet('NS3',
                             '1a',
                             'Paritaprevir',
                             1,
                             2,
                             3,
                             {'R10A': 4, 'A20L': 4, 'R10A+A20L': 8})]
        expected_rules_text = """\
- code: PTV
  genotypes:
  - genotype: 1A
    reference: HCV1A-H77-NS3
    region: NS3
    rules: SCORE FROM ( R10A => 4, A20L => 4 )
  name: Paritaprevir
"""
        rules_file = StringIO()

        write_rules(rule_sets, REFERENCES, rules_file)

        self.assertEqual(expected_rules_text, rules_file.getvalue())

    @skip('HCV mutation combinations not implemented yet.')
    def test_combination_weakened(self):
        rule_sets = [RuleSet('NS3',
                             '1a',
                             'Paritaprevir',
                             1,
                             2,
                             3,
                             {'R10A': 4, 'A20L': 4, 'R10A+A20L': 4, 'T30V': 8})]
        expected_rules_text = """\
- code: PTV
  genotypes:
  - genotype: 1A
    reference: HCV1A-H77-NS3
    region: NS3
    rules: SCORE FROM ( T30V => 8, MIN( R10A => 4, A20L => 4 ) )
  name: Paritaprevir
"""
        rules_file = StringIO()

        write_rules(rule_sets, REFERENCES, rules_file)

        self.assertEqual(expected_rules_text, rules_file.getvalue())

    def test_combination_with_susceptible_components(self):
        rule_sets = [RuleSet('NS3',
                             '1a',
                             'Paritaprevir',
                             1,
                             2,
                             3,
                             {'T30V': 8, 'R10A+A20L': 8})]
        expected_rules_text = """\
- code: PTV
  genotypes:
  - genotype: 1A
    reference: HCV1A-H77-NS3
    region: NS3
    rules: SCORE FROM ( T30V => 8 )
  name: Paritaprevir
"""
        rules_file = StringIO()

        write_rules(rule_sets, REFERENCES, rules_file)

        self.assertEqual(expected_rules_text, rules_file.getvalue())

    def test_combination_bad(self):
        rule_sets = [RuleSet('NS3',
                             '1a',
                             'Paritaprevir',
                             1,
                             2,
                             3,
                             {'R10A': 4, 'A20L': 4, 'V10A+A20L': 4})]
        rules_file = StringIO()

        with self.assertRaisesRegex(
                ValueError,
                r"Components score could not be calculated for NS3 1a V10A\+A20L\."):
            write_rules(rule_sets, REFERENCES, rules_file)

    def test_unknown_drug(self):
        rule_sets = [RuleSet('NS3',
                             '1a',
                             'Paulrevir',
                             1,
                             2,
                             3,
                             {'R10A': 4, 'A20L': 4, 'V10A+A20L': 4}),
                     RuleSet('NS3',
                             '1b',
                             'Paulrevir',
                             1,
                             2,
                             3,
                             {'R10A': 4, 'A20L': 4, 'V10A+A20L': 4})]
        rules_file = StringIO()

        with self.assertRaisesRegex(
                ValueError,
                r"Unknown drug: Paulrevir found in NS3 1a\."):
            write_rules(rule_sets, REFERENCES, rules_file)

    def test_line_wrap(self):
        rule_sets = [RuleSet('NS3', '1a', 'Paritaprevir', 1, 2, 3, {'R1A': 4,
                                                                    'R2A': 4,
                                                                    'R3A': 4,
                                                                    'R4A': 4,
                                                                    'R5A': 4,
                                                                    'R6ACLNV': 4,
                                                                    'R7A': 4})]
        expected_rules_text = """\
- code: PTV
  genotypes:
  - genotype: 1A
    reference: HCV1A-H77-NS3
    region: NS3
    rules: SCORE FROM ( R1A => 4, R2A => 4, R3A => 4, R4A => 4, R5A => 4,
       R6ACLNV => 4, R7A => 4 )
  name: Paritaprevir
"""
        rules_file = StringIO()

        write_rules(rule_sets, REFERENCES, rules_file)

        self.assertEqual(expected_rules_text, rules_file.getvalue())

    def test_text_score(self):
        rule_sets = [RuleSet('NS3', '1a', 'Paritaprevir', 1, 2, 3, {'R1A': 4,
                                                                    'R2A': 'Effect unknown'})]
        expected_rules_text = """\
- code: PTV
  genotypes:
  - genotype: 1A
    reference: HCV1A-H77-NS3
    region: NS3
    rules: SCORE FROM ( R1A => 4, R2A => "Effect unknown" )
  name: Paritaprevir
"""
        rules_file = StringIO()

        write_rules(rule_sets, REFERENCES, rules_file)

        self.assertEqual(expected_rules_text, rules_file.getvalue())

    def test_text_score_with_quotes(self):
        rule_sets = [RuleSet('NS3',
                             '1a',
                             'Paritaprevir',
                             1,
                             2,
                             3,
                             {'R1A': 4, 'R2A': 'with "quotes"'})]
        rules_file = StringIO()

        with self.assertRaisesRegex(
                ValueError,
                r'Text score contained quotes: with "quotes".'):
            write_rules(rule_sets, REFERENCES, rules_file)

    def test_bad_format(self):
        rule_sets = [RuleSet('R1',
                             '1a',
                             'Paritaprevir',
                             1,
                             2,
                             3,
                             {'bad-mutation-format': 4})]
        rules_file = StringIO()

        with self.assertRaisesRegex(
                ValueError,
                'Unable to parse mutation: bad-mutation-format'):
            write_rules(rule_sets, REFERENCES, rules_file)
