from argparse import Namespace
from copy import copy
from unittest import TestCase

from io import StringIO

from openpyxl import Workbook

from micall.utils.hcv_rules_import import load_references, WorksheetReader, \
    MonitoredPositionsReader, FoldRangesReader, RulesWriter

REFERENCES = load_references()


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
        expected_section = Namespace(drug_name='Example',
                                     sheet_name='NS3_GT1b',
                                     not_indicated=True)
        expected_entries = [Namespace(section=expected_section)]

        self.assertReads(expected_entries, worksheets)

    def test_fold_shift_levels(self):
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
                                     upper_fold=100,
                                     lower_range_fold=20,
                                     upper_range_fold=100)
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
                                     upper_fold=100,
                                     lower_range_fold=20,
                                     upper_range_fold=100)
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
                                     upper_fold=100,
                                     lower_range_fold=20,
                                     upper_range_fold=100)
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

    def test_fold_shift_ranges(self):
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
            ['',     '>100x FS, resistance likely'],
            ['',     'Range FS values -'],
            ['',     '2.5-10x FS range, likely susceptible'],
            ['',     '10-100x FS range, resistance possible'],
            ['',     '100-1000x FS range, resistance likely']]
        worksheets = [create_worksheet('NS3_GT1b', row_data)]
        expected_section = Namespace(drug_name='Example',
                                     sheet_name='NS3_GT1b',
                                     lower_fold=20,
                                     upper_fold=100,
                                     lower_range_fold=10,
                                     upper_range_fold=100)
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

    def test_fold_shift_range_limits(self):
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
            ['',     '>100x FS, resistance likely'],
            ['',     'Range FS values -'],
            ['',     '<10x FS range, likely susceptible'],
            ['',     '10-100x FS range, resistance possible'],
            ['',     '>100x FS range, resistance likely']]
        worksheets = [create_worksheet('NS3_GT1b', row_data)]
        expected_section = Namespace(drug_name='Example',
                                     sheet_name='NS3_GT1b',
                                     lower_fold=20,
                                     upper_fold=100,
                                     lower_range_fold=10,
                                     upper_range_fold=100)
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
                                     upper_fold=None,
                                     lower_range_fold=None,
                                     upper_range_fold=None)
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
                             phenotype='likely susceptible'),
                   Namespace(mutation='T40A',
                             section=section,
                             phenotype='resistance possible')]
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

    def test_one_position_two_mutations(self):
        section = Namespace(drug_name='Paritaprevir', sheet_name='NS3_GT1a')
        entries = [Namespace(mutation='WT',
                             section=section,
                             phenotype='likely susceptible'),
                   Namespace(mutation='T40A',
                             section=section,
                             phenotype='resistance possible'),
                   Namespace(mutation='T40R',
                             section=section,
                             phenotype='resistance possible')]
        expected_rules = """\
- code: PTV
  genotypes:
  - genotype: 1A
    reference: HCV1A-H77-NS3
    region: NS3
    rules: SCORE FROM ( T40AR => 4 )
  name: Paritaprevir
"""

        self.assertWrites(expected_rules, entries)

    def test_no_resistance(self):
        section = Namespace(drug_name='Paritaprevir', sheet_name='NS3_GT1a')
        entries = [Namespace(mutation='WT',
                             section=section,
                             phenotype='likely susceptible')]
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

    def test_not_indicated(self):
        section = Namespace(drug_name='Paritaprevir',
                            sheet_name='NS3_GT1a',
                            not_indicated=True)
        entries = [Namespace(section=section)]
        expected_rules = """\
- code: PTV
  genotypes:
  - genotype: 1A
    reference: HCV1A-H77-NS3
    region: NS3
    rules: SCORE FROM ( TRUE => "Not indicated" )
  name: Paritaprevir
"""

        self.assertWrites(expected_rules, entries)

    def test_missing_genotypes(self):
        section1 = Namespace(drug_name='Grazoprevir', sheet_name='NS3_GT1a')
        section2 = Namespace(drug_name='Paritaprevir', sheet_name='NS3_GT2')
        entries = [Namespace(mutation='WT',
                             section=section1,
                             phenotype='likely susceptible'),
                   Namespace(mutation='WT',
                             section=section2,
                             phenotype='likely susceptible')]
        expected_rules = """\
- code: GZR
  genotypes:
  - genotype: 1A
    reference: HCV1A-H77-NS3
    region: NS3
    rules: SCORE FROM ( TRUE => 0 )
  - genotype: '2'
    reference: HCV2-JFH-1-NS3
    region: NS3
    rules: SCORE FROM ( TRUE => "Not indicated" )
  name: Grazoprevir
- code: PTV
  genotypes:
  - genotype: 1A
    reference: HCV1A-H77-NS3
    region: NS3
    rules: SCORE FROM ( TRUE => "Not indicated" )
  - genotype: '2'
    reference: HCV2-JFH-1-NS3
    region: NS3
    rules: SCORE FROM ( TRUE => 0 )
  name: Paritaprevir
"""

        self.assertWrites(expected_rules, entries)

    def test_exclude_zeroes(self):
        section = Namespace(drug_name='Paritaprevir', sheet_name='NS3_GT1a')
        entries = [Namespace(mutation='WT',
                             section=section,
                             phenotype='likely susceptible'),
                   Namespace(mutation='R20A',
                             section=section,
                             phenotype='likely susceptible'),
                   Namespace(mutation='T40A',
                             section=section,
                             phenotype='resistance likely')]
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
                             phenotype='likely susceptible'),
                   Namespace(mutation='T40A',
                             section=section,
                             phenotype='resiSTANce possible')]
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
                             phenotype='likely susceptible'),
                   Namespace(mutation='R20S',
                             section=section,
                             phenotype='bogus resistance'),
                   Namespace(mutation='T40A',
                             section=section,
                             phenotype='resistance possible')]
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
                             phenotype='likely susceptible'),
                   Namespace(mutation='T41R',
                             section=section1,
                             phenotype='resistance likely'),
                   Namespace(mutation='WT',
                             section=section2,
                             phenotype='likely susceptible'),
                   Namespace(mutation='T40A',
                             section=section2,
                             phenotype='resistance possible')]
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
                             phenotype='likely susceptible'),
                   Namespace(mutation='T41R',
                             section=section1,
                             phenotype='resistance likely'),
                   Namespace(mutation='WT',
                             section=section2,
                             phenotype='likely susceptible'),
                   Namespace(mutation='T40A',
                             section=section2,
                             phenotype='resistance possible'),
                   Namespace(mutation='WT',
                             section=section3,
                             phenotype='likely susceptible'),
                   Namespace(mutation='T20A',
                             section=section3,
                             phenotype='resistance possible')]
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
                             phenotype='resistance possible'),
                   Namespace(mutation='T40A',
                             section=section,
                             phenotype='resistance possible')]
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
                             phenotype='likely susceptible'),
                   Namespace(mutation='T30A?',
                             section=section,
                             phenotype='resistance possible'),
                   Namespace(mutation='T40A',
                             section=section,
                             phenotype='resistance possible')]
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
                             phenotype='likely susceptible'),
                   Namespace(mutation='T30A?',
                             section=section,
                             phenotype='resistance possible'),
                   Namespace(mutation='T30W?',
                             section=section,
                             phenotype='resistance possible'),
                   Namespace(mutation='T40A',
                             section=section,
                             phenotype='resistance possible')]
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
                             phenotype='likely susceptible'),
                   Namespace(mutation='E30A',
                             section=section,
                             phenotype='resistance possible'),
                   Namespace(mutation='T40W',
                             section=section,
                             phenotype='resistance possible'),
                   Namespace(mutation='E30A+T40W',
                             section=section,
                             phenotype='resistance likely')]
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
                             phenotype='likely susceptible'),
                   Namespace(mutation='E30A',
                             section=section,
                             phenotype='resistance possible'),
                   Namespace(mutation='E30A+R40W*',
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

    def test_monitored_positions_long(self):
        section = Namespace(drug_name='Paritaprevir',
                            sheet_name='NS3_GT1a',
                            monitored_positions=[80, 81])
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
    rules: SCORE FROM ( Q80K => 4, Q80L => 8, Q80!KLQ => "Effect unknown",
       D81!D => "Effect unknown" )
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

    def test_fold_shifts_match(self):
        section = Namespace(drug_name='Paritaprevir',
                            sheet_name='NS3_GT1a',
                            lower_fold=20.0,
                            upper_fold=100.0)
        entries = [Namespace(mutation='WT',
                             section=section,
                             phenotype='likely susceptible',
                             fold_shift='1x',
                             clinical_ras=None),
                   Namespace(mutation='Q80K',
                             section=section,
                             phenotype='resistance possible',
                             fold_shift='50x',
                             clinical_ras=None),
                   Namespace(mutation='Q80L',
                             section=section,
                             phenotype='resistance likely',
                             fold_shift='200x',
                             clinical_ras=None)]
        expected_rules = """\
- code: PTV
  genotypes:
  - genotype: 1A
    reference: HCV1A-H77-NS3
    region: NS3
    rules: SCORE FROM ( Q80K => 4, Q80L => 8 )
  name: Paritaprevir
"""

        self.assertWrites(expected_rules, entries)

    def test_fold_shifts_differ(self):
        section = Namespace(drug_name='Paritaprevir',
                            sheet_name='NS3_GT1a',
                            lower_fold=20.0,
                            upper_fold=100.0)
        entries = [Namespace(mutation='WT',
                             section=section,
                             phenotype='likely susceptible',
                             fold_shift='1x',
                             clinical_ras=None),
                   Namespace(mutation='Q80K',
                             section=section,
                             phenotype='resistance possible',
                             fold_shift='2x',
                             clinical_ras=None),
                   Namespace(mutation='Q80L',
                             section=section,
                             phenotype='resistance likely',
                             fold_shift='100x',  # Should be > 100.0
                             clinical_ras=None)]
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
Phenotype changes:
  NS3_GT1a: Paritaprevir Q80K likely susceptible => resistance possible
  NS3_GT1a: Paritaprevir Q80L resistance possible => resistance likely
"""

        self.assertWrites(expected_rules, entries)

    def test_fold_shift_no_levels(self):
        section = Namespace(drug_name='Paritaprevir',
                            sheet_name='NS3_GT1a',
                            lower_fold=None,
                            upper_fold=None)
        entries = [Namespace(mutation='WT',
                             section=section,
                             phenotype='likely susceptible',
                             fold_shift='1x',
                             clinical_ras=None),
                   Namespace(mutation='Q80K',
                             section=section,
                             phenotype='resistance possible',
                             fold_shift='2x',
                             clinical_ras=None),
                   Namespace(mutation='Q80L',
                             section=section,
                             phenotype='resistance likely',
                             fold_shift='100x',
                             clinical_ras=None)]
        expected_rules = """\
- code: PTV
  genotypes:
  - genotype: 1A
    reference: HCV1A-H77-NS3
    region: NS3
    rules: SCORE FROM ( Q80K => 4, Q80L => 8 )
  name: Paritaprevir
"""

        self.assertWrites(expected_rules, entries)

    def test_fold_shift_ranges(self):
        section = Namespace(drug_name='Paritaprevir',
                            sheet_name='NS3_GT1a',
                            lower_fold=10,
                            upper_fold=100,
                            lower_range_fold=50,
                            upper_range_fold=100)
        entries = [Namespace(mutation='WT',
                             section=section,
                             phenotype='likely susceptible',
                             fold_shift='1x',
                             clinical_ras=None),
                   Namespace(mutation='Q80K',
                             section=section,
                             phenotype='likely susceptible',
                             fold_shift='2-50x',
                             clinical_ras=None),
                   Namespace(mutation='Q80L',
                             section=section,
                             phenotype='resistance likely',
                             fold_shift='100-1000x',
                             clinical_ras=None)]
        expected_rules = """\
- code: PTV
  genotypes:
  - genotype: 1A
    reference: HCV1A-H77-NS3
    region: NS3
    rules: SCORE FROM ( Q80L => 8 )
  name: Paritaprevir
"""

        self.assertWrites(expected_rules, entries)

    def test_fold_shift_range_spans_levels(self):
        section = Namespace(drug_name='Paritaprevir',
                            sheet_name='NS3_GT1a',
                            lower_fold=10,
                            upper_fold=100,
                            lower_range_fold=50,
                            upper_range_fold=100)
        entries = [Namespace(mutation='WT',
                             section=section,
                             phenotype='likely susceptible',
                             fold_shift='1x',
                             clinical_ras=None),
                   Namespace(mutation='Q80K',
                             section=section,
                             phenotype='likely susceptible',
                             fold_shift='2-50x',
                             clinical_ras=None),
                   Namespace(mutation='Q80L',
                             section=section,
                             phenotype='resistance likely',
                             fold_shift='20-1000x',
                             clinical_ras=None)]
        expected_rules = """\
- code: PTV
  genotypes:
  - genotype: 1A
    reference: HCV1A-H77-NS3
    region: NS3
    rules: SCORE FROM ( Q80L => 8 )
  name: Paritaprevir
"""

        self.assertWrites(expected_rules, entries)

    def test_fold_shift_int(self):
        section = Namespace(drug_name='Paritaprevir',
                            sheet_name='NS3_GT1a',
                            lower_fold=20.0,
                            upper_fold=100.0)
        entries = [Namespace(mutation='WT',
                             section=section,
                             phenotype='likely susceptible',
                             fold_shift=1,
                             clinical_ras=None),
                   Namespace(mutation='Q80K',
                             section=section,
                             phenotype='resistance possible',
                             fold_shift=2,
                             clinical_ras=None),
                   Namespace(mutation='Q80L',
                             section=section,
                             phenotype='resistance likely',
                             fold_shift=100,  # Should be > 100.0
                             clinical_ras=None)]
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
Phenotype changes:
  NS3_GT1a: Paritaprevir Q80K likely susceptible => resistance possible
  NS3_GT1a: Paritaprevir Q80L resistance possible => resistance likely
"""

        self.assertWrites(expected_rules, entries)

    def test_fold_shift_extra_chars(self):
        section = Namespace(drug_name='Paritaprevir',
                            sheet_name='NS3_GT1a',
                            lower_fold=20.0,
                            upper_fold=100.0)
        entries = [Namespace(mutation='WT',
                             section=section,
                             phenotype='likely susceptible',
                             fold_shift='1x',
                             clinical_ras=None),
                   Namespace(mutation='Q80K',
                             section=section,
                             phenotype='resistance possible',
                             fold_shift='<2.5x',
                             clinical_ras=None),
                   Namespace(mutation='Q80L',
                             section=section,
                             phenotype='resistance likely',
                             fold_shift='>1,000x',
                             clinical_ras=None)]
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
Phenotype change: NS3_GT1a: Paritaprevir Q80K \
likely susceptible => resistance possible.
"""

        self.assertWrites(expected_rules, entries)

    def test_fold_shift_and_clinical(self):
        section = Namespace(drug_name='Paritaprevir',
                            sheet_name='NS3_GT1a',
                            lower_fold=20.0,
                            upper_fold=100.0)
        entries = [Namespace(mutation='WT',
                             section=section,
                             phenotype='likely susceptible',
                             fold_shift='1x',
                             clinical_ras=None),
                   Namespace(mutation='Q80K',
                             section=section,
                             phenotype='resistance possible',
                             fold_shift='2x',
                             clinical_ras='Yes'),
                   Namespace(mutation='Q80L',
                             section=section,
                             phenotype='resistance likely',
                             fold_shift='200x',
                             clinical_ras='yes')]
        expected_rules = """\
- code: PTV
  genotypes:
  - genotype: 1A
    reference: HCV1A-H77-NS3
    region: NS3
    rules: SCORE FROM ( Q80K => 4, Q80L => 8 )
  name: Paritaprevir
"""

        self.assertWrites(expected_rules, entries)

    def test_fold_shift_comparison(self):
        section = Namespace(drug_name='Paritaprevir',
                            sheet_name='NS3_GT1a',
                            lower_fold=20.0,
                            upper_fold=100.0)
        entries = [Namespace(mutation='WT',
                             section=section,
                             phenotype='likely susceptible',
                             fold_shift='1x',
                             clinical_ras=None),
                   Namespace(mutation='Q80K',
                             section=section,
                             phenotype='likely susceptible',
                             fold_shift='<20x',
                             clinical_ras=None),
                   Namespace(mutation='Q80L',
                             section=section,
                             phenotype='resistance likely',
                             fold_shift='>100x',
                             clinical_ras=None)]
        expected_rules = """\
- code: PTV
  genotypes:
  - genotype: 1A
    reference: HCV1A-H77-NS3
    region: NS3
    rules: SCORE FROM ( Q80L => 8 )
  name: Paritaprevir
"""

        self.assertWrites(expected_rules, entries)

    def test_fold_shift_invalid(self):
        section = Namespace(drug_name='Paritaprevir',
                            sheet_name='NS3_GT1a',
                            lower_fold=20.0,
                            upper_fold=100.0)
        entries = [Namespace(mutation='WT',
                             section=section,
                             phenotype='likely susceptible',
                             fold_shift='1x',
                             clinical_ras=None),
                   Namespace(mutation='Q80K',
                             section=section,
                             phenotype='resistance possible',
                             fold_shift='2 times',
                             clinical_ras=None),
                   Namespace(mutation='Q80L',
                             section=section,
                             phenotype='resistance likely',
                             fold_shift='100x',  # Should be > 100.0
                             clinical_ras=None)]
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
Invalid fold shift: NS3_GT1a: Paritaprevir Q80K '2 times'.
Phenotype change: NS3_GT1a: Paritaprevir Q80L \
resistance possible => resistance likely.
"""

        self.assertWrites(expected_rules, entries)

    def test_fold_shift_missing(self):
        section = Namespace(drug_name='Paritaprevir',
                            sheet_name='NS3_GT1a',
                            lower_fold=20.0,
                            upper_fold=100.0)
        entries = [Namespace(mutation='WT',
                             section=section,
                             phenotype='likely susceptible',
                             fold_shift='1x',
                             clinical_ras=None),
                   Namespace(mutation='Q80K',
                             section=section,
                             phenotype='resistance possible',
                             fold_shift=None,
                             clinical_ras=None),
                   Namespace(mutation='G90L',
                             section=section,
                             phenotype='resistance likely',
                             fold_shift=None,
                             clinical_ras=None)]
        expected_rules = """\
- code: PTV
  genotypes:
  - genotype: 1A
    reference: HCV1A-H77-NS3
    region: NS3
    rules: SCORE FROM ( Q80K => 4, G90L => 8 )
  name: Paritaprevir
"""
        self.expected_errors = """\
Missing fold shift: NS3_GT1a: Paritaprevir Q80K, G90L.
"""

        self.assertWrites(expected_rules, entries)
