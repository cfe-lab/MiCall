from unittest import TestCase, skip

from io import StringIO

from openpyxl import Workbook

from micall.utils.hcv_rules_import import create_rule_set, RuleSet, write_rules, read_rule_sets, load_references

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

    def test(self):
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
