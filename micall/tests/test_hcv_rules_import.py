from unittest import TestCase

from io import StringIO

from micall.utils.hcv_rules_import import create_rule_set, RuleSet, write_rules


class CreateRuleSetTest(TestCase):
    def test(self):
        sheet_name = 'R1_GT2x'
        drug_name = 'Paritaprevir in GT1a'
        phenotype_column = 23
        expected_region = 'R1'
        expected_genotype = '2x'
        expected_drug_name = 'Paritaprevir'

        rule_set = create_rule_set(sheet_name, drug_name, phenotype_column)

        self.assertEqual(phenotype_column, rule_set.phenotype_column)
        self.assertEqual(expected_region, rule_set.region)
        self.assertEqual(expected_genotype, rule_set.genotype)
        self.assertEqual(expected_drug_name, rule_set.drug_name)


class WriteRulesTest(TestCase):
    def test_single(self):
        rule_sets = [RuleSet('NS3', '1a', 'Paritaprevir', 99, {'R10V': 4})]
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

        write_rules(rule_sets, rules_file)

        self.assertEqual(expected_rules_text, rules_file.getvalue())

    def test_two_genotypes(self):
        rule_sets = [RuleSet('NS3', '1a', 'Paritaprevir', 99, {'R10V': 4}),
                     RuleSet('NS3', '1b', 'Paritaprevir', 99, {'A20L': 8})]
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

        write_rules(rule_sets, rules_file)

        self.assertEqual(expected_rules_text, rules_file.getvalue())

    def test_two_drugs(self):
        rule_sets = [RuleSet('NS3', '1a', 'Paritaprevir', 99, {'R10V': 4}),
                     RuleSet('NS3', '1b', 'Boceprevir', 99, {'A20L': 8})]
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

        write_rules(rule_sets, rules_file)

        self.assertEqual(expected_rules_text, rules_file.getvalue())

    def test_two_mutations(self):
        rule_sets = [RuleSet('NS3', '1a', 'Paritaprevir', 99, {'R10V': 4, 'A20L': 8})]
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

        write_rules(rule_sets, rules_file)

        self.assertEqual(expected_rules_text, rules_file.getvalue())

    def test_two_variants(self):
        rule_sets = [RuleSet('NS3', '1a', 'Paritaprevir', 99, {'R10A': 4, 'R10V': 4})]
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

        write_rules(rule_sets, rules_file)

        self.assertEqual(expected_rules_text, rules_file.getvalue())

    def test_combination_unchanged(self):
        rule_sets = [RuleSet('NS3', '1a', 'Paritaprevir', 99, {'R10A': 4,
                                                              'A20L': 4,
                                                              'R10A+A20L': 8})]
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

        write_rules(rule_sets, rules_file)

        self.assertEqual(expected_rules_text, rules_file.getvalue())

    def test(self):
        rule_sets = [RuleSet('NS3', '1a', 'Paritaprevir', 99, {'R1A': 4,
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

        write_rules(rule_sets, rules_file)

        self.assertEqual(expected_rules_text, rules_file.getvalue())

    def test_bad_format(self):
        rule_sets = [RuleSet('R1',
                             '1a',
                             'Paritaprevir',
                             99,
                             {'bad-mutation-format': 4})]
        rules_file = StringIO()

        with self.assertRaisesRegex(
                ValueError,
                'Unable to parse mutation: bad-mutation-format'):
            write_rules(rule_sets, rules_file)
