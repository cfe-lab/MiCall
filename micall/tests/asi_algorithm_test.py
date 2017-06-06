import os
from io import StringIO
from operator import attrgetter
from unittest import TestCase

from micall.hivdb.asi_algorithm import AsiAlgorithm, translate_complete_to_array, BNFVal


class AsiAlgorithmTest(TestCase):
    def setUp(self):
        self.asi = AsiAlgorithm(os.path.join(os.path.dirname(__file__),
                                             "..",
                                             "hivdb",
                                             "HIVDB_8.3.xml"))

    def test_interpret(self):
        aa_seq = [[amino] for amino in self.asi.stds['RT']]
        aa_seq[40] = ['L']
        compared_attrs = ('code', 'score', 'level', 'level_name')
        expected_drugs = [('3TC', 0.0, 1, 'Susceptible'),
                          ('ABC', 5.0, 1, 'Susceptible'),
                          ('AZT', 15.0, 3, 'Low-Level Resistance'),
                          ('D4T', 15.0, 3, 'Low-Level Resistance'),
                          ('DDI', 10.0, 2, 'Potential Low-Level Resistance'),
                          ('FTC', 0.0, 1, 'Susceptible'),
                          ('TDF', 5.0, 1, 'Susceptible'),
                          ('EFV', 0.0, 1, 'Susceptible'),
                          ('ETR', 0.0, 1, 'Susceptible'),
                          ('NVP', 0.0, 1, 'Susceptible'),
                          ('RPV', 0.0, 1, 'Susceptible')]
        expected_mutation_comments = [
            'M41L is a TAM that usually occurs with T215Y. In combination, '
            'M41L plus T215Y confer intermediate / high-level resistance to '
            'AZT and d4T and contribute to reduced ddI, ABC and TDF '
            'susceptibility.']

        result = self.asi.interpret(aa_seq, 'RT')

        drugs = list(map(attrgetter(*compared_attrs), result.drugs))
        self.assertEqual(expected_drugs, drugs)
        self.assertEqual(expected_mutation_comments, result.mutation_comments)

    def test_protease(self):
        aa_seq = [[amino] for amino in self.asi.stds['PR']]
        aa_seq[23] = ['I']
        compared_attrs = ('code', 'score', 'level', 'level_name')
        expected_drugs = [('ATV/r', 10.0, 2, 'Potential Low-Level Resistance'),
                          ('DRV/r', 0.0, 1, 'Susceptible'),
                          ('FPV/r', 10.0, 2, 'Potential Low-Level Resistance'),
                          ('IDV/r', 15.0, 3, 'Low-Level Resistance'),
                          ('LPV/r', 10.0, 2, 'Potential Low-Level Resistance'),
                          ('NFV', 10.0, 2, 'Potential Low-Level Resistance'),
                          ('SQV/r', 10.0, 2, 'Potential Low-Level Resistance'),
                          ('TPV/r', -5.0, 1, 'Susceptible')]
        expected_mutation_comments = [
            'L24I is a non-polymorphic mutation selected by IDV and LPV. It '
            'contributes reduced susceptibility to each PI except DRV and TPV.']

        result = self.asi.interpret(aa_seq, 'PR')

        drugs = list(map(attrgetter(*compared_attrs), result.drugs))
        self.assertEqual(expected_drugs, drugs)
        self.assertEqual(expected_mutation_comments, result.mutation_comments)

    def test_comment_filter(self):
        aa_seq = [[amino] for amino in self.asi.stds['PR']]
        aa_seq[23] = ['D']
        expected_mutation_comments = [
            'L24I is a non-polymorphic mutation selected by IDV and LPV. It '
            'contributes reduced susceptibility to each PI except DRV and TPV. '
            'L24D is a highly unusual mutation at this position.']

        result = self.asi.interpret(aa_seq, 'PR')

        self.assertEqual(expected_mutation_comments, result.mutation_comments)

    def test_gaps(self):
        aa_seq = [[]] * 23 + [['I']]
        compared_attrs = ('code', 'score', 'level', 'level_name')
        expected_drugs = [('ATV/r', 10.0, 2, 'Potential Low-Level Resistance'),
                          ('DRV/r', 0.0, 1, 'Susceptible'),
                          ('FPV/r', 10.0, 2, 'Potential Low-Level Resistance'),
                          ('IDV/r', 15.0, 3, 'Low-Level Resistance'),
                          ('LPV/r', 10.0, 2, 'Potential Low-Level Resistance'),
                          ('NFV', 10.0, 2, 'Potential Low-Level Resistance'),
                          ('SQV/r', 10.0, 2, 'Potential Low-Level Resistance'),
                          ('TPV/r', -5.0, 1, 'Susceptible')]
        expected_mutation_comments = [
            'L24I is a non-polymorphic mutation selected by IDV and LPV. It '
            'contributes reduced susceptibility to each PI except DRV and TPV.']

        result = self.asi.interpret(aa_seq, 'PR')

        drugs = list(map(attrgetter(*compared_attrs), result.drugs))
        self.assertEqual(expected_drugs, drugs)
        self.assertEqual(expected_mutation_comments, result.mutation_comments)

    def test_multiple_regions(self):
        aa_seq1 = [[]] * 40 + [['L']]
        aa_seq2 = [[]] * 23 + [['I']]
        compared_attrs = ('code', 'score', 'level', 'level_name')
        expected_drugs = [('ATV/r', 10.0, 2, 'Potential Low-Level Resistance'),
                          ('DRV/r', 0.0, 1, 'Susceptible'),
                          ('FPV/r', 10.0, 2, 'Potential Low-Level Resistance'),
                          ('IDV/r', 15.0, 3, 'Low-Level Resistance'),
                          ('LPV/r', 10.0, 2, 'Potential Low-Level Resistance'),
                          ('NFV', 10.0, 2, 'Potential Low-Level Resistance'),
                          ('SQV/r', 10.0, 2, 'Potential Low-Level Resistance'),
                          ('TPV/r', -5.0, 1, 'Susceptible')]
        expected_mutation_comments = [
            'L24I is a non-polymorphic mutation selected by IDV and LPV. It '
            'contributes reduced susceptibility to each PI except DRV and TPV.']

        self.asi.interpret(aa_seq1, 'RT')
        result = self.asi.interpret(aa_seq2, 'PR')

        drugs = list(map(attrgetter(*compared_attrs), result.drugs))
        self.assertEqual(expected_drugs, drugs)
        self.assertEqual(expected_mutation_comments, result.mutation_comments)

    def test_mutations(self):
        aa_seq = [[amino] for amino in self.asi.stds['RT']]
        self.assertEqual(['E'], aa_seq[39])
        aa_seq[39] = ['L']
        self.assertEqual(['M'], aa_seq[40])
        aa_seq[40] = ['L']
        expected_mutations = {'NRTI': ['M41L'], 'NNRTI': []}

        result = self.asi.interpret(aa_seq, 'RT')

        self.assertEqual(expected_mutations, result.mutations)

    def test_drug_classes(self):
        aa_seq = [[]]
        compared_attrs = ('code', 'name', 'drug_class')
        expected_drugs = [('3TC', 'lamivudine', 'NRTI'),
                          ('ABC', 'abacavir', 'NRTI'),
                          ('AZT', 'zidovudine', 'NRTI'),
                          ('D4T', 'stavudine', 'NRTI'),
                          ('DDI', 'didanosine', 'NRTI'),
                          ('FTC', 'emtricitabine', 'NRTI'),
                          ('TDF', 'tenofovir', 'NRTI'),
                          ('EFV', 'efavirenz', 'NNRTI'),
                          ('ETR', 'etravirine', 'NNRTI'),
                          ('NVP', 'nevirapine', 'NNRTI'),
                          ('RPV', 'rilpivirine', 'NNRTI')]

        result = self.asi.interpret(aa_seq, 'RT')

        drugs = list(map(attrgetter(*compared_attrs), result.drugs))
        self.assertEqual(expected_drugs, drugs)


class AsiAlgorithmTinyRulesTest(TestCase):
    def setUp(self):
        xml = "<ALGORITHM><DEFINITIONS/></ALGORITHM>"
        self.asi = AsiAlgorithm(StringIO(xml))

    def test_bnf_default_repr(self):
        expected_repr = "BNFVal(False, False, 0, set())"

        bnf = BNFVal(False)

        self.assertEqual(expected_repr, repr(bnf))

    def test_bnf_residue(self):
        aminos = [['L'], ['R']]
        cond = '2R|'
        expected_cond = '|'
        expected_score = 0
        expected_mutations = {'2R'}
        expected_repr = "BNFVal('|', True, 0, {'2R'})"

        bnr = self.asi.bnf_residue(cond, aminos)

        self.assertEqual(expected_cond, bnr.cond)
        self.assertTrue(bnr.truth)
        self.assertEqual(expected_score, bnr.score)
        self.assertEqual(expected_mutations, bnr.mutations)
        self.assertEqual(expected_repr, repr(bnr))

    def test_bnf_residue_multiple_choices(self):
        aminos = [['L'], ['R']]
        cond = '2RFL|'
        expected_cond = '|'
        expected_score = 0
        expected_mutations = {'2R'}

        bnr = self.asi.bnf_residue(cond, aminos)

        self.assertEqual(expected_cond, bnr.cond)
        self.assertTrue(bnr.truth)
        self.assertEqual(expected_score, bnr.score)
        self.assertEqual(expected_mutations, bnr.mutations)

    def test_bnf_residue_multiple_choices_found(self):
        aminos = [['L'], ['R', 'F', 'K']]
        cond = '2RFL|'
        expected_cond = '|'
        expected_score = 0
        expected_mutations = {'2R', '2F'}

        bnr = self.asi.bnf_residue(cond, aminos)

        self.assertEqual(expected_cond, bnr.cond)
        self.assertTrue(bnr.truth)
        self.assertEqual(expected_score, bnr.score)
        self.assertEqual(expected_mutations, bnr.mutations)

    def test_bnf_residue_negative_fails(self):
        aminos = [['L'], ['R']]
        cond = 'NOT 2R|'
        expected_cond = '|'
        expected_score = 0
        expected_mutations = set()

        bnr = self.asi.bnf_residue(cond, aminos)

        self.assertEqual(expected_cond, bnr.cond)
        self.assertFalse(bnr.truth)
        self.assertEqual(expected_score, bnr.score)
        self.assertEqual(expected_mutations, bnr.mutations)

    def test_bnf_residue_negative_passes(self):
        aminos = [['L'], ['F']]
        cond = 'NOT 2R|'
        expected_cond = '|'
        expected_score = 0
        expected_mutations = set()

        bnr = self.asi.bnf_residue(cond, aminos)

        self.assertEqual(expected_cond, bnr.cond)
        self.assertTrue(bnr.truth)
        self.assertEqual(expected_score, bnr.score)
        self.assertEqual(expected_mutations, bnr.mutations)

    def test_bnf_residue_negative_short(self):
        aminos = [['L']]
        cond = 'NOT 2R|'
        expected_cond = '|'
        expected_score = 0
        expected_mutations = set()

        bnr = self.asi.bnf_residue(cond, aminos)

        self.assertEqual(expected_cond, bnr.cond)
        self.assertFalse(bnr.truth)
        self.assertEqual(expected_score, bnr.score)
        self.assertEqual(expected_mutations, bnr.mutations)

    def test_bnf_residue_negative2_fails(self):
        aminos = [['L'], ['R']]
        cond = '2(NOT R)|'
        expected_cond = '|'
        expected_score = 0
        expected_mutations = set()

        bnr = self.asi.bnf_residue(cond, aminos)

        self.assertEqual(expected_cond, bnr.cond)
        self.assertFalse(bnr.truth)
        self.assertEqual(expected_score, bnr.score)
        self.assertEqual(expected_mutations, bnr.mutations)

    def test_bnf_residue_negative2_passes(self):
        aminos = [['L'], ['F']]
        cond = '2(NOT R)|'
        expected_cond = '|'
        expected_score = 0
        expected_mutations = set()

        bnr = self.asi.bnf_residue(cond, aminos)

        self.assertEqual(expected_cond, bnr.cond)
        self.assertTrue(bnr.truth)
        self.assertEqual(expected_score, bnr.score)
        self.assertEqual(expected_mutations, bnr.mutations)

    def test_bnf_residue_negative2_short(self):
        aminos = [['L']]
        cond = '2(NOT R)|'
        expected_cond = '|'
        expected_score = 0
        expected_mutations = set()

        bnr = self.asi.bnf_residue(cond, aminos)

        self.assertEqual(expected_cond, bnr.cond)
        self.assertFalse(bnr.truth)
        self.assertEqual(expected_score, bnr.score)
        self.assertEqual(expected_mutations, bnr.mutations)

    def test_bnf_exclude_fails(self):
        aminos = [['L'], ['R']]
        cond = 'EXCLUDE 2R|'
        expected_cond = '|'
        expected_score = 0
        expected_mutations = set()

        bnr = self.asi.bnf_excludestatement(cond, aminos)

        self.assertEqual(expected_cond, bnr.cond)
        self.assertFalse(bnr.truth)
        self.assertEqual(expected_score, bnr.score)
        self.assertEqual(expected_mutations, bnr.mutations)

    def test_bnf_exclude_passes(self):
        aminos = [['L'], ['F']]
        cond = 'EXCLUDE 2R|'
        expected_cond = '|'
        expected_score = 0
        expected_mutations = set()

        bnr = self.asi.bnf_excludestatement(cond, aminos)

        self.assertEqual(expected_cond, bnr.cond)
        self.assertTrue(bnr.truth)
        self.assertEqual(expected_score, bnr.score)
        self.assertEqual(expected_mutations, bnr.mutations)

    def test_bnf_exclude_short(self):
        aminos = [['L']]
        cond = 'EXCLUDE 2R|'
        expected_cond = '|'
        expected_score = 0
        expected_mutations = set()

        bnr = self.asi.bnf_excludestatement(cond, aminos)

        self.assertEqual(expected_cond, bnr.cond)
        self.assertFalse(bnr.truth)
        self.assertEqual(expected_score, bnr.score)
        self.assertEqual(expected_mutations, bnr.mutations)

    def test_bnf_select_at_least_less(self):
        aminos = [['L'], ['F'], ['L']]
        cond = 'SELECT ATLEAST 2 FROM (1P, 2F, 3R)|'
        expected_cond = '|'
        expected_score = 0
        expected_mutations = {'2F'}

        bnr = self.asi.bnf_selectstatement(cond, aminos)

        self.assertEqual(expected_cond, bnr.cond)
        self.assertFalse(bnr.truth)
        self.assertEqual(expected_score, bnr.score)
        self.assertEqual(expected_mutations, bnr.mutations)

    def test_bnf_select_at_least_equal(self):
        aminos = [['L'], ['F'], ['R']]
        cond = 'SELECT ATLEAST 2 FROM (1P, 2F, 3R)|'
        expected_cond = '|'
        expected_score = 0
        expected_mutations = {'2F', '3R'}

        bnr = self.asi.bnf_selectstatement(cond, aminos)

        self.assertEqual(expected_cond, bnr.cond)
        self.assertTrue(bnr.truth)
        self.assertEqual(expected_score, bnr.score)
        self.assertEqual(expected_mutations, bnr.mutations)

    def test_bnf_select_at_least_more(self):
        aminos = [['P'], ['F'], ['R']]
        cond = 'SELECT ATLEAST 2 FROM (1P, 2F, 3R)|'
        expected_cond = '|'
        expected_score = 0
        expected_mutations = {'1P', '2F', '3R'}

        bnr = self.asi.bnf_selectstatement(cond, aminos)

        self.assertEqual(expected_cond, bnr.cond)
        self.assertTrue(bnr.truth)
        self.assertEqual(expected_score, bnr.score)
        self.assertEqual(expected_mutations, bnr.mutations)

    def test_bnf_select_at_most_less(self):
        aminos = [['L'], ['F'], ['L']]
        cond = 'SELECT NOTMORETHAN 2 FROM (1P, 2F, 3R)|'
        expected_cond = '|'
        expected_score = 0
        expected_mutations = {'2F'}

        bnr = self.asi.bnf_selectstatement(cond, aminos)

        self.assertEqual(expected_cond, bnr.cond)
        self.assertTrue(bnr.truth)
        self.assertEqual(expected_score, bnr.score)
        self.assertEqual(expected_mutations, bnr.mutations)

    def test_bnf_select_at_most_equal(self):
        aminos = [['L'], ['F'], ['R']]
        cond = 'SELECT NOTMORETHAN 2 FROM (1P, 2F, 3R)|'
        expected_cond = '|'
        expected_score = 0
        expected_mutations = {'2F', '3R'}

        bnr = self.asi.bnf_selectstatement(cond, aminos)

        self.assertEqual(expected_cond, bnr.cond)
        self.assertTrue(bnr.truth)
        self.assertEqual(expected_score, bnr.score)
        self.assertEqual(expected_mutations, bnr.mutations)

    def test_bnf_select_at_most_more(self):
        aminos = [['P'], ['F'], ['R']]
        cond = 'SELECT NOTMORETHAN 2 FROM (1P, 2F, 3R)|'
        expected_cond = '|'
        expected_score = 0
        expected_mutations = {'1P', '2F', '3R'}

        bnr = self.asi.bnf_selectstatement(cond, aminos)

        self.assertEqual(expected_cond, bnr.cond)
        self.assertFalse(bnr.truth)
        self.assertEqual(expected_score, bnr.score)
        self.assertEqual(expected_mutations, bnr.mutations)

    def test_bnf_select_exactly(self):
        aminos = [['P'], ['L'], ['R']]
        cond = 'SELECT EXACTLY 2 FROM (1P, 2F, 3R)|'
        expected_cond = '|'
        expected_score = 0
        expected_mutations = {'1P', '3R'}

        bnr = self.asi.bnf_selectstatement(cond, aminos)

        self.assertEqual(expected_cond, bnr.cond)
        self.assertTrue(bnr.truth)
        self.assertEqual(expected_score, bnr.score)
        self.assertEqual(expected_mutations, bnr.mutations)

    def test_bnf_select_combination_and(self):
        aminos = [['P'], ['L'], ['R']]
        cond = 'SELECT ATLEAST 2 AND NOTMORETHAN 2 FROM (1P, 2F, 3R)|'
        expected_cond = '|'
        expected_score = 0
        expected_mutations = {'1P', '3R'}

        bnr = self.asi.bnf_selectstatement(cond, aminos)

        self.assertEqual(expected_cond, bnr.cond)
        self.assertTrue(bnr.truth)
        self.assertEqual(expected_score, bnr.score)
        self.assertEqual(expected_mutations, bnr.mutations)

    def test_bnf_select_combination_or(self):
        aminos = [['I'], ['L'], ['R']]
        cond = 'SELECT ATLEAST 3 OR NOTMORETHAN 1 FROM (1P, 2F, 3R)|'
        expected_cond = '|'
        expected_score = 0
        expected_mutations = {'3R'}

        bnr = self.asi.bnf_selectstatement(cond, aminos)

        self.assertEqual(expected_cond, bnr.cond)
        self.assertTrue(bnr.truth)
        self.assertEqual(expected_score, bnr.score)
        self.assertEqual(expected_mutations, bnr.mutations)

    def test_bnf_select_nested_parentheses(self):
        aminos = [['P'], ['F'], ['R']]
        cond = 'SELECT EXACTLY 2 FROM (1P, 2(NOT F), 3R)|'
        expected_cond = '|'
        expected_score = 0
        expected_mutations = {'1P', '3R'}

        bnr = self.asi.bnf_selectstatement(cond, aminos)

        self.assertEqual(expected_cond, bnr.cond)
        self.assertTrue(bnr.truth)
        self.assertEqual(expected_score, bnr.score)
        self.assertEqual(expected_mutations, bnr.mutations)

    def test_bnf_scoreitem_false(self):
        aminos = [['L'], ['M']]
        cond = '2R => 30|'
        expected_cond = '|'
        expected_score = 0.0
        expected_mutations = set()

        bnr = self.asi.bnf_scoreitem(cond, aminos)

        self.assertEqual(expected_cond, bnr.cond)
        self.assertFalse(bnr.truth)
        self.assertEqual(expected_score, bnr.score)
        self.assertEqual(expected_mutations, bnr.mutations)

    def test_bnf_scoreitem_true(self):
        aminos = [['L'], ['R']]
        cond = '2R => 30|'
        expected_cond = '|'
        expected_score = 30.0
        expected_mutations = {'2R'}

        bnr = self.asi.bnf_scoreitem(cond, aminos)

        self.assertEqual(expected_cond, bnr.cond)
        self.assertTrue(bnr.truth)
        self.assertEqual(expected_score, bnr.score)
        self.assertEqual(expected_mutations, bnr.mutations)

    def test_bnf_scoreitem_with_leftover(self):
        aminos = [['L'], ['R']]
        cond = '1R => 10, 2R => 30, 10R => 20, 20R => 19|'
        expected_cond = ', 2R => 30, 10R => 20, 20R => 19|'
        expected_score = 0.0
        expected_mutations = set()
        expected_repr = "BNFVal(', 2R => 3...0R => 19|', False, 0, set())"

        bnr = self.asi.bnf_scoreitem(cond, aminos)

        self.assertEqual(expected_cond, bnr.cond)
        self.assertFalse(bnr.truth)
        self.assertEqual(expected_score, bnr.score)
        self.assertEqual(expected_mutations, bnr.mutations)
        self.assertEqual(expected_repr, repr(bnr))

    def test_bnf_scorecondition(self):
        aminos = [['L'], ['R']]
        cond = 'SCORE FROM(1F => 20, 2R => 30)|'
        expected_cond = '|'
        expected_score = 30.0
        expected_mutations = {'2R'}

        bnr = self.asi.bnf_scorecondition(cond, aminos)

        self.assertEqual(expected_cond, bnr.cond)
        self.assertEqual(expected_score, bnr.score)
        self.assertEqual(expected_mutations, bnr.mutations)

    def test_bnf_scorecondition_multiple(self):
        aminos = [['F'], ['R']]
        cond = 'SCORE FROM(1F => 20, 2R => 30)|'
        expected_cond = '|'
        expected_score = 50.0
        expected_mutations = {'1F', '2R'}

        bnr = self.asi.bnf_scorecondition(cond, aminos)

        self.assertEqual(expected_cond, bnr.cond)
        self.assertEqual(expected_score, bnr.score)
        self.assertEqual(expected_mutations, bnr.mutations)

    def test_bnf_scorecondition_and(self):
        aminos = [['F'], ['R']]
        cond = 'SCORE FROM((1F AND 2R) => 10)|'
        expected_cond = '|'
        expected_score = 10.0
        expected_mutations = {'1F', '2R'}

        bnr = self.asi.bnf_scorecondition(cond, aminos)

        self.assertEqual(expected_cond, bnr.cond)
        self.assertEqual(expected_score, bnr.score)
        self.assertEqual(expected_mutations, bnr.mutations)

    def test_bnf_scorecondition_and_false(self):
        aminos = [['F'], ['R']]
        cond = 'SCORE FROM((1F AND 2F) => 10)|'
        expected_cond = '|'
        expected_score = 0
        expected_mutations = set()

        bnr = self.asi.bnf_scorecondition(cond, aminos)

        self.assertEqual(expected_cond, bnr.cond)
        self.assertEqual(expected_score, bnr.score)
        self.assertEqual(expected_mutations, bnr.mutations)

    def test_bnf_scorecondition_or_both(self):
        aminos = [['F'], ['R']]
        cond = 'SCORE FROM((1F OR 2R) => 10)|'
        expected_cond = '|'
        expected_score = 10
        expected_mutations = {'1F', '2R'}

        bnr = self.asi.bnf_scorecondition(cond, aminos)

        self.assertEqual(expected_cond, bnr.cond)
        self.assertEqual(expected_score, bnr.score)
        self.assertEqual(expected_mutations, bnr.mutations)

    def test_bnf_scorecondition_or_one(self):
        aminos = [['F'], ['L']]
        cond = 'SCORE FROM((1F OR 2R) => 10)|'
        expected_cond = '|'
        expected_score = 10
        expected_mutations = {'1F'}

        bnr = self.asi.bnf_scorecondition(cond, aminos)

        self.assertEqual(expected_cond, bnr.cond)
        self.assertEqual(expected_score, bnr.score)
        self.assertEqual(expected_mutations, bnr.mutations)

    def test_bnf_scorecondition_or_neither(self):
        aminos = [['L'], ['L']]
        cond = 'SCORE FROM((1F OR 2R) => 10)|'
        expected_cond = '|'
        expected_score = 0
        expected_mutations = set()

        bnr = self.asi.bnf_scorecondition(cond, aminos)

        self.assertEqual(expected_cond, bnr.cond)
        self.assertEqual(expected_score, bnr.score)
        self.assertEqual(expected_mutations, bnr.mutations)

    def test_bnf_scorecondition_max_both(self):
        aminos = [['F'], ['R']]
        cond = 'SCORE FROM(MAX ( 1F => 10, 2R => 20 ))|'
        expected_cond = '|'
        expected_score = 20
        expected_mutations = {'1F', '2R'}

        bnr = self.asi.bnf_scorecondition(cond, aminos)

        self.assertEqual(expected_cond, bnr.cond)
        self.assertEqual(expected_score, bnr.score)
        self.assertEqual(expected_mutations, bnr.mutations)

    def test_bnf_scorecondition_max_one(self):
        aminos = [['F'], ['L']]
        cond = 'SCORE FROM(MAX ( 1F => 10, 2R => 20 ))|'
        expected_cond = '|'
        expected_score = 10
        expected_mutations = {'1F'}

        bnr = self.asi.bnf_scorecondition(cond, aminos)

        self.assertEqual(expected_cond, bnr.cond)
        self.assertEqual(expected_score, bnr.score)
        self.assertEqual(expected_mutations, bnr.mutations)

    def test_bnf_scorecondition_max_neither(self):
        aminos = [['L'], ['L']]
        cond = 'SCORE FROM(MAX ( 1F => 10, 2R => 20 ))|'
        expected_cond = '|'
        expected_score = 0
        expected_mutations = set()

        bnr = self.asi.bnf_scorecondition(cond, aminos)

        self.assertEqual(expected_cond, bnr.cond)
        self.assertEqual(expected_score, bnr.score)
        self.assertEqual(expected_mutations, bnr.mutations)


class AsiAlgorithmNewRulesTest(TestCase):
    default_drugs = """\
  <DRUG>
    <NAME>ABC</NAME>
    <FULLNAME>abacavir</FULLNAME>
    <RULE>
      <CONDITION><![CDATA[SCORE FROM(41L => 15)]]></CONDITION>
      <ACTIONS>
        <SCORERANGE>
          <USE_GLOBALRANGE/>
        </SCORERANGE>
      </ACTIONS>
    </RULE>
  </DRUG>
"""
    default_comments = ""

    @staticmethod
    def create_asi(drugs=default_drugs, comments=default_comments):
        xml = """\
<ALGORITHM>
  <ALGNAME>HIVDB</ALGNAME>
  <ALGVERSION>fake</ALGVERSION>
  <ALGDATE>2017-03-02</ALGDATE>
  <DEFINITIONS>
    <GENE_DEFINITION>
      <NAME>RT</NAME>
      <DRUGCLASSLIST>NRTI</DRUGCLASSLIST>
    </GENE_DEFINITION>
    <LEVEL_DEFINITION>
      <ORDER>1</ORDER>
      <ORIGINAL>Susceptible</ORIGINAL>
      <SIR>S</SIR>
    </LEVEL_DEFINITION>
    <LEVEL_DEFINITION>
      <ORDER>5</ORDER>
      <ORIGINAL>High-Level Resistance</ORIGINAL>
      <SIR>R</SIR>
    </LEVEL_DEFINITION>
    <DRUGCLASS>
      <NAME>NRTI</NAME>
      <DRUGLIST>ABC</DRUGLIST>
    </DRUGCLASS>
    <GLOBALRANGE><![CDATA[(-INF TO 9 => 1,  10 TO INF => 5)]]></GLOBALRANGE>
    <COMMENT_DEFINITIONS>
      {comments}
    </COMMENT_DEFINITIONS>
  </DEFINITIONS>
  {drugs}
  <MUTATION_COMMENTS>
  </MUTATION_COMMENTS>
</ALGORITHM>
""".format(drugs=drugs, comments=comments)
        return AsiAlgorithm(StringIO(xml))

    def test_interpret(self):
        asi = self.create_asi()
        aa_seq = [[amino] for amino in asi.stds['RT']]
        aa_seq[40] = ['L']
        compared_attrs = ('code', 'score', 'level', 'level_name')
        expected_drugs = [('ABC', 15.0, 5, 'High-Level Resistance')]
        expected_mutation_comments = []

        result = asi.interpret(aa_seq, 'RT')

        drugs = list(map(attrgetter(*compared_attrs), result.drugs))
        self.assertEqual(expected_drugs, drugs)
        self.assertEqual(expected_mutation_comments, result.mutation_comments)

    def test_level_action(self):
        drugs = """\
  <DRUG>
    <NAME>ABC</NAME>
    <FULLNAME>abacavir</FULLNAME>
    <RULE>
      <CONDITION><![CDATA[41L]]></CONDITION>
      <ACTIONS>
        <LEVEL>5</LEVEL>
      </ACTIONS>
    </RULE>
  </DRUG>
"""

        asi = self.create_asi(drugs=drugs)
        aa_seq = [[amino] for amino in asi.stds['RT']]
        aa_seq[40] = ['L']
        compared_attrs = ('code', 'score', 'level', 'level_name')
        expected_drugs = [('ABC', 0.0, 5, 'High-Level Resistance')]
        expected_mutation_comments = []

        result = asi.interpret(aa_seq, 'RT')

        drugs = list(map(attrgetter(*compared_attrs), result.drugs))
        self.assertEqual(expected_drugs, drugs)
        self.assertEqual(expected_mutation_comments, result.mutation_comments)

    def test_custom_score_range(self):
        drugs = """\
  <DRUG>
    <NAME>ABC</NAME>
    <FULLNAME>abacavir</FULLNAME>
    <RULE>
      <CONDITION><![CDATA[SCORE FROM(41L => 15)]]></CONDITION>
      <ACTIONS>
        <SCORERANGE>(-INF TO 20 => 1, 20 TO INF => 5)</SCORERANGE>
      </ACTIONS>
    </RULE>
  </DRUG>
"""

        asi = self.create_asi(drugs=drugs)
        aa_seq = [[amino] for amino in asi.stds['RT']]
        aa_seq[40] = ['L']
        compared_attrs = ('code', 'score', 'level')
        expected_drugs = [('ABC', 15.0, 1)]
        expected_mutation_comments = []

        result = asi.interpret(aa_seq, 'RT')

        drugs = list(map(attrgetter(*compared_attrs), result.drugs))
        self.assertEqual(expected_drugs, drugs)
        self.assertEqual(expected_mutation_comments, result.mutation_comments)

    def test_drug_comment(self):
        drugs = """\
  <DRUG>
    <NAME>ABC</NAME>
    <FULLNAME>abacavir</FULLNAME>
    <RULE>
      <CONDITION><![CDATA[41L]]></CONDITION>
      <ACTIONS>
        <COMMENT ref="RT41L"/>
      </ACTIONS>
    </RULE>
  </DRUG>
"""
        comments = """\
      <COMMENT_STRING id="RT41L">
        <TEXT><![CDATA[Example comment.]]></TEXT>
        <SORT_TAG>1</SORT_TAG>
      </COMMENT_STRING>
"""

        asi = self.create_asi(drugs=drugs, comments=comments)
        aa_seq = [[amino] for amino in asi.stds['RT']]
        aa_seq[40] = ['L']
        compared_attrs = ('code', 'score', 'level', 'level_name', 'comments')
        expected_drugs = [('ABC', 0.0, 1, 'Susceptible', ['Example comment.'])]
        expected_mutation_comments = []

        result = asi.interpret(aa_seq, 'RT')

        drugs = list(map(attrgetter(*compared_attrs), result.drugs))
        self.assertEqual(expected_drugs, drugs)
        self.assertEqual(expected_mutation_comments, result.mutation_comments)


class TranslateSequenceTest(TestCase):
    def test_translate(self):
        nucs = 'CCCATTAGT'
        expected_aminos = [['P'], ['I'], ['S']]

        aminos = translate_complete_to_array(nucs)

        self.assertEqual(expected_aminos, aminos)

    def test_translate_mixture(self):
        nucs = 'CCCWTTAGT'
        expected_aminos = [['P'], ['F', 'I'], ['S']]

        aminos = translate_complete_to_array(nucs)

        sorted_aminos = [sorted(choices) for choices in aminos]
        self.assertEqual(expected_aminos, sorted_aminos)

    def test_translate_synonymous_mixture(self):
        nucs = 'CCCTTYAGT'
        expected_aminos = [['P'], ['F'], ['S']]

        aminos = translate_complete_to_array(nucs)

        sorted_aminos = [sorted(choices) for choices in aminos]
        self.assertEqual(expected_aminos, sorted_aminos)
