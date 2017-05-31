import os
from operator import attrgetter
from unittest import TestCase

from micall.hivdb.asi_algorithm import AsiAlgorithm, translate_complete_to_array


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
        expected_mutations = ['M41L']

        result = self.asi.interpret(aa_seq, 'RT')

        self.assertEqual(expected_mutations, result.mutations)

    def test_bnf_residue(self):
        aminos = [['L'], ['R']]
        cond = '2R'
        expected_cond = ''
        expected_score = 0
        expected_mutations = {'2R'}

        bnr = self.asi.bnf_residue(cond, aminos)

        self.assertEqual(expected_cond, bnr.cond)
        self.assertTrue(bnr.truth)
        self.assertEqual(expected_score, bnr.score)
        self.assertEqual(expected_mutations, bnr.mutations)

    def test_bnf_residue_multiple_choices(self):
        aminos = [['L'], ['R']]
        cond = '2RFL'
        expected_cond = ''
        expected_score = 0
        expected_mutations = {'2R'}

        bnr = self.asi.bnf_residue(cond, aminos)

        self.assertEqual(expected_cond, bnr.cond)
        self.assertTrue(bnr.truth)
        self.assertEqual(expected_score, bnr.score)
        self.assertEqual(expected_mutations, bnr.mutations)

    def test_bnf_residue_multiple_choices_found(self):
        aminos = [['L'], ['R', 'F', 'K']]
        cond = '2RFL'
        expected_cond = ''
        expected_score = 0
        expected_mutations = {'2R', '2F'}

        bnr = self.asi.bnf_residue(cond, aminos)

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
