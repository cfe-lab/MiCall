from operator import attrgetter
from unittest import TestCase

from micall.hivdb.asi_algorithm import AsiAlgorithm, translate_complete_to_array


class AsiAlgorithmTest(TestCase):
    def test_interpret(self):
        algorithm_hivdb = AsiAlgorithm("micall/hivdb/HIVDB_8.3.xml")
        aa_seq = [[amino] for amino in algorithm_hivdb.rt_std]
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

        result = algorithm_hivdb.interpret(aa_seq, 'RT')

        drugs = list(map(attrgetter(*compared_attrs), result.drugs))
        self.assertEqual(expected_drugs, drugs)
        self.assertEqual(expected_mutation_comments, result.mutation_comments)

    def test_protease(self):
        self.maxDiff = None
        algorithm_hivdb = AsiAlgorithm("micall/hivdb/HIVDB_8.3.xml")
        aa_seq = [[amino] for amino in algorithm_hivdb.pr_std]
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

        result = algorithm_hivdb.interpret(aa_seq, 'PR')

        drugs = list(map(attrgetter(*compared_attrs), result.drugs))
        self.assertEqual(expected_drugs, drugs)
        self.assertEqual(expected_mutation_comments, result.mutation_comments)

    def test_gaps(self):
        self.maxDiff = None
        algorithm_hivdb = AsiAlgorithm("micall/hivdb/HIVDB_8.3.xml")
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

        result = algorithm_hivdb.interpret(aa_seq, 'PR')

        drugs = list(map(attrgetter(*compared_attrs), result.drugs))
        self.assertEqual(expected_drugs, drugs)
        self.assertEqual(expected_mutation_comments, result.mutation_comments)

    # def test_mutations(self):
    #     algorithm_hivdb = AsiAlgorithm("micall/hivdb/HIVDB_8.3.xml")
    #     aa_seq = [[amino] for amino in algorithm_hivdb.rt_std]
    #     self.assertNotEqual('L', aa_seq[40])
    #     aa_seq[40] = ['L']
    #     self.assertEqual('R', aa_seq[41])
    #     aa_seq[41] = ['L']
    #     expected_mutations = ['RT41L']
    #
    #     result = algorithm_hivdb.interpret(aa_seq, 'RT')
    #
    #     self.assertEqual(expected_mutations, result.mutations)


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
