from operator import attrgetter
from unittest import TestCase

from micall.hivdb.asi_algorithm import AsiAlgorithm, translate_complete_to_array


class AsiAlgorithmTest(TestCase):
    def test_interpret(self):
        algorithm_hivdb = AsiAlgorithm("micall/hivdb/HIVDB_8.3.xml")
        aa_seq = ['A'] * 40 + ['L'] + ['A']
        compared_attrs = ('code', 'score', 'level', 'comments')
        expected_drugs = [('3TC', 0.0, 1, []),
                          ('ABC', 5.0, 1, []),
                          ('AZT', 15.0, 3, []),
                          ('D4T', 15.0, 1, []),
                          ('DDI', 10.0, 1, []),
                          ('FTC', 0.0, 1, []),
                          ('TDF', 5.0, 1, []),
                          ('EFV', 0.0, 1, []),
                          ('ETR', 0.0, 1, []),
                          ('NVP', 0.0, 1, []),
                          ('RPV', 0.0, 1, [])]
        expected_mutation_comments = [
            'M41L is a TAM that usually occurs with T215Y. In combination, '
            'M41L plus T215Y confer intermediate / high-level resistance to '
            'AZT and d4T and contribute to reduced ddI, ABC and TDF '
            'susceptibility.']

        result = algorithm_hivdb.interpret(aa_seq, 'RT')
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
        expected_aminos = [['P'], ['I', 'F'], ['S']]

        aminos = translate_complete_to_array(nucs)

        self.assertEqual(expected_aminos, aminos)
