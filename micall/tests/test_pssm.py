from unittest.case import TestCase

from micall.g2p.pssm_lib import Pssm


class PssmTest(TestCase):
    def test_single_sequence(self):
        pssm = Pssm()
        nucs = ('TGTACAAGACCCAACAACAATACAAGAAAAAGTATACATATAGGACCAGGGAG'
                'AGCATTTTATGCAACAGGAGAAATAATAGGAGATATAAGACAAGCACATTGT')
        expected_aa = [['C'], ['T'], ['R'], ['P'], ['N'], ['-'], ['N'], ['N'],
                       ['T'], ['-'], ['-'], ['R'], ['K'], ['S'], ['I'], ['H'],
                       ['I'], ['-'], ['-'], ['-'], ['G'], ['P'], ['G'], ['R'],
                       ['-'], ['-'], ['-'], ['A'], ['F'], ['Y'], ['A'], ['T'],
                       ['-'], ['-'], ['-'], ['-'], ['G'], ['E'], ['I'], ['I'],
                       ['G'], ['D'], ['I'], ['-'], ['-'], ['R'], ['Q'], ['A'],
                       ['H'], ['C']]
        expected_score = 0.067753

        score, aligned_aa = pssm.run_g2p(nucs)

        self.assertEqual(expected_aa, aligned_aa)
        self.assertAlmostEqual(expected_score, score, places=5)

    def test_multiple_sequences(self):
        pssm = Pssm()
        nucs = [('TGTACAAGACCCAACAACAATACAAGAAAAAGTATACATATAGGACCAGGGAG'
                 'AGCATTTTATGCAACAGGAGAAATAATAGGAGATATAAGACAAGCACATTGT'),
                ('TGTACAAGTCCCAACAACAATACAAGAAAAAGTATACATATAGGACCAGGGAG'
                 'AGCATTTTATGCAACAGGAGAAATAATAGGAGATATAAGACAAGCACATTGT')]
        expected_aa = None  # Not returned when submitting more than one seq.
        expected_scores = [0.06775, 0.06486]

        scores, aligned_aa = pssm.run_g2p(nucs)

        rounded_scores = [round(score, 5) for score in scores]
        self.assertEqual(expected_aa, aligned_aa)
        self.assertEqual(expected_scores, rounded_scores)

    def test_ambiguous_sequence(self):
        pssm = Pssm()
        nucs = ('TGTACAAGWCCCAACAACAATACAAGAAAAAGTATACATATAGGACCAGGGAG'
                'AGCATTTTATGCAACAGGAGAAATAATAGGAGATATAAGACAAGCACATTGT')
        expected_aa = [['C'], ['T'], ['R', 'S'], ['P'], ['N'], ['-'], ['N'], ['N'],
                       ['T'], ['-'], ['-'], ['R'], ['K'], ['S'], ['I'], ['H'],
                       ['I'], ['-'], ['-'], ['-'], ['G'], ['P'], ['G'], ['R'],
                       ['-'], ['-'], ['-'], ['A'], ['F'], ['Y'], ['A'], ['T'],
                       ['-'], ['-'], ['-'], ['-'], ['G'], ['E'], ['I'], ['I'],
                       ['G'], ['D'], ['I'], ['-'], ['-'], ['R'], ['Q'], ['A'],
                       ['H'], ['C']]
        # Average of two possible scores (see test_multiple_sequences).
        expected_score = (0.06775 + 0.06486) / 2

        score, aligned_aa = pssm.run_g2p(nucs)

        self.assertEqual(expected_aa, aligned_aa)
        self.assertAlmostEqual(expected_score, score, places=5)
