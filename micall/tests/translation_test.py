import unittest
from micall.core import aln2counts


class TranslateTest(unittest.TestCase):
    def testSingleCodon(self):
        nucs = 'TTT'
        expected_aminos = 'F'

        aminos = aln2counts.translate(nucs)

        self.assertEqual(expected_aminos, aminos)

    def testPartialCodon(self):
        nucs = 'TTTCC'
        expected_aminos = 'F'

        aminos = aln2counts.translate(nucs)

        self.assertEqual(expected_aminos, aminos)

    def testTwoCodons(self):
        nucs = 'TTTCCT'
        expected_aminos = 'FP'

        aminos = aln2counts.translate(nucs)

        self.assertEqual(expected_aminos, aminos)

    def testLowerCase(self):
        nucs = 'TttCCT'
        expected_aminos = 'FP'

        aminos = aln2counts.translate(nucs)

        self.assertEqual(expected_aminos, aminos)

    def testOffset(self):
        nucs = "TTTCCT"
        offset = 3
        expected_aminos = "-FP"

        aminos = aln2counts.translate(nucs, offset)

        self.assertEqual(expected_aminos, aminos)

    def testSingleDashAmbiguous(self):
        nucs = '-TT'
        expected_aminos = '?'

        aminos = aln2counts.translate(nucs)

        self.assertEqual(expected_aminos, aminos)

    def testSingleDashUnambiguous(self):
        nucs = 'CG-'  # CGA, CGC, CGG, CGT all map to R
        expected_aminos = 'R'

        aminos = aln2counts.translate(nucs)

        self.assertEqual(expected_aminos, aminos)

    def testTwoDashes(self):
        nucs = '--T'
        expected_aminos = '?'

        aminos = aln2counts.translate(nucs)

        self.assertEqual(expected_aminos, aminos)

    def testThreeDashes(self):
        nucs = '---'
        expected_aminos = '-'

        aminos = aln2counts.translate(nucs)

        self.assertEqual(expected_aminos, aminos)

    def testAmbiguousBasesThatAreSynonyms(self):
        nucs = 'TTY'  # TTC or TTT: both map to F
        expected_aminos = 'F'

        aminos = aln2counts.translate(nucs)

        self.assertEqual(expected_aminos, aminos)

    def testTwoAmbiguousBasesThatAreSynonyms(self):
        nucs = 'MGR'  # CGA, CGG, AGA, or AGG: all map to R
        expected_aminos = 'R'

        aminos = aln2counts.translate(nucs)

        self.assertEqual(expected_aminos, aminos)

    def testTwoAmbiguousBasesThatAreNotSynonyms(self):
        nucs = 'RGR'  # GGA, GGG, AGA, or AGG: map to G and R, respectively
        expected_aminos = '?'

        aminos = aln2counts.translate(nucs)

        self.assertEqual(expected_aminos, aminos)

    def testMixturesNotTranslated(self):
        nucs = 'TTY'  # TTC or TTT: both map to F
        expected_aminos = '?'

        aminos = aln2counts.translate(nucs, translate_mixtures=False)

        self.assertEqual(expected_aminos, aminos)
