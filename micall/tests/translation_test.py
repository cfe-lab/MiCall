import unittest
from micall.utils.translation import translate


class TranslateTest(unittest.TestCase):
    def testSingleCodon(self):
        nucs = 'TTT'
        expected_aminos = 'F'

        aminos = translate(nucs)

        self.assertEqual(expected_aminos, aminos)

    def testPartialCodon(self):
        nucs = 'TTTCC'
        expected_aminos = 'F'

        aminos = translate(nucs)

        self.assertEqual(expected_aminos, aminos)

    def testTwoCodons(self):
        nucs = 'TTTCCT'
        expected_aminos = 'FP'

        aminos = translate(nucs)

        self.assertEqual(expected_aminos, aminos)

    def testLowerCase(self):
        nucs = 'TttCCT'
        expected_aminos = 'FP'

        aminos = translate(nucs)

        self.assertEqual(expected_aminos, aminos)

    def testOffset(self):
        nucs = "TTTCCT"
        offset = 3
        expected_aminos = "-FP"

        aminos = translate(nucs, offset)

        self.assertEqual(expected_aminos, aminos)

    def testSingleDashAmbiguous(self):
        nucs = '-TT'
        expected_aminos = '?'

        aminos = translate(nucs)

        self.assertEqual(expected_aminos, aminos)

    def testSingleDashUnambiguous(self):
        nucs = 'CG-'  # CGA, CGC, CGG, CGT all map to R
        expected_aminos = 'R'

        aminos = translate(nucs)

        self.assertEqual(expected_aminos, aminos)

    def testTwoDashes(self):
        nucs = '--T'
        expected_aminos = '?'

        aminos = translate(nucs)

        self.assertEqual(expected_aminos, aminos)

    def testThreeDashes(self):
        nucs = '---'
        expected_aminos = '-'

        aminos = translate(nucs)

        self.assertEqual(expected_aminos, aminos)

    def testAmbiguousBasesThatAreSynonyms(self):
        nucs = 'TTY'  # TTC or TTT: both map to F
        expected_aminos = 'F'

        aminos = translate(nucs)

        self.assertEqual(expected_aminos, aminos)

    def testTwoAmbiguousBasesThatAreSynonyms(self):
        nucs = 'MGR'  # CGA, CGG, AGA, or AGG: all map to R
        expected_aminos = 'R'

        aminos = translate(nucs)

        self.assertEqual(expected_aminos, aminos)

    def testTwoAmbiguousBasesThatAreNotSynonyms(self):
        nucs = 'RGR'  # GGA, GGG, AGA, or AGG: map to G and R, respectively
        expected_aminos = '?'

        aminos = translate(nucs)

        self.assertEqual(expected_aminos, aminos)

    def testMixturesNotTranslated(self):
        nucs = 'TTY'  # TTC or TTT: both map to F
        expected_aminos = '?'

        aminos = translate(nucs, translate_mixtures=False)

        self.assertEqual(expected_aminos, aminos)
