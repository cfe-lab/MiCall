import unittest
from micall.utils.translation import translate, reverse_and_complement


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

    def testListAmbiguousOverridesMixturesNotTranslated(self):
        nucs = 'TTY'
        expected_aminos = 'F'

        aminos = translate(nucs, translate_mixtures=False, list_ambiguous=True)

        self.assertEqual(expected_aminos, aminos)

    def testAmbiguousAminosListed(self):
        nucs = 'TTM'  # TTA or TTC: map to L or F
        expected_aminos = '[FL]'

        aminos = translate(nucs, list_ambiguous=True)

        self.assertEqual(expected_aminos, aminos)

    def testReturnList(self):
        nucs = 'CGATTM'  # TTA or TTC: map to L or F
        expected_aminos = [['R'], ['F', 'L']]

        aminos = translate(nucs, return_list=True)

        self.assertEqual(expected_aminos, aminos)

    def testReturnListWithoutMixtures(self):
        """ Don't know why you would use this combination, but stay sane. """

        nucs = 'CGATTM'  # TTA or TTC: map to L or F
        expected_aminos = [['R'], ['?']]

        aminos = translate(nucs, return_list=True, translate_mixtures=False)

        self.assertEqual(expected_aminos, aminos)

    def testStatisticsUnambiguous(self):
        nucs = 'TTATTCTTTTTA'
        expected_aminos = 'LFFL'
        stats = {}
        expected_stats = dict(length=4, ambiguous=0, max_aminos=1)

        aminos = translate(nucs, stats=stats, list_ambiguous=True)

        self.assertEqual(expected_aminos, aminos)
        self.assertEqual(expected_stats, stats)

    def testStatisticsBlank(self):
        nucs = ''
        expected_aminos = ''
        stats = {}
        expected_stats = dict(length=0, ambiguous=0, max_aminos=0)

        aminos = translate(nucs, stats=stats, list_ambiguous=True)

        self.assertEqual(expected_aminos, aminos)
        self.assertEqual(expected_stats, stats)

    def testStatisticsAmbiguous(self):
        nucs = 'TTMTTCNTTTTA'
        expected_aminos = '[FL]F[FILV]L'
        stats = {}
        expected_stats = dict(length=4, ambiguous=2, max_aminos=4)

        aminos = translate(nucs, stats=stats, list_ambiguous=True)

        self.assertEqual(expected_aminos, aminos)
        self.assertEqual(expected_stats, stats)


class ReverseAndComplementTest(unittest.TestCase):
    def testSimple(self):
        fwd = 'ACTG'
        expected = 'CAGT'

        rev = reverse_and_complement(fwd)

        self.assertEqual(expected, rev)
