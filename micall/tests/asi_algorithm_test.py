from operator import attrgetter
from unittest import TestCase

from micall.hivdb.asi_algorithm import AsiAlgorithm, translate_complete_to_array


class AsiAlgorithmTest(TestCase):
    def test_read(self):
        algorithm_hivdb = AsiAlgorithm("micall/hivdb/HIVDB_8.3.xml")
        rt_seq = 'CCCATTAGTCCTATTGAAACTGTACCAGTAAAATTAAAGCCAGGAATGGATGGCCCAAAGGTYAARCAATGGCCATTGACAGAAGAAAAAATAAAAGCATTAGTAGAAATTTGTACAGAAATGGAAAAGGAAGGRAAGATTTCAAAAATTGGACCTGAAAATCCATACAATACTCCAGTATTTGCCATAAAGAAAAAAGACAGTACTAAATGGAGAAAATTAGTAGATTTCAGAGAACTTAATAARAGAACTCAAGACTTCTGGGAAGTTCAATTAGGAATACCACATCCYGCAGGGTTAAAAAAGAAMAAGTCAGTAACAGTACTRGATGTGGGTGATGCATATTTTTCAGTTCCCTTATATGAAGACTTCAGGAAGTATACTGCATTCACCATACCTAGYACAAACAATGAGACACCAGGGATTAGATATCAGTACAATGTGCTGCCACAAGGATGGAAAGGATCACCAGCAATATTCCAAAGTAGCATGATAAAAATCTTAGAGCCTTTCAGAAAACAAAATCCAGARATAGTCATCTATCAATACGTGGATGATTTGTATGTAGSATCTGACTTAGAAATAGGGCAGCATAGAACAAAGATAGAGGAACTGAGAGCACATCTRTTRAAGTGGGGATTTACCACACCAGACAAAAAACATCAGAAAGAGCCTCCATTCCTTTGGATGGGTTATGAACTCCATCCTGATAAATGGACR'
        aa_seq = translate_complete_to_array(rt_seq)
        expected_drugs = [('3TC', 60.0, 1, []),
                          ('ABC', 15.0, 3, []),
                          ('AZT', -10.0, 1, []),
                          ('D4T', -10.0, 1, []),
                          ('DDI', 10.0, 1, []),
                          ('FTC', 60.0, 1, []),
                          ('TDF', -10.0, 1, []),
                          ('EFV', 105.0, 1, []),
                          ('ETR', 10.0, 1, []),
                          ('NVP', 120.0, 1, []),
                          ('RPV', 15.0, 1, [])]
        expected_mutation_comments = [
            'K103N is a non-polymorphic mutation that causes high-level resistance to NVP and EFV.',
            'M184V/I cause high-level in vitro resistance to 3TC and FTC and low-level resistance '
            'to ddI and ABC. However, M184V/I are not contraindications to continued treatment with '
            '3TC or FTC because they increase susceptibility to AZT, TDF and d4T and are associated '
            'with clinically significant reductions in HIV-1 replication.',
            'G190A is a non-polymorphic mutation that causes high-level resistance to NVP and '
            'intermediate resistance to EFV. It has a low weight in the Tibotec ETR genotypic '
            'susceptibility score but does not appear to be selected by ETR or RPV or to reduce '
            'their in vitro susceptibility in the absence of other NNRTI-resistance mutations.']

        result = algorithm_hivdb.interpret(aa_seq, 'RT')
        drugs = list(map(attrgetter('code', 'score', 'level', 'comments'), result.drugs))
        self.assertEqual(expected_drugs, drugs)
        self.assertEqual(expected_mutation_comments, result.mutation_comments)
