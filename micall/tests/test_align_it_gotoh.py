import unittest
from gotoh import align_it, align_it_aa


class AlignItGotohTest(unittest.TestCase):
    def test_empty_ref(self):
        seqa = ""
        seqb = "CTCGGCTTGCTGAAGCGCGCACGGCAAGAGGCGAG" #HIV1B-sl1

        gip=4
        gep=1
        use_terminal_gap_penalty=0

        [result_seqa, result_seqb, result_score] = align_it(seqa, seqb, gip, gep, use_terminal_gap_penalty)

        expected_seqa = "-"*len(seqb)
        expected_seqb = seqb
        expected_score = gip + gep*len(seqb)

        self.assertEqual(expected_seqa, result_seqa)
        self.assertEqual(expected_seqb, result_seqb)
        self.assertEqual(expected_score, result_score)


    def test_empty_seq(self):
        seqa = "CTCGGCTTGCTGAAGCGCGCACGGCAAGAGGCGAG" #HIV1B-sl1
        seqb = ""

        gip=4
        gep=1
        use_terminal_gap_penalty=0

        [result_seqa, result_seqb, result_score] = align_it(seqa, seqb, gip, gep, use_terminal_gap_penalty)

        expected_seqa = seqa
        expected_seqb = "-"*len(seqa)
        expected_score = gip + gep*len(seqa)

        self.assertEqual(expected_seqa, result_seqa)
        self.assertEqual(expected_seqb, result_seqb)
        self.assertEqual(expected_score, result_score)


    def test_empty_all(self):
        seqa = ""
        seqb = ""

        gip=4
        gep=1
        use_terminal_gap_penalty=1

        [result_seqa, result_seqb, result_score] = align_it(seqa, seqb, gip, gep, use_terminal_gap_penalty)

        expected_seqa = ""
        expected_seqb = ""
        expected_score = 0

        self.assertEqual(expected_seqa, result_seqa)
        self.assertEqual(expected_seqb, result_seqb)
        self.assertEqual(expected_score, result_score)


    def test_align_it(self):
        seqa = "CTCGGCTTGCTGAAGCGCGCACGGCAAGAGGCGAG" # HIV1B-sl1
        seqb = "CTAGCGGAGGCTAG" # HIV1B-sl3

        gip=4
        gep=1
        use_terminal_gap_penalty=1

        [result_seqa, result_seqb, result_score] = align_it(seqa, seqb, gip, gep, use_terminal_gap_penalty)

        expected_seqa = "CTCGGCTTGCTGAAGCGCGCACGGCAAGAGGCGAG"
        expected_seqb = "---------CT--AGCG----------GAGGCTAG"
        expected_score = 41

        self.assertEqual(expected_seqa, result_seqa)
        self.assertEqual(expected_seqb, result_seqb)
        self.assertEqual(expected_score, result_score)


    def test_align_it_zero_penalty(self):
        seqa = "CTCGGCTTGCTGAAGCGCGCACGGCAAGAGGCGAG" # HIV1B-sl1
        seqb = "CTAGCGGAGGCTAG" # HIV1B-sl3

        gip=0
        gep=0
        use_terminal_gap_penalty=1

        [result_seqa, result_seqb, result_score] = align_it(seqa, seqb, gip, gep, use_terminal_gap_penalty)

        expected_seqa = "CTCGGCTTGCTGAAGCGCGCACGGC-AAGAGGCGAG"
        expected_seqb = "CT----------A-GCG-G-A-GGCTA-G-------"
        expected_score = 65

        self.assertEqual(expected_seqa, result_seqa)
        self.assertEqual(expected_seqb, result_seqb)
        self.assertEqual(expected_score, result_score)


    def test_align_it_aa(self):
        # HIV1B-tat
        seqa = "MEPVDPRLEPWKHPGSQPKTACTNCYCKKCCFHCQVCFITKALGISYGRKKRRQRRRAHQNSQTHQASLSKQ*"
        # HIV1B-vpr
        seqb = "MEQAPEDQGPQREPHNEWTLELLEELKNEAVRHFPRIWLHGLGQHIYETYGDTWAGVEAIIRILQQLLFIHFRIGCRHSRIGVTRQRRARNGASRS*"

        gip=4
        gep=1
        use_terminal_gap_penalty=1

        [result_seqa, result_seqb, result_score] = align_it_aa(seqa, seqb, gip, gep, use_terminal_gap_penalty)

        expected_seqa = "ME--PVD--P-RLEP---W--------K-----H-P-----G-SQP--KTACTNCY---C----------KKCC--F-HCQVCF-ITKALG-----ISYG--RKKRRQRRRAHQNSQTHQASLSKQ*"
        expected_seqb = "MEQAPEDQGPQR-EPHNEWTLELLEELKNEAVRHFPRIWLHGLGQHIYET-----YGDTWAGVEAIIRILQQ--LLFIH----FRI----GCRHSRI--GVT----RQRR-AR-NG----ASRS--*"
        expected_score = 49

        self.assertEqual(expected_seqa, result_seqa)
        self.assertEqual(expected_seqb, result_seqb)
        self.assertEqual(expected_score, result_score)


    def test_align_it_aa_zero_penalty(self):
        # HIV1B-tat
        seqa = "MEPVDPRLEPWKHPGSQPKTACTNCYCKKCCFHCQVCFITKALGISYGRKKRRQRRRAHQNSQTHQASLSKQ*"
        # HIV1B-vpr
        seqb = "MEQAPEDQGPQREPHNEWTLELLEELKNEAVRHFPRIWLHGLGQHIYETYGDTWAGVEAIIRILQQLLFIHFRIGCRHSRIGVTRQRRARNGASRS*"

        gip=0
        gep=0
        use_terminal_gap_penalty=1

        [result_seqa, result_seqb, result_score] = align_it_aa(seqa, seqb, gip, gep, use_terminal_gap_penalty)

        expected_seqa = "ME--PV-D--P-RLEP---W--------K-----H-P-----GS--QPK----TACTNCYC-----------KKCC------F-HCQVCFITKALGISYG-RKK--R----RQRRRAHQ-NSQTHQ-ASL-SKQ*"
        expected_seqb = "MEQAP-EDQGPQR-EPHNEWTLELLEELKNEAVRHFPRIWLHG-LGQ--HIYET-----Y-GDTWAGVEAIIR---ILQQLLFIH----F--R---I--GCR--HSRIGVTRQRR-A--RN-----GAS-RS--*"
        expected_score = 260

        self.assertEqual(expected_seqa, result_seqa)
        self.assertEqual(expected_seqb, result_seqb)
        self.assertEqual(expected_score, result_score)