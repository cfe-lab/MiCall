import os
from io import StringIO
from operator import attrgetter
from unittest import TestCase

from micall.hivdb.asi_algorithm import AsiAlgorithm


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
        aa_seq = [[]] * 23 + [['I']] + [[]] * 75
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
        aa_seq1 = [[]] * 40 + [['L']] + [[]] * 399
        aa_seq2 = [[]] * 23 + [['I']] + [[]] * 75
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
        aa_seq = [[]] * 440
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

    def test_get_gene_positions(self):
        drugs = """\
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
        expected_positions = {41}
        asi = self.create_asi(drugs=drugs)

        positions = asi.get_gene_positions('RT')

        self.assertEqual(expected_positions, positions)

    def test_score_list_positions(self):
        drugs = """\
          <DRUG>
            <NAME>ABC</NAME>
            <FULLNAME>abacavir</FULLNAME>
            <RULE>
              <CONDITION><![CDATA[SCORE FROM(41L => 15, 99F => 30)]]></CONDITION>
              <ACTIONS>
                <SCORERANGE>
                  <USE_GLOBALRANGE/>
                </SCORERANGE>
              </ACTIONS>
            </RULE>
          </DRUG>
        """
        expected_positions = {41, 99}
        asi = self.create_asi(drugs=drugs)

        positions = asi.get_gene_positions('RT')

        self.assertEqual(expected_positions, positions)

    def test_and_positions(self):
        drugs = """\
          <DRUG>
            <NAME>ABC</NAME>
            <FULLNAME>abacavir</FULLNAME>
            <RULE>
              <CONDITION><![CDATA[SCORE FROM((41L AND 99F) => 15)]]></CONDITION>
              <ACTIONS>
                <SCORERANGE>
                  <USE_GLOBALRANGE/>
                </SCORERANGE>
              </ACTIONS>
            </RULE>
          </DRUG>
        """
        expected_positions = {41, 99}
        asi = self.create_asi(drugs=drugs)

        positions = asi.get_gene_positions('RT')

        self.assertEqual(expected_positions, positions)

    def test_or_position(self):
        drugs = """\
          <DRUG>
            <NAME>ABC</NAME>
            <FULLNAME>abacavir</FULLNAME>
            <RULE>
              <CONDITION><![CDATA[SCORE FROM((41L OR 99F) => 15)]]></CONDITION>
              <ACTIONS>
                <SCORERANGE>
                  <USE_GLOBALRANGE/>
                </SCORERANGE>
              </ACTIONS>
            </RULE>
          </DRUG>
        """
        expected_positions = {41, 99}
        asi = self.create_asi(drugs=drugs)

        positions = asi.get_gene_positions('RT')

        self.assertEqual(expected_positions, positions)

    def test_max_position(self):
        drugs = """\
          <DRUG>
            <NAME>ABC</NAME>
            <FULLNAME>abacavir</FULLNAME>
            <RULE>
              <CONDITION><![CDATA[SCORE FROM(MAX (41A => 30, 99C => 60))]]></CONDITION>
              <ACTIONS>
                <SCORERANGE>
                  <USE_GLOBALRANGE/>
                </SCORERANGE>
              </ACTIONS>
            </RULE>
          </DRUG>
        """
        expected_positions = {41, 99}
        asi = self.create_asi(drugs=drugs)

        positions = asi.get_gene_positions('RT')

        self.assertEqual(expected_positions, positions)

    def test_negative_positions(self):
        drugs = """\
          <DRUG>
            <NAME>ABC</NAME>
            <FULLNAME>abacavir</FULLNAME>
            <RULE>
              <CONDITION><![CDATA[SCORE FROM(41A => -10)]]></CONDITION>
              <ACTIONS>
                <SCORERANGE>
                  <USE_GLOBALRANGE/>
                </SCORERANGE>
              </ACTIONS>
            </RULE>
          </DRUG>
        """
        expected_positions = {41}
        asi = self.create_asi(drugs=drugs)

        positions = asi.get_gene_positions('RT')

        self.assertEqual(expected_positions, positions)

    def test(self):
        drugs = """\
          <DRUG>
            <NAME>ABC</NAME>
            <FULLNAME>abacavir</FULLNAME>
            <RULE>
              <CONDITION><![CDATA[SCORE FROM(MAX((41A AND 99F) => 5, 62K => 10))]]></CONDITION>
              <ACTIONS>
                <SCORERANGE>
                  <USE_GLOBALRANGE/>
                </SCORERANGE>
              </ACTIONS>
            </RULE>
          </DRUG>
        """
        expected_positions = {41, 99, 62}
        asi = self.create_asi(drugs=drugs)

        positions = asi.get_gene_positions('RT')

        self.assertEqual(expected_positions, positions)

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


class AsiAlgorithmJsonRulesTest(TestCase):
    default_drugs = [{"name": "madeupivir",
                      "genotypes": [{"rules": "SCORE FROM(41L => 4)",
                                     "genotype": "1A",
                                     "region": "NS5a",
                                     "reference": "HCV1A-H77-NS5a"},
                                    {"rules": "SCORE FROM(42L => 8)",
                                     "genotype": "1B",
                                     "region": "NS5a",
                                     "reference": "HCV1B-Con1-NS5a"}],
                      "code": "MDP"}]

    def test_interpret(self):
        asi = AsiAlgorithm(rules_yaml=self.default_drugs, genotype='1A')
        aa_seq = [['A']] * 40 + [['L']] + [['A']] * 407
        compared_attrs = ('code', 'score', 'level', 'level_name')
        expected_drugs = [('MDP', 4.0, 2, 'Resistance Possible')]
        expected_mutation_comments = []

        result = asi.interpret(aa_seq, 'HCV1A-H77-NS5a')
        drugs = list(map(attrgetter(*compared_attrs), result.drugs))
        self.assertEqual(expected_drugs, drugs)
        self.assertEqual(expected_mutation_comments, result.mutation_comments)
