from unittest.mock import patch, call
from unittest import TestCase
from io import StringIO
from micall.core.coverage_plots import coverage_plot, plot_concordance
from micall.core.project_config import ProjectConfig


class CoveragePlotsTest(TestCase):
    def setUp(self):
        self.addTypeEqualityFunc(str, self.assertMultiLineEqual)
        config_json = StringIO("""\
{
  "projects": {
    "R1": {
      "max_variants": 0,
      "regions": [
        {
          "coordinate_region": "R1",
          "coordinate_region_length": 3,
          "key_positions": [],
          "min_coverage1": 10,
          "min_coverage2": 50,
          "min_coverage3": 100,
          "seed_region_names": [
            "R1-seed"
          ]
        }
      ]
    },
    "R1-and-R2": {
      "max_variants": 0,
      "regions": [
        {
          "coordinate_region": "R1",
          "coordinate_region_length": 3,
          "key_positions": [
            {
              "end_pos": null,
              "start_pos": 1
            },
            {
              "end_pos": null,
              "start_pos": 3
            }
          ],
          "min_coverage1": 10,
          "min_coverage2": 50,
          "min_coverage3": 100,
          "seed_region_names": [
            "R1-seed"
          ]
        },
        {
          "coordinate_region": "R2",
          "coordinate_region_length": 1,
          "key_positions": [],
          "min_coverage1": 10,
          "min_coverage2": 50,
          "min_coverage3": 100,
          "seed_region_names": [
            "R2-seed"
          ]
        }
      ]
    }
  }
}
""")
        self.config = ProjectConfig()
        self.config.load(config_json)

    @patch('matplotlib.pyplot.savefig')
    @patch('micall.core.project_config.ProjectConfig.loadScoring')
    def test_simple(self, config_mock, savefig_mock):
        config_mock.return_value = self.config
        amino_csv = StringIO("""\
seed,region,q-cutoff,query.aa.pos,refseq.aa.pos,\
A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,*,X,partial,del,ins,clip,v3_overlap
R1-seed,R1,15,100,1,0,5,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
R1-seed,R1,15,101,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,5,0,0,0,0,0,0,0,0,0,0
R1-seed,R1,15,102,3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,5,0,0,0,0,0,0,0,0,0,0,0,0
""")
        expected_scores = """\
project,region,seed,q.cut,min.coverage,which.key.pos,off.score,on.score
R1,R1,R1-seed,15,5,1,-1,1
R1-and-R2,R1,R1-seed,15,5,1,-1,1
"""
        scores_csv = StringIO()
        amino_csv.name = 'E1234.amino.csv'
        expected_calls = [call('E1234.R1.R1.png'),
                          call('E1234.R1.R1.details.png'),
                          call('E1234.R1-and-R2.R1.png'),
                          call('E1234.R1-and-R2.R1.details.png')]

        coverage_plot(amino_csv,
                      coverage_scores_csv=scores_csv,
                      coverage_maps_prefix='E1234')

        self.assertEqual(expected_calls, savefig_mock.mock_calls)
        self.assertEqual(expected_scores, scores_csv.getvalue())

    @patch('matplotlib.pyplot.savefig')
    def test_concordance_plot(self, savefig_mock):
        concordance_detailed_csv = StringIO("""\
reference,region,concordance,coverage,position
R1A-seed,R1A,0,0,10.0
R1A-seed,R1A_second,1.0,1.0,10.0
R1A-seed,R1A_second,1.0,1.0,11.0
R1A-seed,R1A_second,0.95,0.95,12.0
R1A-seed,R1A_second,0.9,0.9,13.0
""")

        expected_calls = [call('test/path/concordance.R1A-seed.R1A.png'),
                          call('test/path/concordance.R1A-seed.R1A_second.png')]

        plot_concordance(concordance_detailed_csv,
                         plot_path="test/path",)

        self.assertEqual(expected_calls, savefig_mock.mock_calls)
