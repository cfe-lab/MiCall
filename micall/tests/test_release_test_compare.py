from collections import defaultdict
from io import StringIO
from unittest import TestCase, skip

from release_test_compare import compare_sample, SampleFiles, Sample, \
    MiseqRun, Scenarios, ConsensusDistance, group_samples_file, \
    group_nucs_file, compare_consensus, map_consensus_sequences


class CompareSampleTest(TestCase):
    def test_empty(self):
        sample = Sample(MiseqRun(target_path='run1/Results/versionX'),
                        'sample42',
                        SampleFiles(),
                        SampleFiles())
        expected_report = ''

        report, _, _ = compare_sample(sample)

        self.assertEqual(expected_report, report)

    def test_x4_big_change(self):
        sample = Sample(MiseqRun(target_path='run1/Results/versionX'),
                        'sample42',
                        SampleFiles(g2p_summary=[dict(X4pct='50.00')]),
                        SampleFiles(g2p_summary=[dict(X4pct='60.00')]))
        expected_report = 'run1:sample42 G2P: 50.00 => 60.00\n'

        report, _, _ = compare_sample(sample)

        self.assertEqual(expected_report, report)

    def test_x4_small_change(self):
        sample = Sample(MiseqRun(target_path='run1/Results/versionX'),
                        'sample42',
                        SampleFiles(g2p_summary=[dict(X4pct='50.00')]),
                        SampleFiles(g2p_summary=[dict(X4pct='51.00')]))
        expected_report = ''

        report, _, _ = compare_sample(sample)

        self.assertEqual(expected_report, report)

    def test_blank(self):
        sample = Sample(MiseqRun(target_path='run1/Results/versionX'),
                        'sample42',
                        SampleFiles(g2p_summary=[dict(X4pct='50.00')]),
                        SampleFiles(g2p_summary=[dict(X4pct='')]))
        expected_report = 'run1:sample42 G2P: 50.00 => \n'

        report, _, _ = compare_sample(sample)

        self.assertEqual(expected_report, report)

    def test_other_difference_with_blanks(self):
        sample = Sample(MiseqRun(target_path='run1/Results/versionX'),
                        'sample42',
                        SampleFiles(g2p_summary=[dict(X4pct='', other='x')]),
                        SampleFiles(g2p_summary=[dict(X4pct='', other='y')]))
        expected_report = ''

        report, _, _ = compare_sample(sample)

        self.assertEqual(expected_report, report)

    def test_same_final(self):
        sample = Sample(MiseqRun(target_path='run1/Results/versionX'),
                        'sample42',
                        SampleFiles(g2p_summary=[dict(X4pct='2.99', final='X4')]),
                        SampleFiles(g2p_summary=[dict(X4pct='3.01', final='X4')]))
        expected_report = ''

        report, _, _ = compare_sample(sample)

        self.assertEqual(expected_report, report)

    def test_different_final(self):
        sample = Sample(MiseqRun(target_path='run1/Results/versionX'),
                        'sample42',
                        SampleFiles(g2p_summary=[dict(X4pct='1.99', final='R5')]),
                        SampleFiles(g2p_summary=[dict(X4pct='2.01', final='X4')]))
        expected_report = 'run1:sample42 G2P: R5 1.99 => X4 2.01\n'

        report, _, _ = compare_sample(sample)

        self.assertEqual(expected_report, report)

    def test_same_coverage(self):
        sample = Sample(MiseqRun(target_path='run1/Results/versionX'),
                        'sample42',
                        SampleFiles(coverage_scores=[{'project': 'HIV',
                                                      'region': 'PR',
                                                      'on.score': '4'}]),
                        SampleFiles(coverage_scores=[{'project': 'HIV',
                                                      'region': 'PR',
                                                      'on.score': '4'}]))
        expected_report = ''

        report, _, _ = compare_sample(sample)

        self.assertEqual(expected_report, report)

    def test_single_coverage(self):
        sample = Sample(MiseqRun(target_path='run1/Results/versionX'),
                        'sample42',
                        SampleFiles(coverage_scores=[{'project': 'HIV',
                                                      'region': 'PR',
                                                      'on.score': '4'}]),
                        SampleFiles(coverage_scores=[{'project': 'HIV',
                                                      'region': 'PR',
                                                      'on.score': '3'}]))
        expected_report = 'run1:sample42 coverage: HIV PR 4 => 3\n'

        report, _, _ = compare_sample(sample)

        self.assertEqual(expected_report, report)

    def test_multiple_coverage(self):
        sample = Sample(MiseqRun(target_path='run1/Results/versionX'),
                        'sample42',
                        SampleFiles(coverage_scores=[{'project': 'HIV',
                                                      'region': 'PR',
                                                      'on.score': '3'},
                                                     {'project': 'HIV',
                                                      'region': 'RT',
                                                      'on.score': '2'}]),
                        SampleFiles(coverage_scores=[{'project': 'HIV',
                                                      'region': 'RT',
                                                      'on.score': '3'},
                                                     {'project': 'HIV',
                                                      'region': 'PR',
                                                      'on.score': '4'}]))
        expected_report = ('run1:sample42 coverage: HIV PR 3 => 4\n'
                           'run1:sample42 coverage: HIV RT 2 => 3\n')

        report, _, _ = compare_sample(sample)

        self.assertEqual(expected_report, report)

    def test_missing_coverage(self):
        sample = Sample(MiseqRun(target_path='run1/Results/versionX'),
                        'sample42',
                        SampleFiles(coverage_scores=None),
                        SampleFiles(coverage_scores=[{'project': 'HIV',
                                                      'region': 'PR',
                                                      'on.score': '2'}]))
        expected_report = 'run1:sample42 coverage: HIV PR - => 2\n'

        report, _, _ = compare_sample(sample)

        self.assertEqual(expected_report, report)

    def test_missing_coverage_same_remap_counts(self):
        sample = Sample(MiseqRun(target_path='run1/Results/versionX'),
                        'sample42',
                        SampleFiles(coverage_scores=[{'project': 'HCV',
                                                      'region': 'HCV1-E1',
                                                      'seed': 'HCV1-seed',
                                                      'on.score': '2'},
                                                     {'project': 'HCV',
                                                      'region': 'HCV2-E1',
                                                      'seed': 'HCV2-seed',
                                                      'on.score': '3'}],
                                    remap_counts=[{'type': 'remap-1 HCV1-seed'},
                                                  {'type': 'remap-1 HCV2-seed'}]),
                        SampleFiles(coverage_scores=[{'project': 'HCV',
                                                      'region': 'HCV1-E1',
                                                      'seed': 'HCV1-seed',
                                                      'on.score': '2'}],
                                    remap_counts=[{'type': 'remap-1 HCV1-seed'},
                                                  {'type': 'remap-1 HCV2-seed'}]))
        expected_report = 'run1:sample42 coverage: HCV HCV2-E1 3 => -\n'
        expected_scenario_counts = {}

        report, scenario_counts, _ = compare_sample(sample)

        self.assertEqual(expected_report, report)
        self.assertEqual(expected_scenario_counts, scenario_counts)

    def test_missing_coverage_different_remap_counts(self):
        sample = Sample(MiseqRun(target_path='run1/Results/versionX'),
                        'sample42',
                        SampleFiles(coverage_scores=[{'project': 'HCV',
                                                      'region': 'HCV1-E1',
                                                      'seed': 'HCV1-seed',
                                                      'on.score': '2'},
                                                     {'project': 'HCV',
                                                      'region': 'HCV2-E1',
                                                      'seed': 'HCV2-seed',
                                                      'on.score': '3'}],
                                    remap_counts=[{'type': 'remap-1 HCV1-seed'},
                                                  {'type': 'remap-1 HCV2-seed'}]),
                        SampleFiles(coverage_scores=[{'project': 'HCV',
                                                      'region': 'HCV1-E1',
                                                      'seed': 'HCV1-seed',
                                                      'on.score': '2'}],
                                    remap_counts=[{'type': 'remap-1 HCV1-seed'}]))
        expected_report = ''
        expected_scenario_counts = {Scenarios.REMAP_COUNTS_CHANGED: [
            '  run1:sample42 coverage: HCV HCV2-E1 3 => -\n']}

        report, scenario_counts, _ = compare_sample(sample, Scenarios.REMAP_COUNTS_CHANGED)

        self.assertEqual(expected_report, report)
        self.assertEqual(expected_scenario_counts, scenario_counts)

    def test_missing_coverage_different_remap_counts_no_scenarios(self):
        sample = Sample(MiseqRun(target_path='run1/Results/versionX'),
                        'sample42',
                        SampleFiles(coverage_scores=[{'project': 'HCV',
                                                      'region': 'HCV1-E1',
                                                      'seed': 'HCV1-seed',
                                                      'on.score': '2'},
                                                     {'project': 'HCV',
                                                      'region': 'HCV2-E1',
                                                      'seed': 'HCV2-seed',
                                                      'on.score': '3'}],
                                    remap_counts=[{'type': 'remap-1 HCV1-seed'},
                                                  {'type': 'remap-1 HCV2-seed'}]),
                        SampleFiles(coverage_scores=[{'project': 'HCV',
                                                      'region': 'HCV1-E1',
                                                      'seed': 'HCV1-seed',
                                                      'on.score': '2'}],
                                    remap_counts=[{'type': 'remap-1 HCV1-seed'}]))
        expected_report = 'run1:sample42 coverage: HCV HCV2-E1 3 => -\n'
        expected_scenario_counts = {}

        report, scenario_counts, _ = compare_sample(sample, Scenarios.NONE)

        self.assertEqual(expected_report, report)
        self.assertEqual(expected_scenario_counts, scenario_counts)

    def test_missing_coverage_to_low(self):
        sample = Sample(MiseqRun(target_path='run1/Results/versionX'),
                        'sample42',
                        SampleFiles(coverage_scores=None),
                        SampleFiles(coverage_scores=[{'project': 'HIV',
                                                      'region': 'PR',
                                                      'on.score': '1'}]))
        expected_report = ''

        report, _, _ = compare_sample(sample)

        self.assertEqual(expected_report, report)

    def test_different_remap_counts(self):
        sample = Sample(MiseqRun(target_path='run1/Results/versionX'),
                        'sample42',
                        SampleFiles(coverage_scores=[{'project': 'HIV',
                                                      'region': 'PR',
                                                      'seed': 'HIV1-seed',
                                                      'on.score': '3'}],
                                    remap_counts=[{'type': 'remap-1 HIV1-seed'}]),
                        SampleFiles(coverage_scores=[{'project': 'HIV',
                                                      'region': 'PR',
                                                      'seed': 'HIV1-seed',
                                                      'on.score': '4'}],
                                    remap_counts=[{'type': 'remap-1 HIV1-seed'},
                                                  {'type': 'remap-2 HIV1-seed'}]))
        expected_report = ''
        expected_scenario_counts = {Scenarios.REMAP_COUNTS_CHANGED: [
            '  run1:sample42 coverage: HIV PR 3 => 4\n']}

        report, scenario_counts, _ = compare_sample(sample, Scenarios.REMAP_COUNTS_CHANGED)

        self.assertEqual(expected_report, report)
        self.assertEqual(expected_scenario_counts, scenario_counts)

    def test_remap_counts_ignore_raw_and_prelim(self):
        sample = Sample(MiseqRun(target_path='run1/Results/versionX'),
                        'sample42',
                        SampleFiles(coverage_scores=[{'project': 'HIV',
                                                      'region': 'PR',
                                                      'seed': 'HIV1-seed',
                                                      'on.score': '3'}],
                                    remap_counts=[{'type': 'raw'},
                                                  {'type': 'remap-1 HIV1-seed'}]),
                        SampleFiles(coverage_scores=[{'project': 'HIV',
                                                      'region': 'PR',
                                                      'seed': 'HIV1-seed',
                                                      'on.score': '4'}],
                                    remap_counts=[{'type': 'prelim HIV1-seed'},
                                                  {'type': 'remap-1 HIV1-seed'},
                                                  {'type': 'remap-2 HIV1-seed'}]))
        expected_report = ''
        expected_scenario_counts = {Scenarios.REMAP_COUNTS_CHANGED: [
            '  run1:sample42 coverage: HIV PR 3 => 4\n']}

        report, scenario_counts, _ = compare_sample(sample, Scenarios.REMAP_COUNTS_CHANGED)

        self.assertEqual(expected_report, report)
        self.assertEqual(expected_scenario_counts, scenario_counts)

    def test_remap_counts_ignore_final(self):
        sample = Sample(MiseqRun(target_path='run1/Results/versionX'),
                        'sample42',
                        SampleFiles(coverage_scores=[{'project': 'HIV',
                                                      'region': 'PR',
                                                      'seed': 'HIV1-seed',
                                                      'on.score': '3'}],
                                    remap_counts=[{'type': 'remap-1 HIV1-seed'},
                                                  {'type': 'remap-final HIV1-seed'}]),
                        SampleFiles(coverage_scores=[{'project': 'HIV',
                                                      'region': 'PR',
                                                      'seed': 'HIV1-seed',
                                                      'on.score': '4'}],
                                    remap_counts=[{'type': 'remap-1 HIV1-seed'},
                                                  {'type': 'remap-2 HIV1-seed'},
                                                  {'type': 'remap-final HIV1-seed'}]))
        expected_report = ''
        expected_scenario_counts = {Scenarios.REMAP_COUNTS_CHANGED: [
            '  run1:sample42 coverage: HIV PR 3 => 4\n']}

        report, scenario_counts, _ = compare_sample(sample, Scenarios.REMAP_COUNTS_CHANGED)

        self.assertEqual(expected_report, report)
        self.assertEqual(expected_scenario_counts, scenario_counts)

    def test_remap_counts_ignore_other_seeds(self):
        sample = Sample(MiseqRun(target_path='run1/Results/versionX'),
                        'sample42',
                        SampleFiles(coverage_scores=[{'project': 'HIV',
                                                      'region': 'PR',
                                                      'seed': 'HIV1-seed',
                                                      'on.score': '3'}],
                                    remap_counts=[{'type': 'remap-1 HIV1-seed'},
                                                  {'type': 'remap-1 HIV2-seed'},
                                                  {'type': 'remap-2 HIV2-seed'}]),
                        SampleFiles(coverage_scores=[{'project': 'HIV',
                                                      'region': 'PR',
                                                      'seed': 'HIV1-seed',
                                                      'on.score': '4'}],
                                    remap_counts=[{'type': 'remap-1 HIV1-seed'},
                                                  {'type': 'remap-2 HIV1-seed'}]))
        expected_report = ''
        expected_scenario_counts = {Scenarios.REMAP_COUNTS_CHANGED: [
            '  run1:sample42 coverage: HIV PR 3 => 4\n']}

        report, scenario_counts, _ = compare_sample(sample, Scenarios.REMAP_COUNTS_CHANGED)

        self.assertEqual(expected_report, report)
        self.assertEqual(expected_scenario_counts, scenario_counts)

    def test_equivalent_remap_counts(self):
        sample = Sample(MiseqRun(target_path='run1/Results/versionX'),
                        'sample42',
                        SampleFiles(coverage_scores=[{'project': 'HIV',
                                                      'region': 'PR',
                                                      'seed': 'HIV1-seed',
                                                      'on.score': '3'}],
                                    remap_counts=[{'type': 'remap-1 HIV1-seed',
                                                   'count': '20'}]),
                        SampleFiles(coverage_scores=[{'project': 'HIV',
                                                      'region': 'PR',
                                                      'seed': 'HIV1-seed',
                                                      'on.score': '4'}],
                                    remap_counts=[{'type': 'remap-1 HIV1-seed',
                                                   'count': '21'}]))
        expected_report = 'run1:sample42 coverage: HIV PR 3 => 4\n'
        expected_scenario_counts = {}

        report, scenario_counts, _ = compare_sample(sample)

        self.assertEqual(expected_report, report)
        self.assertEqual(expected_scenario_counts, scenario_counts)

    def test_consensus_change(self):
        sample = Sample(MiseqRun(target_path='run1/Results/versionX'),
                        'sample42',
                        SampleFiles(coverage_scores=[{'region': 'R1',
                                                      'seed': 'R1-seed',
                                                      'on.score': '4'}],
                                    nuc_limits={'R1-seed': [('R1', 101, 108)]},
                                    consensus=[{'region': 'R1-seed',
                                                'consensus-percent-cutoff': 'MAX',
                                                'offset': '100',
                                                'sequence': 'ACACACGT'}]),
                        SampleFiles(coverage_scores=[{'region': 'R1',
                                                      'seed': 'R1-seed',
                                                      'on.score': '4'}],
                                    nuc_limits={'R1-seed': [('R1', 101, 108)]},
                                    consensus=[{'region': 'R1-seed',
                                                'consensus-percent-cutoff': 'MAX',
                                                'offset': '100',
                                                'sequence': 'ACACACGG'}]))
        expected_report = ('run1:sample42 consensus: R1-seed R1 MAX\n'
                           '- ACACACGT\n'
                           '?        ^\n'
                           '+ ACACACGG\n'
                           '?        ^\n')
        expected_consensus_distances = [ConsensusDistance(region='R1',
                                                          cutoff='MAX',
                                                          distance=1,
                                                          pct_diff=12.5)]

        report, _, consensus_distances = compare_sample(sample)

        self.assertEqual(expected_report, report)
        self.assertEqual(expected_consensus_distances, consensus_distances)

    def test_consensus_change_diff(self):
        source_seqs = {('R1-seed', 'R1', 'MAX'): 'ACACAC'}
        target_seqs = {('R1-seed', 'R1', 'MAX'): 'ACACAT'}
        sample = Sample(MiseqRun(target_path='run1/Results/versionX'),
                        'sample42',
                        SampleFiles(),
                        SampleFiles())
        expected_diffs = ['run1:sample42 consensus: R1-seed R1 MAX',
                          '- ACACAC',
                          '?      ^',
                          '+ ACACAT',
                          '?      ^']
        expected_scenarios = {}
        diffs = []
        scenarios = defaultdict(list)

        compare_consensus(
            sample,
            source_seqs,
            target_seqs,
            diffs,
            Scenarios.NONE,
            scenarios)

        self.assertEqual(expected_diffs, diffs)
        self.assertEqual(expected_scenarios, scenarios)

    def test_consensus_change_scenario(self):
        source_seqs = {('R1-seed', 'R1', 'MAX'): 'ACACAC'}
        target_seqs = {('R1-seed', 'R1', 'MAX'): 'ACACAT'}
        sample = Sample(MiseqRun(target_path='run1/Results/versionX'),
                        'sample42',
                        SampleFiles(),
                        SampleFiles())
        expected_diffs = []
        expected_scenarios = {Scenarios.MAIN_CONSENSUS_CHANGED: ['.']}
        diffs = []
        scenarios = defaultdict(list)

        compare_consensus(
            sample,
            source_seqs,
            target_seqs,
            diffs,
            Scenarios.MAIN_CONSENSUS_CHANGED | Scenarios.OTHER_CONSENSUS_CHANGED,
            scenarios)

        self.assertEqual(expected_diffs, diffs)
        self.assertEqual(expected_scenarios, scenarios)

    def test_other_consensus_change(self):
        source_seqs = {('R1-seed', 'R1', '0.250'): 'ACACAC'}
        target_seqs = {('R1-seed', 'R1', '0.250'): 'ACACAT'}
        sample = Sample(MiseqRun(target_path='run1/Results/versionX'),
                        'sample42',
                        SampleFiles(),
                        SampleFiles())
        expected_diffs = ['run1:sample42 consensus: R1-seed R1 0.250',
                          '- ACACAC',
                          '?      ^',
                          '+ ACACAT',
                          '?      ^']
        expected_scenarios = {}
        diffs = []
        scenarios = defaultdict(list)

        compare_consensus(
            sample,
            source_seqs,
            target_seqs,
            diffs,
            Scenarios.NONE,
            scenarios)

        self.assertEqual(expected_diffs, diffs)
        self.assertEqual(expected_scenarios, scenarios)

    def test_other_consensus_change_scenario(self):
        source_seqs = {('R1-seed', 'R1', '0.250'): 'ACACAC'}
        target_seqs = {('R1-seed', 'R1', '0.250'): 'ACACAT'}
        sample = Sample(MiseqRun(target_path='run1/Results/versionX'),
                        'sample42',
                        SampleFiles(),
                        SampleFiles())
        expected_diffs = []
        expected_scenarios = {Scenarios.OTHER_CONSENSUS_CHANGED: ['.']}
        diffs = []
        scenarios = defaultdict(list)

        compare_consensus(
            sample,
            source_seqs,
            target_seqs,
            diffs,
            Scenarios.MAIN_CONSENSUS_CHANGED | Scenarios.OTHER_CONSENSUS_CHANGED,
            scenarios)

        self.assertEqual(expected_diffs, diffs)
        self.assertEqual(expected_scenarios, scenarios)

    def test_hla_consensus_change_scenario(self):
        source_seqs = {('HLA1-seed', 'HLA-1', '0.250'): 'ACACAC'}
        target_seqs = {('HLA1-seed', 'HLA-1', '0.250'): 'ACACAT'}
        sample = Sample(MiseqRun(target_path='run1/Results/versionX'),
                        'sample42',
                        SampleFiles(),
                        SampleFiles())
        expected_diffs = []
        expected_scenarios = {Scenarios.MAIN_CONSENSUS_CHANGED: ['.']}
        diffs = []
        scenarios = defaultdict(list)

        compare_consensus(
            sample,
            source_seqs,
            target_seqs,
            diffs,
            Scenarios.MAIN_CONSENSUS_CHANGED | Scenarios.OTHER_CONSENSUS_CHANGED,
            scenarios)

        self.assertEqual(expected_diffs, diffs)
        self.assertEqual(expected_scenarios, scenarios)

    def test_same_consensus(self):
        source_seqs = {('R1-seed', 'R1', 'MAX'): 'ACACAC'}
        target_seqs = {('R1-seed', 'R1', 'MAX'): 'ACACAC'}
        sample = Sample(MiseqRun(target_path='run1/Results/versionX'),
                        'sample42',
                        SampleFiles(),
                        SampleFiles())
        expected_diffs = []
        expected_scenarios = {}
        diffs = []
        scenarios = defaultdict(list)

        compare_consensus(
            sample,
            source_seqs,
            target_seqs,
            diffs,
            Scenarios.MAIN_CONSENSUS_CHANGED | Scenarios.OTHER_CONSENSUS_CHANGED,
            scenarios)

        self.assertEqual(expected_diffs, diffs)
        self.assertEqual(expected_scenarios, scenarios)

    def test_map_consensus_sequences(self):
        files = SampleFiles(coverage_scores=[{'region': 'R1',
                                              'seed': 'R1-seed',
                                              'on.score': '4'}],
                            nuc_limits={'R1-seed': [('R1', 101, 108)]},
                            consensus=[{'region': 'R1-seed',
                                        'consensus-percent-cutoff': 'MAX',
                                        'offset': '100',
                                        'sequence': 'ACACACGG'}])
        expected_sequences = {('R1-seed', 'R1', 'MAX'): 'ACACACGG'}

        sequences = map_consensus_sequences(files)

        self.assertEqual(expected_sequences, sequences)

    def test_consensus_ignores_trailing_dashes(self):
        files = SampleFiles(coverage_scores=[{'region': 'R1',
                                              'seed': 'R1-seed',
                                              'on.score': '4'}],
                            nuc_limits={'R1-seed': [('R1', 101, 110)]},
                            consensus=[{'region': 'R1-seed',
                                        'consensus-percent-cutoff': 'MAX',
                                        'offset': '100',
                                        'sequence': 'ACACACGG--'}])
        expected_sequences = {('R1-seed', 'R1', 'MAX'): 'ACACACGG'}

        sequences = map_consensus_sequences(files)

        self.assertEqual(expected_sequences, sequences)

    def test_one_consensus_changes(self):
        source_seqs = {('R1-seed', 'R1', 'MAX'): 'ACACACGT',
                       ('R2-seed', 'R2', 'MAX'): 'ACACACGT'}
        target_seqs = {('R1-seed', 'R1', 'MAX'): 'ACACACGT',
                       ('R2-seed', 'R2', 'MAX'): 'ACACAMGT'}
        sample = Sample(MiseqRun(target_path='run1/Results/versionX'),
                        'sample42',
                        SampleFiles(),
                        SampleFiles())
        expected_diffs = ['run1:sample42 consensus: R2-seed R2 MAX',
                          '- ACACACGT',
                          '?      ^',
                          '+ ACACAMGT',
                          '?      ^']
        expected_consensus_distances = [ConsensusDistance(region='R1',
                                                          cutoff='MAX',
                                                          distance=0,
                                                          pct_diff=0),
                                        ConsensusDistance(region='R2',
                                                          cutoff='MAX',
                                                          distance=1,
                                                          pct_diff=12.5)]
        diffs = []
        scenarios = defaultdict(list)

        consensus_distances = compare_consensus(
            sample,
            source_seqs,
            target_seqs,
            diffs,
            Scenarios.NONE,
            scenarios)

        self.assertEqual(expected_diffs, diffs)
        self.assertEqual(expected_consensus_distances, consensus_distances)

    def test_consensus_trailing_change(self):
        source_seqs = {('R1-seed', 'R1', 'MAX'): 'ACTTAC------GTAC'}
        target_seqs = {('R1-seed', 'R1', 'MAX'): 'ACTTAC'}
        sample = Sample(MiseqRun(target_path='run1/Results/versionX'),
                        'sample42',
                        SampleFiles(),
                        SampleFiles())
        expected_diffs = ['run1:sample42 consensus: R1-seed R1 MAX',
                          '- ACTTAC------GTAC',
                          '+ ACTTAC']
        expected_consensus_distances = [ConsensusDistance(region='R1',
                                                          cutoff='MAX',
                                                          distance=4,
                                                          pct_diff=25)]
        diffs = []
        scenarios = defaultdict(list)

        consensus_distances = compare_consensus(
            sample,
            source_seqs,
            target_seqs,
            diffs,
            Scenarios.NONE,
            scenarios)

        self.assertEqual(expected_diffs, diffs)
        self.assertEqual(expected_consensus_distances, consensus_distances)

    def test_consensus_missing(self):
        source_seqs = {('R1-seed', 'R1', 'MAX'): 'ACTTAC'}
        target_seqs = {}
        sample = Sample(MiseqRun(target_path='run1/Results/versionX'),
                        'sample42',
                        SampleFiles(),
                        SampleFiles())
        expected_diffs = ['run1:sample42 consensus: R1-seed R1 MAX',
                          '- ACTTAC']
        expected_consensus_distances = []
        diffs = []
        scenarios = defaultdict(list)

        consensus_distances = compare_consensus(
            sample,
            source_seqs,
            target_seqs,
            diffs,
            Scenarios.NONE,
            scenarios)

        self.assertEqual(expected_diffs, diffs)
        self.assertEqual(expected_consensus_distances, consensus_distances)

    def test_consensus_added(self):
        source_seqs = {}
        target_seqs = {('R1-seed', 'R1', 'MAX'): 'ACTTAC'}
        sample = Sample(MiseqRun(target_path='run1/Results/versionX'),
                        'sample42',
                        SampleFiles(),
                        SampleFiles())
        expected_diffs = ['run1:sample42 consensus: R1-seed R1 MAX',
                          '+ ACTTAC']
        expected_consensus_distances = []
        diffs = []
        scenarios = defaultdict(list)

        consensus_distances = compare_consensus(
            sample,
            source_seqs,
            target_seqs,
            diffs,
            Scenarios.NONE,
            scenarios)

        self.assertEqual(expected_diffs, diffs)
        self.assertEqual(expected_consensus_distances, consensus_distances)

    def test_hiv_seed_changed(self):
        sample = Sample(MiseqRun(target_path='run1/Results/versionX'),
                        'sample42',
                        SampleFiles(consensus=[{'region': 'HIV1-X',
                                                'consensus-percent-cutoff': 'MAX',
                                                'offset': '100',
                                                'sequence': 'ACACAC'}]),
                        SampleFiles(consensus=[{'region': 'HIV1-Y',
                                                'consensus-percent-cutoff': 'MAX',
                                                'offset': '200',
                                                'sequence': 'ACACAC'}]))
        expected_report = ''

        report, _, _ = compare_sample(sample)

        self.assertEqual(expected_report, report)

    def test_consensus_needs_coverage(self):
        files = SampleFiles(coverage_scores=[{'region': 'R1',
                                              'seed': 'R1-seed',
                                              'on.score': '1'},
                                             {'region': 'R2',
                                              'seed': 'R2-seed',
                                              'on.score': '4'}],
                            nuc_limits={'R1-seed': [('R1', 101, 108)],
                                        'R2-seed': [('R2', 101, 108)]},
                            consensus=[{'region': 'R1-seed',
                                        'consensus-percent-cutoff': 'MAX',
                                        'offset': '100',
                                        'sequence': 'ACACACGG'},
                                       {'region': 'R2-seed',
                                        'consensus-percent-cutoff': 'MAX',
                                        'offset': '100',
                                        'sequence': 'ACACACGG'}])
        expected_sequences = {('R2-seed', 'R2', 'MAX'): 'ACACACGG'}

        sequences = map_consensus_sequences(files)

        self.assertEqual(expected_sequences, sequences)

    def test_map_consensus_multiple_regions(self):
        files = SampleFiles(coverage_scores=[{'region': 'R1a',
                                              'seed': 'R1-seed',
                                              'on.score': '4'},
                                             {'region': 'R1b',
                                              'seed': 'R1-seed',
                                              'on.score': '4'}],
                            nuc_limits={'R1-seed': [('R1a', 99, 103),
                                                    ('R1b', 109, 113)]},
                            consensus=[{'region': 'R1-seed',
                                        'consensus-percent-cutoff': 'MAX',
                                        'offset': '98',
                                        'sequence': 'ACACATGTGTAGAGATT'}])
        expected_sequences = {('R1-seed', 'R1a', 'MAX'): 'ACACA',
                              ('R1-seed', 'R1b', 'MAX'): 'AGAGA'}

        sequences = map_consensus_sequences(files)

        self.assertEqual(expected_sequences, sequences)

    def test_map_consensus_offset_past_start(self):
        files = SampleFiles(coverage_scores=[{'region': 'R1',
                                              'seed': 'R1-seed',
                                              'on.score': '4'}],
                            nuc_limits={'R1-seed': [('R1', 96, 105)]},
                            consensus=[{'region': 'R1-seed',
                                        'consensus-percent-cutoff': 'MAX',
                                        'offset': '100',
                                        'sequence': 'ACACATGTGT'}])
        expected_sequences = {('R1-seed', 'R1', 'MAX'): '-----ACACA'}

        sequences = map_consensus_sequences(files)

        self.assertEqual(expected_sequences, sequences)

    def test_map_consensus_no_consensus(self):
        files = SampleFiles(coverage_scores=[{'region': 'R1',
                                              'seed': 'R1-seed',
                                              'on.score': '4'}],
                            nuc_limits={'R1-seed': [('R1', 96, 105)]},
                            consensus=None)
        expected_sequences = {}

        sequences = map_consensus_sequences(files)

        self.assertEqual(expected_sequences, sequences)

    def test_consensus_groups_hcv(self):
        files = SampleFiles(coverage_scores=[{'region': 'HCV1A-H77-NS4b',
                                              'seed': 'HCV1a-seed',
                                              'on.score': '4'}],
                            nuc_limits={'HCV1a-seed': [('HCV1A-H77-NS4b', 101, 108)]},
                            consensus=[{'region': 'HCV1a-seed',
                                        'consensus-percent-cutoff': 'MAX',
                                        'offset': '100',
                                        'sequence': 'ACACACGG'}])
        expected_sequences = {('HCV1a-seed', 'HCV-NS4b', 'MAX'): 'ACACACGG'}

        sequences = map_consensus_sequences(files)

        self.assertEqual(expected_sequences, sequences)


class GroupSamplesTest(TestCase):
    def test_group_samples_file(self):
        output_file = StringIO("""\
sample,a,b
s1,a11,b11
s1,a12,b12
s2,a2,b2
""")
        expected_groups = dict(s1=[dict(sample='s1', a='a11', b='b11'),
                                   dict(sample='s1', a='a12', b='b12')],
                               s2=[dict(sample='s2', a='a2', b='b2')])

        groups = group_samples_file(output_file)

        self.assertEqual(expected_groups, groups)

    def test_group_nucs_single(self):
        output_file = StringIO("""\
sample,seed,region,query.nuc.pos
s1,R1-seed,R1,100
s1,R1-seed,R1,110
""")
        expected_groups = dict(s1={'R1-seed': [('R1', 100, 110)]})

        groups = group_nucs_file(output_file)

        self.assertEqual(expected_groups, groups)

    def test_group_nucs_two_samples(self):
        output_file = StringIO("""\
sample,seed,region,query.nuc.pos
s1,R1-seed,R1,100
s1,R1-seed,R1,110
s2,R1-seed,R1,101
s2,R1-seed,R1,111
""")
        expected_groups = dict(s1={'R1-seed': [('R1', 100, 110)]},
                               s2={'R1-seed': [('R1', 101, 111)]})

        groups = group_nucs_file(output_file)

        self.assertEqual(expected_groups, groups)

    def test_group_nucs_two_regions(self):
        output_file = StringIO("""\
sample,seed,region,query.nuc.pos
s1,R1-seed,R1a,100
s1,R1-seed,R1a,110
s1,R1-seed,R1b,201
s1,R1-seed,R1b,211
""")
        expected_groups = dict(s1={'R1-seed': [('R1a', 100, 110),
                                               ('R1b', 201, 211)]})

        groups = group_nucs_file(output_file)

        self.assertEqual(expected_groups, groups)

    def test_group_nucs_blank(self):
        output_file = StringIO("""\
sample,seed,region,query.nuc.pos
s1,R1-seed,R1,100
s1,R1-seed,R1,110
s1,R1-seed,R1,
""")
        expected_groups = dict(s1={'R1-seed': [('R1', 100, 110)]})

        groups = group_nucs_file(output_file)

        self.assertEqual(expected_groups, groups)
