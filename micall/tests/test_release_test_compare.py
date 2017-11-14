from unittest import TestCase

from release_test_compare import compare_sample, SampleFiles, Sample, MiseqRun, Scenarios


class CompareSampleTest(TestCase):
    def test_empty(self):
        sample = Sample(MiseqRun(target_path='run1/Results/versionX'),
                        'sample42',
                        SampleFiles(),
                        SampleFiles())
        expected_report = ''

        report, _ = compare_sample(sample)

        self.assertEqual(expected_report, report)

    def test_x4_big_change(self):
        sample = Sample(MiseqRun(target_path='run1/Results/versionX'),
                        'sample42',
                        SampleFiles(g2p_summary=[dict(X4pct='50.00')]),
                        SampleFiles(g2p_summary=[dict(X4pct='60.00')]))
        expected_report = 'run1:sample42 G2P: 50.00 => 60.00\n'

        report, _ = compare_sample(sample)

        self.assertEqual(expected_report, report)

    def test_x4_small_change(self):
        sample = Sample(MiseqRun(target_path='run1/Results/versionX'),
                        'sample42',
                        SampleFiles(g2p_summary=[dict(X4pct='50.00')]),
                        SampleFiles(g2p_summary=[dict(X4pct='51.00')]))
        expected_report = ''

        report, _ = compare_sample(sample)

        self.assertEqual(expected_report, report)

    def test_blank(self):
        sample = Sample(MiseqRun(target_path='run1/Results/versionX'),
                        'sample42',
                        SampleFiles(g2p_summary=[dict(X4pct='50.00')]),
                        SampleFiles(g2p_summary=[dict(X4pct='')]))
        expected_report = 'run1:sample42 G2P: 50.00 => \n'

        report, _ = compare_sample(sample)

        self.assertEqual(expected_report, report)

    def test_other_difference_with_blanks(self):
        sample = Sample(MiseqRun(target_path='run1/Results/versionX'),
                        'sample42',
                        SampleFiles(g2p_summary=[dict(X4pct='', other='x')]),
                        SampleFiles(g2p_summary=[dict(X4pct='', other='y')]))
        expected_report = ''

        report, _ = compare_sample(sample)

        self.assertEqual(expected_report, report)

    def test_same_final(self):
        sample = Sample(MiseqRun(target_path='run1/Results/versionX'),
                        'sample42',
                        SampleFiles(g2p_summary=[dict(X4pct='2.99', final='X4')]),
                        SampleFiles(g2p_summary=[dict(X4pct='3.01', final='X4')]))
        expected_report = ''

        report, _ = compare_sample(sample)

        self.assertEqual(expected_report, report)

    def test_different_final(self):
        sample = Sample(MiseqRun(target_path='run1/Results/versionX'),
                        'sample42',
                        SampleFiles(g2p_summary=[dict(X4pct='1.99', final='R5')]),
                        SampleFiles(g2p_summary=[dict(X4pct='2.01', final='X4')]))
        expected_report = 'run1:sample42 G2P: R5 1.99 => X4 2.01\n'

        report, _ = compare_sample(sample)

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

        report, _ = compare_sample(sample)

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

        report, _ = compare_sample(sample)

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

        report, _ = compare_sample(sample)

        self.assertEqual(expected_report, report)

    def test_missing_coverage(self):
        sample = Sample(MiseqRun(target_path='run1/Results/versionX'),
                        'sample42',
                        SampleFiles(coverage_scores=None),
                        SampleFiles(coverage_scores=[{'project': 'HIV',
                                                      'region': 'PR',
                                                      'on.score': '2'}]))
        expected_report = 'run1:sample42 coverage: HIV PR - => 2\n'

        report, _ = compare_sample(sample)

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

        report, scenario_counts = compare_sample(sample)

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

        report, scenario_counts = compare_sample(sample, Scenarios.REMAP_COUNTS_CHANGED)

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

        report, scenario_counts = compare_sample(sample, Scenarios.NONE)

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

        report, _ = compare_sample(sample)

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

        report, scenario_counts = compare_sample(sample, Scenarios.REMAP_COUNTS_CHANGED)

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

        report, scenario_counts = compare_sample(sample, Scenarios.REMAP_COUNTS_CHANGED)

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

        report, scenario_counts = compare_sample(sample, Scenarios.REMAP_COUNTS_CHANGED)

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

        report, scenario_counts = compare_sample(sample, Scenarios.REMAP_COUNTS_CHANGED)

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

        report, scenario_counts = compare_sample(sample)

        self.assertEqual(expected_report, report)
        self.assertEqual(expected_scenario_counts, scenario_counts)

    def test_removed_key_position(self):
        sample = Sample(MiseqRun(target_path='run1/Results/versionX'),
                        'sample42',
                        SampleFiles(coverage_scores=[{'project': 'HIV',
                                                      'region': 'RT',
                                                      'which.key.pos': '318',
                                                      'on.score': '1'}]),
                        SampleFiles(coverage_scores=[{'project': 'HIV',
                                                      'region': 'RT',
                                                      'which.key.pos': '230',
                                                      'on.score': '4'}]))
        expected_report = ''
        expected_scenario_counts = {Scenarios.V78_KEY_POS_REMOVED_RT318: [
            '  run1:sample42 coverage: HIV RT 1 => 4\n']}

        report, scenario_counts = compare_sample(
            sample,
            Scenarios.V78_KEY_POS_REMOVED_RT318)

        self.assertEqual(expected_report, report)
        self.assertEqual(expected_scenario_counts, scenario_counts)

    def test_consensus_change(self):
        sample = Sample(MiseqRun(target_path='run1/Results/versionX'),
                        'sample42',
                        SampleFiles(consensus=[{'region': 'R1',
                                                'consensus-percent-cutoff': 'MAX',
                                                'offset': '100',
                                                'sequence': 'ACACAC'}]),
                        SampleFiles(consensus=[{'region': 'R1',
                                                'consensus-percent-cutoff': 'MAX',
                                                'offset': '100',
                                                'sequence': 'ACACAT'}]))
        expected_report = ('run1:sample42 consensus: R1 MAX\n'
                           '- 100 ACACAC\n'
                           '?          ^\n'
                           '+ 100 ACACAT\n'
                           '?          ^\n')

        report, _ = compare_sample(sample)

        self.assertEqual(expected_report, report)

    def test_consensus_change_scenario(self):
        sample = Sample(MiseqRun(target_path='run1/Results/versionX'),
                        'sample42',
                        SampleFiles(consensus=[{'region': 'R1',
                                                'consensus-percent-cutoff': 'MAX',
                                                'offset': '100',
                                                'sequence': 'ACACAC'}]),
                        SampleFiles(consensus=[{'region': 'R1',
                                                'consensus-percent-cutoff': 'MAX',
                                                'offset': '100',
                                                'sequence': 'ACACAT'}]))
        expected_report = ''
        expected_scenarios = {Scenarios.MAIN_CONSENSUS_CHANGED: ['.']}

        report, scenarios = compare_sample(
            sample,
            Scenarios.MAIN_CONSENSUS_CHANGED | Scenarios.OTHER_CONSENSUS_CHANGED)

        self.assertEqual(expected_report, report)
        self.assertEqual(expected_scenarios, scenarios)

    def test_other_consensus_change(self):
        sample = Sample(MiseqRun(target_path='run1/Results/versionX'),
                        'sample42',
                        SampleFiles(consensus=[{'region': 'R1',
                                                'consensus-percent-cutoff': '0.250',
                                                'offset': '100',
                                                'sequence': 'ACACAC'}]),
                        SampleFiles(consensus=[{'region': 'R1',
                                                'consensus-percent-cutoff': '0.250',
                                                'offset': '100',
                                                'sequence': 'ACACAT'}]))
        expected_report = ('run1:sample42 consensus: R1 0.250\n'
                           '- 100 ACACAC\n'
                           '?          ^\n'
                           '+ 100 ACACAT\n'
                           '?          ^\n')

        report, _ = compare_sample(sample)

        self.assertEqual(expected_report, report)

    def test_other_consensus_change_scenario(self):
        sample = Sample(MiseqRun(target_path='run1/Results/versionX'),
                        'sample42',
                        SampleFiles(consensus=[{'region': 'R1',
                                                'consensus-percent-cutoff': '0.250',
                                                'offset': '100',
                                                'sequence': 'ACACAC'}]),
                        SampleFiles(consensus=[{'region': 'R1',
                                                'consensus-percent-cutoff': '0.250',
                                                'offset': '100',
                                                'sequence': 'ACACAT'}]))
        expected_report = ''
        expected_scenarios = {Scenarios.OTHER_CONSENSUS_CHANGED: ['.']}

        report, scenarios = compare_sample(
            sample,
            Scenarios.MAIN_CONSENSUS_CHANGED | Scenarios.OTHER_CONSENSUS_CHANGED)

        self.assertEqual(expected_report, report)
        self.assertEqual(expected_scenarios, scenarios)

    def test_hla_consensus_change_scenario(self):
        sample = Sample(MiseqRun(target_path='run1/Results/versionX'),
                        'sample42',
                        SampleFiles(consensus=[{'region': 'HLA-1',
                                                'consensus-percent-cutoff': '0.250',
                                                'offset': '100',
                                                'sequence': 'ACACAC'}]),
                        SampleFiles(consensus=[{'region': 'HLA-1',
                                                'consensus-percent-cutoff': '0.250',
                                                'offset': '100',
                                                'sequence': 'ACACAT'}]))
        expected_report = ''
        expected_scenarios = {Scenarios.MAIN_CONSENSUS_CHANGED: ['.']}

        report, scenarios = compare_sample(
            sample,
            Scenarios.MAIN_CONSENSUS_CHANGED | Scenarios.OTHER_CONSENSUS_CHANGED)

        self.assertEqual(expected_report, report)
        self.assertEqual(expected_scenarios, scenarios)

    def test_same_consensus(self):
        sample = Sample(MiseqRun(target_path='run1/Results/versionX'),
                        'sample42',
                        SampleFiles(consensus=[{'region': 'R1',
                                                'consensus-percent-cutoff': 'MAX',
                                                'offset': '100',
                                                'sequence': 'ACACAC'}]),
                        SampleFiles(consensus=[{'region': 'R1',
                                                'consensus-percent-cutoff': 'MAX',
                                                'offset': '100',
                                                'sequence': 'ACACAC'}]))
        expected_report = ''

        report, _ = compare_sample(sample)

        self.assertEqual(expected_report, report)

    def test_consensus_ignores_trailing_dashes(self):
        sample = Sample(MiseqRun(target_path='run1/Results/versionX'),
                        'sample42',
                        SampleFiles(consensus=[{'region': 'R1',
                                                'consensus-percent-cutoff': 'MAX',
                                                'offset': '100',
                                                'sequence': 'ACACAC'}]),
                        SampleFiles(consensus=[{'region': 'R1',
                                                'consensus-percent-cutoff': 'MAX',
                                                'offset': '100',
                                                'sequence': 'ACACAC---'}]))
        expected_report = ''

        report, _ = compare_sample(sample)

        self.assertEqual(expected_report, report)

    def test_one_consensus_changes(self):
        sample = Sample(MiseqRun(target_path='run1/Results/versionX'),
                        'sample42',
                        SampleFiles(consensus=[{'region': 'R1',
                                                'consensus-percent-cutoff': 'MAX',
                                                'offset': '100',
                                                'sequence': 'ACACAC'},
                                               {'region': 'R2',
                                                'consensus-percent-cutoff': 'MAX',
                                                'offset': '100',
                                                'sequence': 'ACACAC'}]),
                        SampleFiles(consensus=[{'region': 'R1',
                                                'consensus-percent-cutoff': 'MAX',
                                                'offset': '100',
                                                'sequence': 'ACACAC'},
                                               {'region': 'R2',
                                                'consensus-percent-cutoff': 'MAX',
                                                'offset': '100',
                                                'sequence': 'ACACAM'}]))
        expected_report = ('run1:sample42 consensus: R2 MAX\n'
                           '- 100 ACACAC\n'
                           '?          ^\n'
                           '+ 100 ACACAM\n'
                           '?          ^\n')

        report, _ = compare_sample(sample)

        self.assertEqual(expected_report, report)

    def test_consensus_missing(self):
        sample = Sample(MiseqRun(target_path='run1/Results/versionX'),
                        'sample42',
                        SampleFiles(consensus=[{'region': 'R1',
                                                'consensus-percent-cutoff': 'MAX',
                                                'offset': '100',
                                                'sequence': 'ACACAC'}]),
                        SampleFiles(consensus=[]))
        expected_report = ('run1:sample42 consensus: R1 MAX\n'
                           '- 100 ACACAC\n')

        report, _ = compare_sample(sample)

        self.assertEqual(expected_report, report)

    def test_consensus_added(self):
        sample = Sample(MiseqRun(target_path='run1/Results/versionX'),
                        'sample42',
                        SampleFiles(consensus=[]),
                        SampleFiles(consensus=[{'region': 'R1',
                                                'consensus-percent-cutoff': 'MAX',
                                                'offset': '100',
                                                'sequence': 'ACACAC'}]))
        expected_report = ('run1:sample42 consensus: R1 MAX\n'
                           '+ 100 ACACAC\n')

        report, _ = compare_sample(sample)

        self.assertEqual(expected_report, report)

    def test_duplicate_consensus(self):
        sample = Sample(MiseqRun(target_path='run1/Results/versionX'),
                        'sample42',
                        SampleFiles(consensus=[{'region': 'R1',
                                                'consensus-percent-cutoff': 'MAX',
                                                'offset': '100',
                                                'sequence': 'ACACAC'}]),
                        SampleFiles(consensus=[{'region': 'R1',
                                                'consensus-percent-cutoff': 'MAX',
                                                'offset': '100',
                                                'sequence': 'ACACAC'},
                                               {'region': 'R1',
                                                'consensus-percent-cutoff': 'MAX',
                                                'offset': '100',
                                                'sequence': 'ACACAC'}]))
        expected_report = 'run1:sample42 duplicate consensus: R1 MAX\n'

        report, _ = compare_sample(sample)

        self.assertEqual(expected_report, report)

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

        report, _ = compare_sample(sample)

        self.assertEqual(expected_report, report)
