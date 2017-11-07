from unittest import TestCase

from release_test_compare import compare_sample, SampleFiles, Sample, MiseqRun


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

    def test_missing_coverage_seed(self):
        sample = Sample(MiseqRun(target_path='run1/Results/versionX'),
                        'sample42',
                        SampleFiles(coverage_scores=[{'project': 'HIV',
                                                      'region': 'PR',
                                                      'seed': 'HIV1-seed',
                                                      'on.score': '2'},
                                                     {'project': 'HCV',
                                                      'region': 'E2',
                                                      'seed': 'HCV1-seed',
                                                      'on.score': '2'}],
                                    remap_counts=[{'type': 'remap-1 HIV1-seed'},
                                                  {'type': 'remap-1 HCV1-seed'}]),
                        SampleFiles(coverage_scores=[{'project': 'HIV',
                                                      'region': 'PR',
                                                      'seed': 'HIV1-seed',
                                                      'on.score': '2'}],
                                    remap_counts=[{'type': 'remap-1 HIV1-seed'}]))
        expected_report = 'run1:sample42 coverage: HCV E2 2 => -\n'

        report, _ = compare_sample(sample)

        self.assertEqual(expected_report, report)

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
        expected_scenario_counts = {'different remap counts': 1}

        report, scenario_counts = compare_sample(sample)

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
        expected_scenario_counts = {'different remap counts': 1}

        report, scenario_counts = compare_sample(sample)

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
        expected_scenario_counts = {'different remap counts': 1}

        report, scenario_counts = compare_sample(sample)

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
        expected_scenario_counts = {'different remap counts': 1}

        report, scenario_counts = compare_sample(sample)

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
        expected_scenario_counts = {'removed key pos 318': 1}

        report, scenario_counts = compare_sample(sample)

        self.assertEqual(expected_report, report)
        self.assertEqual(expected_scenario_counts, scenario_counts)
