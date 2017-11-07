from unittest import TestCase

from release_test_compare import compare_sample, SampleFiles, Sample, MiseqRun


class CompareSampleTest(TestCase):
    def test_empty(self):
        sample = Sample(MiseqRun(target_path='run1/Results/versionX'),
                        'sample42',
                        SampleFiles(),
                        SampleFiles())
        expected_report = ''

        report = compare_sample(sample)

        self.assertEqual(expected_report, report)

    def test_x4_big_change(self):
        sample = Sample(MiseqRun(target_path='run1/Results/versionX'),
                        'sample42',
                        SampleFiles(g2p_summary=[dict(X4pct='50.00')]),
                        SampleFiles(g2p_summary=[dict(X4pct='60.00')]))
        expected_report = 'run1:sample42 G2P: 50.00 => 60.00\n'

        report = compare_sample(sample)

        self.assertEqual(expected_report, report)

    def test_x4_small_change(self):
        sample = Sample(MiseqRun(target_path='run1/Results/versionX'),
                        'sample42',
                        SampleFiles(g2p_summary=[dict(X4pct='50.00')]),
                        SampleFiles(g2p_summary=[dict(X4pct='51.00')]))
        expected_report = ''

        report = compare_sample(sample)

        self.assertEqual(expected_report, report)

    def test_blank(self):
        sample = Sample(MiseqRun(target_path='run1/Results/versionX'),
                        'sample42',
                        SampleFiles(g2p_summary=[dict(X4pct='50.00')]),
                        SampleFiles(g2p_summary=[dict(X4pct='')]))
        expected_report = 'run1:sample42 G2P: 50.00 => \n'

        report = compare_sample(sample)

        self.assertEqual(expected_report, report)

    def test(self):
        sample = Sample(MiseqRun(target_path='run1/Results/versionX'),
                        'sample42',
                        SampleFiles(g2p_summary=[dict(X4pct='', other='x')]),
                        SampleFiles(g2p_summary=[dict(X4pct='', other='y')]))
        expected_report = ''

        report = compare_sample(sample)

        self.assertEqual(expected_report, report)
