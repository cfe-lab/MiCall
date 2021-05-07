from collections import defaultdict
from io import StringIO
from unittest import TestCase

import typing

from release_test_compare import compare_sample, SampleFiles, Sample, \
    MiseqRun, Scenarios, ConsensusDistance, group_samples_file, \
    group_nucs_file, compare_consensus


def make_nuc_rows(consensus_seq: str) -> typing.List[typing.Tuple[str, dict]]:
    """ Make a set of nuc.csv rows to represent an expected consensus. """

    rows = []
    for consensus_nuc in consensus_seq:
        row = {nuc: '100' if nuc == consensus_nuc else '0'
               for nuc in 'ACGT'}
        row['coverage'] = '100'
        rows.append((consensus_nuc, row))
    return rows


# noinspection DuplicatedCode
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
        source_seqs = {('R1-seed', 'R1'): make_nuc_rows('ACACACGT')}
        target_seqs = {('R1-seed', 'R1'): make_nuc_rows('ACACACGG')}
        coverage_scores = [{'seed': 'R1-seed',
                            'region': 'R1',
                            'project': 'R1',
                            'on.score': '4'}]
        sample = Sample(MiseqRun(target_path='run1/Results/versionX'),
                        'sample42',
                        SampleFiles(region_consensus=source_seqs,
                                    coverage_scores=coverage_scores),
                        SampleFiles(region_consensus=target_seqs,
                                    coverage_scores=coverage_scores))
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
        source_seqs = {('R1-seed', 'R1'): make_nuc_rows('ACACAC')}
        target_seqs = {('R1-seed', 'R1'): make_nuc_rows('ACACAT')}
        coverage_scores = [{'seed': 'R1-seed',
                            'region': 'R1',
                            'project': 'R1',
                            'on.score': '4'}]
        sample = Sample(MiseqRun(target_path='run1/Results/versionX'),
                        'sample42',
                        SampleFiles(region_consensus=source_seqs,
                                    coverage_scores=coverage_scores),
                        SampleFiles(region_consensus=target_seqs,
                                    coverage_scores=coverage_scores))
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
            diffs,
            Scenarios.NONE,
            scenarios)

        self.assertEqual(expected_diffs, diffs)
        self.assertEqual(expected_scenarios, scenarios)

    def test_consensus_change_scenario(self):
        source_seqs = {('R1-seed', 'R1'): make_nuc_rows('ACACAC')}
        target_seqs = {('R1-seed', 'R1'): make_nuc_rows('ACACAT')}
        coverage_scores = [{'seed': 'R1-seed',
                            'region': 'R1',
                            'project': 'R1',
                            'on.score': '4'}]
        sample = Sample(MiseqRun(target_path='run1/Results/versionX'),
                        'sample42',
                        SampleFiles(region_consensus=source_seqs,
                                    coverage_scores=coverage_scores),
                        SampleFiles(region_consensus=target_seqs,
                                    coverage_scores=coverage_scores))
        expected_diffs = []
        expected_scenarios = {Scenarios.MAIN_CONSENSUS_CHANGED: ['.']}
        diffs = []
        scenarios = defaultdict(list)

        compare_consensus(
            sample,
            diffs,
            Scenarios.MAIN_CONSENSUS_CHANGED | Scenarios.OTHER_CONSENSUS_CHANGED,
            scenarios)

        self.assertEqual(expected_diffs, diffs)
        self.assertEqual(expected_scenarios, scenarios)

    def test_same_consensus(self):
        source_seqs = {('R1-seed', 'R1'): make_nuc_rows('ACACAC')}
        target_seqs = {('R1-seed', 'R1'): make_nuc_rows('ACACAC')}
        coverage_scores = [{'seed': 'R1-seed',
                            'region': 'R1',
                            'project': 'R1',
                            'on.score': '4'}]
        sample = Sample(MiseqRun(target_path='run1/Results/versionX'),
                        'sample42',
                        SampleFiles(region_consensus=source_seqs,
                                    coverage_scores=coverage_scores),
                        SampleFiles(region_consensus=target_seqs,
                                    coverage_scores=coverage_scores))
        expected_diffs = []
        expected_scenarios = {}
        diffs = []
        scenarios = defaultdict(list)

        compare_consensus(
            sample,
            diffs,
            Scenarios.MAIN_CONSENSUS_CHANGED | Scenarios.OTHER_CONSENSUS_CHANGED,
            scenarios)

        self.assertEqual(expected_diffs, diffs)
        self.assertEqual(expected_scenarios, scenarios)

    def test_one_consensus_changes(self):
        source_seqs = {('R1-seed', 'R1'): make_nuc_rows('ACACACGT'),
                       ('R2-seed', 'R2'): make_nuc_rows('ACACACGT')}
        target_seqs = {('R1-seed', 'R1'): make_nuc_rows('ACACACGT'),
                       ('R2-seed', 'R2'): make_nuc_rows('ACACAMGT')}
        coverage_scores = [{'seed': 'R1-seed',
                            'region': 'R1',
                            'project': 'R1',
                            'on.score': '4'},
                           {'seed': 'R2-seed',
                            'region': 'R2',
                            'project': 'R2',
                            'on.score': '4'}]
        sample = Sample(MiseqRun(target_path='run1/Results/versionX'),
                        'sample42',
                        SampleFiles(region_consensus=source_seqs,
                                    coverage_scores=coverage_scores),
                        SampleFiles(region_consensus=target_seqs,
                                    coverage_scores=coverage_scores))
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
            diffs,
            Scenarios.NONE,
            scenarios)

        self.assertEqual(expected_diffs, diffs)
        self.assertEqual(expected_consensus_distances, consensus_distances)

    def test_consensus_trailing_change(self):
        source_seqs = {('R1-seed', 'R1'): make_nuc_rows('ACTTAC------GTAC')}
        target_seqs = {('R1-seed', 'R1'): make_nuc_rows('ACTTAC')}
        coverage_scores = [{'seed': 'R1-seed',
                            'region': 'R1',
                            'project': 'R1',
                            'on.score': '4'}]
        sample = Sample(MiseqRun(target_path='run1/Results/versionX'),
                        'sample42',
                        SampleFiles(region_consensus=source_seqs,
                                    coverage_scores=coverage_scores),
                        SampleFiles(region_consensus=target_seqs,
                                    coverage_scores=coverage_scores))
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
            diffs,
            Scenarios.NONE,
            scenarios)

        self.assertEqual(expected_diffs, diffs)
        self.assertEqual(expected_consensus_distances, consensus_distances)

    def test_consensus_missing(self):
        source_seqs = {('R1-seed', 'R1'): make_nuc_rows('ACTTAC')}
        target_seqs = {}
        coverage_scores = [{'seed': 'R1-seed',
                            'region': 'R1',
                            'project': 'R1',
                            'on.score': '4'}]
        sample = Sample(MiseqRun(target_path='run1/Results/versionX'),
                        'sample42',
                        SampleFiles(region_consensus=source_seqs,
                                    coverage_scores=coverage_scores),
                        SampleFiles(region_consensus=target_seqs,
                                    coverage_scores=coverage_scores))
        expected_diffs = ['run1:sample42 consensus: R1-seed R1 MAX',
                          '- ACTTAC']
        expected_consensus_distances = []
        diffs = []
        scenarios = defaultdict(list)

        consensus_distances = compare_consensus(
            sample,
            diffs,
            Scenarios.NONE,
            scenarios)

        self.assertEqual(expected_diffs, diffs)
        self.assertEqual(expected_consensus_distances, consensus_distances)

    def test_consensus_added(self):
        source_seqs = {}
        target_seqs = {('R1-seed', 'R1'): make_nuc_rows('ACTTAC')}
        coverage_scores = [{'seed': 'R1-seed',
                            'region': 'R1',
                            'project': 'R1',
                            'on.score': '4'}]
        sample = Sample(MiseqRun(target_path='run1/Results/versionX'),
                        'sample42',
                        SampleFiles(region_consensus=source_seqs,
                                    coverage_scores=coverage_scores),
                        SampleFiles(region_consensus=target_seqs,
                                    coverage_scores=coverage_scores))
        expected_diffs = ['run1:sample42 consensus: R1-seed R1 MAX',
                          '+ ACTTAC']
        expected_consensus_distances = []
        diffs = []
        scenarios = defaultdict(list)

        consensus_distances = compare_consensus(
            sample,
            diffs,
            Scenarios.NONE,
            scenarios)

        self.assertEqual(expected_diffs, diffs)
        self.assertEqual(expected_consensus_distances, consensus_distances)


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

    # noinspection DuplicatedCode
    def test_group_nucs_single(self):
        output_file = StringIO("""\
sample,seed,region,query.nuc.pos,coverage
s1,R1-seed,R1,100,10
s1,R1-seed,R1,110,10
""")
        expected_groups = dict(s1={('R1-seed', 'R1'): [
            ('x', {'coverage': '10',
                   'query.nuc.pos': '100',
                   'region': 'R1',
                   'sample': 's1',
                   'seed': 'R1-seed'}),
            ('x', {'coverage': '10',
                   'query.nuc.pos': '110',
                   'region': 'R1',
                   'sample': 's1',
                   'seed': 'R1-seed'})]})

        groups = group_nucs_file(output_file)

        self.assertEqual(expected_groups, groups)

    def test_group_nucs_with_coverage(self):
        output_file = StringIO("""\
sample,seed,region,query.nuc.pos,A,C,G,T,del,ins,coverage
s1,R1-seed,R1,100,100,0,0,0,0,0,100
s1,R1-seed,R1,101,99,1,0,0,0,0,100
s1,R1-seed,R1,102,1,99,0,0,0,0,100
s1,R1-seed,R1,103,0,99,0,0,0,0,99
""")
        expected_groups = dict(s1={('R1-seed', 'R1'): [
            ('A', {'A': '100',
                   'C': '0',
                   'G': '0',
                   'T': '0',
                   'coverage': '100',
                   'del': '0',
                   'ins': '0',
                   'query.nuc.pos': '100',
                   'region': 'R1',
                   'sample': 's1',
                   'seed': 'R1-seed'}),
            ('A', {'A': '99',
                   'C': '1',
                   'G': '0',
                   'T': '0',
                   'coverage': '100',
                   'del': '0',
                   'ins': '0',
                   'query.nuc.pos': '101',
                   'region': 'R1',
                   'sample': 's1',
                   'seed': 'R1-seed'}),
            ('C', {'A': '1',
                   'C': '99',
                   'G': '0',
                   'T': '0',
                   'coverage': '100',
                   'del': '0',
                   'ins': '0',
                   'query.nuc.pos': '102',
                   'region': 'R1',
                   'sample': 's1',
                   'seed': 'R1-seed'}),
            ('x', {'A': '0',
                   'C': '99',
                   'G': '0',
                   'T': '0',
                   'coverage': '99',
                   'del': '0',
                   'ins': '0',
                   'query.nuc.pos': '103',
                   'region': 'R1',
                   'sample': 's1',
                   'seed': 'R1-seed'})]})

        groups = group_nucs_file(output_file)

        self.assertEqual(expected_groups, groups)

    # noinspection DuplicatedCode
    def test_group_nucs_two_samples(self):
        output_file = StringIO("""\
sample,seed,region,query.nuc.pos,coverage
s1,R1-seed,R1,100,10
s1,R1-seed,R1,110,10
s2,R1-seed,R1,101,10
s2,R1-seed,R1,111,10
""")
        expected_groups = dict(
            s1={('R1-seed', 'R1'): [('x', {'coverage': '10',
                                           'query.nuc.pos': '100',
                                           'region': 'R1',
                                           'sample': 's1',
                                           'seed': 'R1-seed'}),
                                    ('x', {'coverage': '10',
                                           'query.nuc.pos': '110',
                                           'region': 'R1',
                                           'sample': 's1',
                                           'seed': 'R1-seed'})]},
            s2={('R1-seed', 'R1'): [('x', {'coverage': '10',
                                           'query.nuc.pos': '101',
                                           'region': 'R1',
                                           'sample': 's2',
                                           'seed': 'R1-seed'}),
                                    ('x', {'coverage': '10',
                                           'query.nuc.pos': '111',
                                           'region': 'R1',
                                           'sample': 's2',
                                           'seed': 'R1-seed'})]})

        groups = group_nucs_file(output_file)

        self.assertEqual(expected_groups, groups)

    def test_group_nucs_two_regions(self):
        output_file = StringIO("""\
sample,seed,region,query.nuc.pos,coverage
s1,R1-seed,R1a,100,10
s1,R1-seed,R1a,110,10
s1,R1-seed,R1b,201,10
s1,R1-seed,R1b,211,10
""")
        expected_groups = dict(
            s1={('R1-seed', 'R1a'): [('x', {'coverage': '10',
                                            'query.nuc.pos': '100',
                                            'region': 'R1a',
                                            'sample': 's1',
                                            'seed': 'R1-seed'}),
                                     ('x', {'coverage': '10',
                                            'query.nuc.pos': '110',
                                            'region': 'R1a',
                                            'sample': 's1',
                                            'seed': 'R1-seed'})],
                ('R1-seed', 'R1b'): [('x', {'coverage': '10',
                                            'query.nuc.pos': '201',
                                            'region': 'R1b',
                                            'sample': 's1',
                                            'seed': 'R1-seed'}),
                                     ('x', {'coverage': '10',
                                            'query.nuc.pos': '211',
                                            'region': 'R1b',
                                            'sample': 's1',
                                            'seed': 'R1-seed'})]})

        groups = group_nucs_file(output_file)

        self.assertEqual(expected_groups, groups)

    def test_group_nucs_blank(self):
        output_file = StringIO("""\
sample,seed,region,query.nuc.pos,coverage
s1,R1-seed,R1,100,10
s1,R1-seed,R1,110,10
s1,R1-seed,R1,,0
""")
        expected_groups = dict(s1={('R1-seed', 'R1'): [
            ('x', {'coverage': '10',
                   'query.nuc.pos': '100',
                   'region': 'R1',
                   'sample': 's1',
                   'seed': 'R1-seed'}),
            ('x', {'coverage': '10',
                   'query.nuc.pos': '110',
                   'region': 'R1',
                   'sample': 's1',
                   'seed': 'R1-seed'}),
            ('x', {'coverage': '0',
                   'query.nuc.pos': '',
                   'region': 'R1',
                   'sample': 's1',
                   'seed': 'R1-seed'})]})

        groups = group_nucs_file(output_file)

        self.assertEqual(expected_groups, groups)
