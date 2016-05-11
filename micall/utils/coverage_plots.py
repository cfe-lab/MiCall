import os
import argparse
import itertools
from collections import Counter
from csv import DictReader, DictWriter
from matplotlib import pyplot as plt, patches
from operator import itemgetter

from micall.core import project_config, aln2counts


def coverage_plot(amino_csv, coverage_scores_csv, path_prefix=None):
    """ Generate coverage plots.

    @param amino_csv: an open file object that holds amino acid frequencies
    @param coverage_scores_csv: an open file object to write the coverage scores
    @param path_prefix: path prefix for coverage maps. Defaults to the path
    of amino_csv without the extension.
    @return: a list of full paths to the image files.
    """
    # imports project information from JSON
    if path_prefix is not None:
        path, prefix = os.path.split(path_prefix)
    else:
        path, filename = os.path.split(amino_csv.name)
        prefix = filename.split('.')[0]
    projects = project_config.ProjectConfig.loadScoring()
    reader = DictReader(amino_csv)
    writer = DictWriter(coverage_scores_csv,
                        ['project',
                         'region',
                         'seed',
                         'q.cut',
                         'min.coverage',
                         'which.key.pos',
                         'off.score',
                         'on.score'],
                        lineterminator=os.linesep)
    writer.writeheader()
    paths = []

    for (seed, region), group in itertools.groupby(reader, itemgetter('seed',
                                                                      'region')):
        counts = Counter()
        for row in group:
            pos = int(row['refseq.aa.pos'])
            total = sum([int(row[aa]) for aa in aln2counts.AMINO_ALPHABET])
            counts[pos] = total
            qcut = row['q-cutoff']
        # use region to retrieve coordinate reference
        for project_region in projects.getProjectRegions(seed, region):
            min_coverage = None
            min_coverage_pos = None
            project_name = project_region['project_name']
            region_length = project_region['coordinate_region_length']
            x = range(1, region_length+1)
            y = [counts[pos] for pos in x]

            _fig, ax = plt.subplots(figsize=(4, 3), dpi=100)
            key_positions = project_region['key_positions']
            if not key_positions:
                key_positions.append(dict(start_pos=1,
                                          end_pos=region_length))
            for key_pos in key_positions:
                start = key_pos['start_pos']
                end = key_pos['end_pos']
                if end is None:
                    end = start
                for pos in range(start, end+1):
                    count = counts[pos]
                    if min_coverage is None or count < min_coverage:
                        min_coverage = count
                        min_coverage_pos = pos
                start -= 0.5
                end += 0.5
                ax.add_patch(patches.Rectangle(xy=(start, 50),
                                               width=end-start,
                                               height=150,
                                               fc='grey',
                                               ec='grey'))
            if min_coverage <= project_region['min_coverage1']:
                coverage_score_on = 1
            elif min_coverage <= project_region['min_coverage2']:
                coverage_score_on = 2
            elif min_coverage <= project_region['min_coverage3']:
                coverage_score_on = 3
            else:
                coverage_score_on = 4
            max_coverage = max(y)
            if max_coverage == 0:
                coverage_score_off = 0
            elif max_coverage <= 10:
                coverage_score_off = -1
            elif max_coverage <= 100:
                coverage_score_off = -2
            else:
                coverage_score_off = -3
            plt.step(x, y, linewidth=2, where='mid')
            plt.xlim([1, region_length])
            plt.ylim([0.1, 1000000])
            plt.yscale('log')
            plt.xlabel('Reference coordinates (AA)', fontsize=9)
            plt.ylabel('Coverage', fontsize=9)
            plt.tight_layout()
            figname = '{}.{}.{}.png'.format(prefix, project_name, region)
            dest = os.path.join(path, figname)
            paths.append(dest)
            plt.savefig(dest)  # write image to file
            plt.clf()
            row = {'project': project_name,
                   'region': region,
                   'seed': seed,
                   'q.cut': qcut,
                   'min.coverage': min_coverage,
                   'which.key.pos': min_coverage_pos,
                   'off.score': coverage_score_off,
                   'on.score': coverage_score_on}
            writer.writerow(row)

    return paths  # locations of image files


def parse_args():
    parser = argparse.ArgumentParser(description='Generate coverage plots from MiCall outputs.')
    parser.add_argument('amino_csv', type=argparse.FileType('rU'),
                        help='<input> CSV containing amino acid frequency outputs.')
    return parser.parse_args()


def main():
    args = parse_args()
    coverage_plot(amino_csv=args.amino_csv)

if __name__ == '__main__':
    # note, must be called from project root if executing directly
    # i.e., python micall/g2p/sam_g2p.py -h
    main()
elif __name__ == '__live_coding__':
    import unittest
    from micall.tests.coverage_plots_test import CoveragePlotsTest

    suite = unittest.TestSuite()
    suite.addTest(CoveragePlotsTest("test_simple"))
    test_results = unittest.TextTestRunner().run(suite)

    print(test_results.errors)
    print(test_results.failures)
