import os
import argparse
import errno
import itertools
from collections import Counter
from csv import DictReader, DictWriter
from matplotlib import pyplot as plt, patches
from matplotlib.ticker import ScalarFormatter
from operator import itemgetter
import tarfile

from micall.core import project_config, aln2counts


def coverage_plot(amino_csv,
                  coverage_scores_csv,
                  coverage_maps_path=None,
                  coverage_maps_prefix=None,
                  filetype='png'):
    """ Generate coverage plots.

    @param amino_csv: an open file object that holds amino acid frequencies
    @param coverage_scores_csv: an open file object to write the coverage scores
    @param coverage_maps_path: path for coverage maps. Defaults to the path
    of amino_csv.
    @param coverage_maps_prefix: file name prefix for coverage maps. Full name
    will be "path/prefix.SEED.REGION.filetype" or "path/SEED.REGION.filetype"
    if there is no prefix.
    @param filetype: controls which file type will be saved, must be supported
    by matplotlib (probably png, pdf, ps, eps and svg).
    @return: a list of full paths to the image files.
    """
    # imports project information from JSON
    if coverage_maps_path is None:
        coverage_maps_path, _ = os.path.split(amino_csv.name)
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

    MAX_COVERAGE = 1000000
    axis_formatter = ScalarFormatter()
    axis_formatter.set_powerlimits((-8, 8))
    _fig, ax = plt.subplots(figsize=(4, 3), dpi=100)
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
            left_margin = -region_length / 50.0
            plt.xlim([left_margin, region_length])
            plt.ylim([0.5, MAX_COVERAGE])
            plt.yscale('log')
            ax.yaxis.set_major_formatter(axis_formatter)
            ax.add_patch(patches.Rectangle(xy=(left_margin, 0),
                                           width=-left_margin,
                                           height=10,
                                           fc='black',
                                           ec='black'))
            ax.add_patch(patches.Rectangle(xy=(left_margin, 10),
                                           width=-left_margin,
                                           height=40,
                                           fc='red',
                                           ec='red'))
            ax.add_patch(patches.Rectangle(xy=(left_margin, 50),
                                           width=-left_margin,
                                           height=50,
                                           fc='yellow',
                                           ec='yellow'))
            ax.add_patch(patches.Rectangle(xy=(left_margin, 50),
                                           width=-left_margin,
                                           height=50,
                                           fc='yellow',
                                           ec='yellow'))
            ax.add_patch(patches.Rectangle(xy=(left_margin, 100),
                                           width=-left_margin,
                                           height=MAX_COVERAGE-100,
                                           fc='lightgreen',
                                           ec='lightgreen'))
            plt.hlines(100, 1, region_length, linestyles='dashed')
            plt.xlabel('Reference coordinates (AA)', fontsize=9)
            plt.ylabel('Coverage', fontsize=9)
            plt.tight_layout()
            figname_parts = [project_name, region, filetype]
            if coverage_maps_prefix:
                figname_parts.insert(0, coverage_maps_prefix)
            figname = '.'.join(figname_parts)
            dest = os.path.join(coverage_maps_path, figname)
            paths.append(dest)
            plt.savefig(dest)  # write image to file
            plt.cla()  # clear the axis, but don't remove the axis itself.
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


def make_tar_path(tar_path):
    try:
        os.makedirs(tar_path)
    except OSError as exc:
        if exc.errno != errno.EEXIST:
            raise


def parse_args():
    parser = argparse.ArgumentParser(description='Generate coverage plots from MiCall outputs.')
    parser.add_argument('amino_csv', type=argparse.FileType('rU'),
                        help='<input> CSV containing amino acid frequency outputs.')
    parser.add_argument('coverage_scores_csv', type=argparse.FileType('wb'),
                        help='<output> CSV coverage scores.')
    parser.add_argument('coverage_maps_tar',
                        type=argparse.FileType('wb'),
                        help='<output> tar file of coverage maps.')
    return parser.parse_args()


def main():
    args = parse_args()
    coverage_maps_path, _ = os.path.splitext(args.coverage_maps_tar.name)
    make_tar_path(coverage_maps_path)
    coverage_plot(amino_csv=args.amino_csv,
                  coverage_scores_csv=args.coverage_scores_csv,
                  coverage_maps_path=coverage_maps_path)
    with tarfile.open(fileobj=args.coverage_maps_tar, mode='w') as tar:
        for image_name in os.listdir(coverage_maps_path):
            image_path = os.path.join(coverage_maps_path, image_name)
            tar.add(image_path)

if __name__ == '__main__':
    # note, must be called from project root if executing directly
    # i.e., python micall/g2p/sam_g2p.py -h
    main()
elif __name__ == '__live_coding__':
    is_unit_test = False
    if is_unit_test:
        import unittest
        from micall.tests.coverage_plots_test import CoveragePlotsTest

        suite = unittest.TestSuite()
        suite.addTest(CoveragePlotsTest("test_simple"))
        test_results = unittest.TextTestRunner().run(suite)

        print(test_results.errors)
        print(test_results.failures)
    else:
        amino_path = '../tests/working/1234A-V3LOOP_S1.amino.csv'
        coverage_scores_path = '../tests/working/1234A-V3LOOP_S1.coverage_scores.csv'
        coverage_maps_path = '../tests/working/coverage_maps'
        coverage_maps_prefix = '1234A-V3LOOP_S1'
        make_tar_path(coverage_maps_path)

        with open(amino_path, 'rU') as amino_csv, \
                open(coverage_scores_path, 'wb') as coverage_scores_csv:
            coverage_plot(amino_csv, coverage_scores_csv, coverage_maps_path, coverage_maps_prefix)
