#!/usr/bin/env python3.4

import os
import argparse
from collections import Counter
from csv import DictReader, DictWriter
import errno
import itertools
from operator import itemgetter
import tarfile

from matplotlib.ticker import FuncFormatter

from micall.core import project_config, aln2counts

# NOTE: this must be performed BEFORE pyplot is imported
# http://stackoverflow.com/a/3054314/4794
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt, patches  # @IgnorePep8


def coverage_plot(amino_csv,
                  coverage_scores_csv,
                  coverage_maps_path=None,
                  coverage_maps_prefix=None,
                  filetype='png',
                  excluded_projects=None):
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
    @param excluded_projects: a list of project names to exclude
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
    fontsize = 8
    axis_formatter = FuncFormatter(lambda x, p: format(int(x), ','))
    _fig, ax = plt.subplots(figsize=(4, 3), dpi=100)
    for (seed, region), group in itertools.groupby(reader, itemgetter('seed',
                                                                      'region')):
        coverage_counts = Counter()
        stop_counts = Counter()
        low_quality_counts = Counter()
        deletion_counts = Counter()
        partial_counts = Counter()
        insertion_counts = Counter()
        clipping_counts = Counter()
        for row in group:
            pos = int(row['refseq.aa.pos'])
            deletion_count = int(row['del'])
            total = sum([int(row[aa]) for aa in aln2counts.AMINO_ALPHABET])
            total += deletion_count
            coverage_counts[pos] = total
            stop_counts[pos] = int(row['*'])
            low_quality_counts[pos] = int(row['X'])
            deletion_counts[pos] = deletion_count
            partial_counts[pos] = int(row['partial'])
            insertion_counts[pos] = int(row['ins'])
            clipping_counts[pos] = int(row['clip'])
            qcut = row['q-cutoff']
        # use region to retrieve coordinate reference
        for project_region in projects.getProjectRegions(
                seed,
                region,
                excluded_projects=excluded_projects):
            min_coverage = None
            min_coverage_pos = None
            project_name = project_region['project_name']
            region_length = project_region['coordinate_region_length']
            x = range(1, region_length+1)
            y_coverage = [coverage_counts[pos] for pos in x]
            y_stops = [stop_counts[pos] for pos in x]
            y_low_quality = [low_quality_counts[pos] for pos in x]
            y_deletions = [deletion_counts[pos] for pos in x]
            y_partials = [partial_counts[pos] for pos in x]
            y_insertions = [insertion_counts[pos] for pos in x]
            y_clipping = [clipping_counts[pos] for pos in x]

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
                    count = coverage_counts[pos]
                    if min_coverage is None or count < min_coverage:
                        min_coverage = count
                        min_coverage_pos = pos
                start -= 0.5
                end += 0.5
                ax.add_patch(patches.Rectangle(xy=(start, 50),
                                               width=end-start,
                                               height=150,
                                               fc='black',
                                               ec='grey',
                                               zorder=50,
                                               alpha=.5))
            if min_coverage <= project_region['min_coverage1']:
                coverage_score_on = 1
            elif min_coverage <= project_region['min_coverage2']:
                coverage_score_on = 2
            elif min_coverage <= project_region['min_coverage3']:
                coverage_score_on = 3
            else:
                coverage_score_on = 4
            max_coverage = max(y_coverage)
            if max_coverage == 0:
                coverage_score_off = 0
            elif max_coverage <= 10:
                coverage_score_off = -1
            elif max_coverage <= 100:
                coverage_score_off = -2
            else:
                coverage_score_off = -3
            plt.step(x, y_coverage, linewidth=2, where='mid', label='coverage', zorder=100)
            left_margin = -region_length / 25.0
            plt.xlim([left_margin, region_length])
            plt.ylim([0.5, MAX_COVERAGE])
            plt.yscale('log')
            ax.yaxis.set_major_formatter(axis_formatter)
            plt.tick_params(axis='both', labelsize=fontsize)
            ax.add_patch(patches.Rectangle(xy=(left_margin*0.5, 0),
                                           width=-left_margin*0.4,
                                           height=10,
                                           fc='black',
                                           ec='black'))
            ax.add_patch(patches.Rectangle(xy=(left_margin*0.5, 10),
                                           width=-left_margin*0.4,
                                           height=40,
                                           fc='red',
                                           ec='red'))
            ax.add_patch(patches.Rectangle(xy=(left_margin*0.5, 50),
                                           width=-left_margin*0.4,
                                           height=50,
                                           fc='yellow',
                                           ec='yellow'))
            ax.add_patch(patches.Rectangle(xy=(left_margin*0.5, 100),
                                           width=-left_margin*0.4,
                                           height=MAX_COVERAGE-100,
                                           fc='lightgreen',
                                           ec='lightgreen'))
            plt.plot((1, region_length), (100, 100), 'k--', zorder=51)
            plt.xlabel('Reference coordinates (AA)', fontsize=9)
            plt.ylabel('Read count', fontsize=9)
            plt.tight_layout()
            figname_parts = [project_name, region, filetype]
            if coverage_maps_prefix:
                figname_parts.insert(0, coverage_maps_prefix)
            paths.append(save_figure(coverage_maps_path, figname_parts))
            plt.step(x, y_deletions, where='mid', label='deletions', zorder=99)
            plt.step(x, y_stops, where='mid', label='stop codons', zorder=98)
            plt.step(x, y_partials, where='mid', label='partial dels', zorder=97)
            plt.step(x, y_clipping, where='mid', label='soft clipped', zorder=96)
            plt.step(x, y_insertions, where='mid', label='insertions', zorder=95)
            plt.step(x, y_low_quality, where='mid', label='low quality', zorder=94)
            plt.legend(loc='best', fontsize=fontsize, fancybox=True)
            figname_parts.insert(-1, 'details')
            paths.append(save_figure(coverage_maps_path, figname_parts))
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


def save_figure(coverage_maps_path, figname_parts):
    """ Write the current figure to a file.

    :param str coverage_maps_path: the folder to write the file in
    :param list figname_parts: will be joined together with dots to make the
        file name
    :return: the file name it was written to
    """
    figname = '.'.join(figname_parts)
    dest = os.path.join(coverage_maps_path, figname)
    plt.savefig(dest)  # write image to file
    return dest


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
            archive_path = os.path.join('coverage_maps', image_name)
            tar.add(image_path, archive_path)

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
        amino_path = '../tests/working/2110A-V3LOOP_S13.amino.csv'
        coverage_scores_path = '../tests/working/2110A-V3LOOP_S13.coverage_scores.csv'
        coverage_maps_path = '../tests/working/coverage_maps'
        coverage_maps_prefix = '2110A-V3LOOP_S13'
        make_tar_path(coverage_maps_path)

        with open(amino_path, 'rU') as amino_csv, \
                open(coverage_scores_path, 'wb') as coverage_scores_csv:
            coverage_plot(amino_csv, coverage_scores_csv, coverage_maps_path, coverage_maps_prefix)
