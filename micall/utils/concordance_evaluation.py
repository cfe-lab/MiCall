import os
import csv
from collections import defaultdict
import argparse
from numpy import std
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt   # noqa


def plot_histo(fig, counts, xlabel, fig_path, min_number=10):
    if len(counts) < min_number:
        return None, None
    plt.hist(x=counts, bins='auto', color='#0504aa', alpha=0.7, rwidth=0.85)
    plt.xlabel(xlabel)
    num_samples = len(counts)
    plt.ylabel(f'Frequency (total {num_samples})')
    average = sum(counts) / len(counts)
    standarddev = std(counts)
    plt.figtext(0.2, 0.9, f"Average: {average:.2f}, Standard deviation {standarddev:.2f}")
    plt.savefig(fig_path)
    plt.cla()
    fig.texts = []  # remove figure text
    return average, standarddev


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('results_folder', help="Folder in which the result folders live")
    parser.add_argument('subfolder', help="Subfolder in which the concordance files live")
    parser.add_argument('--denovo', default=None, help="Subfolder in which the denovo concordance files live")
    parser.add_argument('--min_counts', default=10, help="Minimum number of data points to create a histogram")
    args = parser.parse_args()

    result_directories = next(os.walk(args.results_folder))[1]
    folders_to_read = [os.path.join(args.results_folder, result_dir, args.subfolder)
                       for result_dir in result_directories]
    if args.denovo is not None:
        denovo_folders = [os.path.join(args.results_folder, result_dir, args.denovo)
                          for result_dir in result_directories]
        folders_to_read = folders_to_read + denovo_folders

    fig, ax = plt.subplots(figsize=(4, 3), dpi=100)
    plt.subplots_adjust(bottom=0.2, left=0.2)

    # {seed: {region: [values]}}
    all_seed_stats = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: [])))
    # {coord: {region: [values]}}
    all_coord_stats = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: [])))

    for folder in folders_to_read:
        try:
            with open(os.path.join(folder, 'concordance_seed.csv')) as concordance_seed_csv, \
                    open(os.path.join(folder, 'concordance.csv')) as concordance_csv:
                concordance_reader = csv.DictReader(concordance_csv)
                concordance_seed_reader = csv.DictReader(concordance_seed_csv)
                for row in concordance_reader:
                    concordance = float(row['pct_concordance'])
                    coverage = float(row['pct_covered'])
                    reference = row['reference']
                    region = row['region']
                    if coverage > 0:
                        all_coord_stats[reference][region]['concordance'].append(concordance)
                        all_coord_stats[reference][region]['coverage'].append(coverage)
                        all_coord_stats[reference][region]['covconc'].append(coverage*concordance / 100)
                for row in concordance_seed_reader:
                    concordance = float(row['pct_concordance'])
                    coverage = float(row['pct_covered'])
                    seed = row['seed_name']
                    region = row['region']
                    if coverage > 0:
                        all_seed_stats[seed][region]['concordance'].append(concordance)
                        all_seed_stats[seed][region]['coverage'].append(coverage)
                        all_seed_stats[seed][region]['covconc'].append(coverage * concordance / 100)
        except FileNotFoundError:
            print(f"No concordance files found in {folder}")


    with open(os.path.join(args.results_folder, 'concordance_stats.csv'), 'w') as stats_csv:
        columns = ['reference', 'region', 'average_concordance', 'std_concordance',
                   'average_covered', 'std_covered', 'average_covconc', 'std_covconc']
        writer = csv.DictWriter(stats_csv, columns, lineterminator=os.linesep)
        for coord, regions in all_coord_stats.items():
            for region, stats in regions.items():
                counts = stats['concordance']
                xlabel = 'Concordance'
                fig_path = os.path.join(args.results_folder, f'{coord}.{region}.coord_concordance.png')
                av_concordance, std_concordance = plot_histo(fig, counts, xlabel, fig_path, min_number=args.min_counts)
                counts = stats['coverage']
                xlabel = '% Covered'
                fig_path = os.path.join(args.results_folder, f'{coord}.{region}.coord_covered.png')
                av_coverage, std_coverage = plot_histo(fig, counts, xlabel, fig_path, min_number=args.min_counts)
                counts = stats['covconc']
                xlabel = '% Covered and Concordant'
                fig_path = os.path.join(args.results_folder, f'{coord}.{region}.coord_coveredconcordant.png')
                av_convconc, std_covconc = plot_histo(fig, counts, xlabel, fig_path, min_number=args.min_counts)
                if av_concordance is not None:  # then the others are not None either
                    row = dict(reference=coord,
                               region=region,
                               average_concordance=av_concordance,
                               std_concordance=std_concordance,
                               average_covered=av_coverage,
                               std_covered=std_coverage,
                               average_covconc=av_convconc,
                               std_covconc=std_covconc)
                    writer.writerow(row)

    with open(os.path.join(args.results_folder, 'concordance_stats_seed.csv'), 'w') as stats_seed_csv:
        columns_seed = ['seed', 'region', 'average_concordance', 'std_concordance',
                        'average_covered', 'std_covered', 'average_covconc', 'std_covconc']
        writer = csv.DictWriter(stats_seed_csv, columns_seed, lineterminator=os.linesep)
        for seed, regions in all_seed_stats.items():
            for region, stats in regions.items():
                counts = stats['concordance']
                xlabel = 'Concordance'
                fig_path = os.path.join(args.results_folder, f'{seed}.{region}.seed_concordance.png')
                av_concordance, std_concordance = plot_histo(fig, counts, xlabel, fig_path, min_number=args.min_counts)
                counts = stats['coverage']
                xlabel = '% Covered'
                fig_path = os.path.join(args.results_folder, f'{seed}.{region}.seed_covered.png')
                av_coverage, std_coverage = plot_histo(fig, counts, xlabel, fig_path, min_number=args.min_counts)
                counts = stats['covconc']
                xlabel = '% Covered and Concordant'
                fig_path = os.path.join(args.results_folder, f'{seed}.{region}.seed_coveredconcordant.png')
                av_convconc, std_covconc = plot_histo(fig, counts, xlabel, fig_path, min_number=args.min_counts)
                if av_concordance is not None:  # then the others are not None either
                    row = dict(seed=seed,
                               region=region,
                               average_concordance=av_concordance,
                               std_concordance=std_concordance,
                               average_covered=av_coverage,
                               std_covered=std_coverage,
                               average_covconc=av_convconc,
                               std_covconc=std_covconc)
                    writer.writerow(row)


if __name__ == '__main__':
    main()
