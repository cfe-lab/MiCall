from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
from csv import DictReader
from pathlib import Path

from micall.monitor.kive_watcher import get_version_key


def parse_args():
    parser = ArgumentParser(description='Find samples with good coverage.',
                            formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('--runs',
                        default='/media/raw_data/MiSeq/runs',
                        type=Path,
                        help='Runs folder to scan for samples.')
    return parser.parse_args()


def main():
    args = parse_args()
    target_seeds = ('HCV-6n',)
    best_coverage = 1
    runs: Path = args.runs
    for i, flag_path in enumerate(sorted(runs.glob('*/needsprocessing'),
                                         reverse=True)):
        run_path: Path = flag_path.parent
        try:
            versions = sorted((run_path / "Results").glob('version_*'),
                              key=get_version_key)
        except FileNotFoundError:
            continue
        if not versions:
            continue
        latest_version = versions[-1]
        coverage_path = latest_version / "coverage_scores.csv"
        if not coverage_path.exists():
            print('Missing file:', coverage_path)
            continue
        good_rows = []
        with coverage_path.open() as f:
            reader = DictReader(f)
            for row in reader:
                seed = row['seed']
                if seed not in target_seeds:
                    continue
                coverage = int(row['on.score'])
                if coverage >= best_coverage:
                    best_coverage = coverage
                    good_rows.append(row)
        if good_rows or i % 100 == 0:
            print(run_path.name, latest_version.name)
            for row in good_rows:
                print('*',
                      row['sample'],
                      row['region'],
                      row['seed'],
                      row['on.score'],
                      row['min.coverage'],
                      row['which.key.pos'])


if __name__ == '__main__':
    main()
