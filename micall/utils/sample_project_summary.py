from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
import logging
from collections import defaultdict, Counter
from pathlib import Path

from micall.utils.sample_sheet_parser import sample_sheet_parser

logger = logging.getLogger(__name__)


def parse_args():
    parser = ArgumentParser(
        description="Report how many samples used each project.",
        formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('start_folder',
                        nargs='?',
                        default='/media/raw_data/MiSeq/runs',
                        help='a run folder, or a parent of many run folders',
                        type=Path)
    return parser.parse_args()


class Scanner:
    def __init__(self):
        self.project_counts = Counter()
        self.latest_dates = defaultdict(str)

    def process_run(self, run_folder: Path):
        if not (run_folder / 'needsprocessing').exists():
            return False
        if (run_folder / 'errorprocessing').exists():
            return True
        sample_sheet_path = run_folder / 'SampleSheet.csv'
        folder_name = run_folder.name
        with sample_sheet_path.open() as f:
            try:
                run_info = sample_sheet_parser(f)
            except Exception:
                raise RuntimeError(f'Failed to process run {folder_name}.')
        project_groups = defaultdict(list)
        for sample_info in run_info['DataSplit']:
            sample_number = sample_info['sample_number']
            project = sample_info['project']
            project_groups[sample_number].append(project)
        for project_names in project_groups.values():
            project_names.sort()
            name_tuple = tuple(project_names)
            self.project_counts[name_tuple] += 1
            self.latest_dates[name_tuple] = max(self.latest_dates[name_tuple], folder_name)
        top_projects = self.project_counts.most_common(3)
        summary = ', '.join(f'({", ".join(project_codes)}): {count}'
                            for project_codes, count in top_projects)
        logger.debug('After %s, top counts are: %s', folder_name, summary)
        return True

    def process_runs(self, runs_folder: Path):
        for file_path in sorted(runs_folder.iterdir()):
            if file_path.is_dir():
                # noinspection PyBroadException
                try:
                    self.process_run(file_path)
                except Exception:
                    logger.warning('Run %s failed.', file_path.name, exc_info=True)

    def report(self):
        for project_codes, count in self.project_counts.most_common():
            latest_date = self.latest_dates[project_codes]
            logger.info('%s: %d up to %s', ', '.join(project_codes), count, latest_date)
        logger.info('Total samples: %d', sum(self.project_counts.values()))


def main():
    logging.basicConfig(level=logging.DEBUG,
                        format='%(asctime)s[%(levelname)s]%(name)s: %(message)s')
    logger.info('Starting.')
    args = parse_args()
    scanner = Scanner()
    if not scanner.process_run(args.start_folder):
        scanner.process_runs(args.start_folder)
    scanner.report()
    logger.info('Done.')

if __name__ == '__main__':
    main()
