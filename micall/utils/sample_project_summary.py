from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
import logging
from collections import defaultdict, Counter
from pathlib import Path

from micall.utils.check_sample_sheet import check_sample_name_consistency
from micall.utils.list_fastq_files import list_fastq_file_names
from micall.utils.sample_sheet_parser import UnknownSamplesInOverrides, read_sample_sheet_and_overrides

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
    verbosity_group = parser.add_mutually_exclusive_group()
    verbosity_group.add_argument('--verbose', action='store_true', help='Increase output verbosity.')
    verbosity_group.add_argument('--no-verbose', action='store_true', help='Normal output verbosity.', default=True)
    verbosity_group.add_argument('--debug', action='store_true', help='Maximum output verbosity.')
    verbosity_group.add_argument('--quiet', action='store_true', help='Minimize output verbosity.')
    return parser.parse_args()


class Scanner:
    def __init__(self):
        self.project_counts: Counter[tuple[str, ...]] = Counter()
        self.latest_dates = defaultdict(str)

    def process_run(self, run_folder: Path):
        if not (run_folder / 'needsprocessing').exists():
            return False
        if (run_folder / 'errorprocessing').exists():
            return True
        sample_sheet_path = run_folder / 'SampleSheet.csv'
        folder_name = run_folder.name
        try:
            run_info = read_sample_sheet_and_overrides(sample_sheet_path)
        except UnknownSamplesInOverrides as e:
            samples = '\n\t'.join(e.samples)
            logger.error('Run %s has unknown samples in overrides:\n\t%s', folder_name, samples)
            return True
        except Exception:
            raise RuntimeError(f'Failed to process run {folder_name}.')
        file_names = list_fastq_file_names(run_folder, "*_R1_*.fastq.gz", fallback_to_run_path=False)
        check_sample_name_consistency(sample_sheet_path, file_names, run_folder)
        project_groups = defaultdict(list)
        data_split = run_info.get('DataSplit')
        if data_split is None:
            raise RuntimeError(f'No DataSplit section in sample sheet for run {folder_name}.')
        assert hasattr(data_split, '__iter__'), \
            f'DataSplit section is not iterable in sample sheet for run {folder_name}.'
        for sample_info in data_split:
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
            print('%s: %d up to %s' % (', '.join(project_codes), count, latest_date))
        print('Total samples:', sum(self.project_counts.values()))


def main():
    args = parse_args()
    if args.quiet:
        logger.setLevel(logging.ERROR)
    elif args.verbose:
        logger.setLevel(logging.INFO)
    elif args.debug:
        logger.setLevel(logging.DEBUG)
    else:
        logger.setLevel(logging.WARN)

    logging.basicConfig(level=logger.level,
                        format='%(asctime)s[%(levelname)s]%(name)s: %(message)s')

    logger.info('Starting.')
    scanner = Scanner()
    if not scanner.process_run(args.start_folder):
        scanner.process_runs(args.start_folder)
    scanner.report()
    logger.info('Done.')

if __name__ == '__main__':
    main()
