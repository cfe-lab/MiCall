import re
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
import logging
from csv import DictReader
from pathlib import Path

from micall.utils.sample_sheet_parser import read_sample_sheet_and_overrides

logger = logging.getLogger(__name__)


def parse_args():
    parser = ArgumentParser(
        description="Look for samples that didn't get processed.",
        formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('start_folder',
                        nargs='?',
                        default='/media/raw_data/MiSeq/runs',
                        help='a run folder, or a parent of many run folders',
                        type=Path)
    parser.add_argument('--skip_mid_hcv',
                        action='store_true',
                        help="Don't report missing samples with the MidHCV project.")
    return parser.parse_args()


def process_run(run_folder: Path, skip_mid_hcv: bool):
    if not (run_folder / 'needsprocessing').exists():
        return False
    if (run_folder / 'errorprocessing').exists():
        return True
    sample_sheet_path = run_folder / 'SampleSheet.csv'
    try:
        run_info = read_sample_sheet_and_overrides(sample_sheet_path)
    except Exception:
        raise RuntimeError(f'Failed to process run {run_folder.name}.')
    sample_names = set(run_info['Data'])
    if skip_mid_hcv:
        sample_names = {sample_name
                        for sample_name in sample_names
                        if not re.match(r'.*MidHCV_S\d+$', sample_name)}
    cascade_path = run_folder / 'Results' / 'version_7.9' / 'cascade.csv'
    with cascade_path.open() as f:
        reader = DictReader(f)
        cascade_samples = {row['sample'] for row in reader}
    missing_samples = sample_names - cascade_samples
    if missing_samples:
        logger.error('Missing samples in run %s: %s',
                     run_folder.name,
                     sorted(missing_samples))
    return True


def process_runs(runs_folder: Path, skip_mid_hcv: bool):
    for file_path in sorted(runs_folder.iterdir()):
        if file_path.is_dir():
            # noinspection PyBroadException
            try:
                process_run(file_path, skip_mid_hcv)
            except Exception:
                logger.warning('Run %s failed.', file_path.name, exc_info=True)


def main():
    logging.basicConfig(level=logging.INFO,
                        format='%(asctime)s[%(levelname)s]%(name)s: %(message)s')
    logger.info('Starting.')
    args = parse_args()
    if not process_run(args.start_folder, args.skip_mid_hcv):
        process_runs(args.start_folder, args.skip_mid_hcv)
    logger.info('Done.')


if __name__ == '__main__':
    main()
