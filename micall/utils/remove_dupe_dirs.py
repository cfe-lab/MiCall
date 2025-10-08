import shutil
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
from pathlib import Path


def parse_args():
    parser = ArgumentParser(
        description='Remove folders that have already been zipped.',
        formatter_class=ArgumentDefaultsHelpFormatter)
    # noinspection PyTypeChecker
    parser.add_argument(
        'raw_data',
        type=Path,
        nargs='?',
        default=str(Path.home() / 'data' / 'RAW_DATA'),
        help='Raw data folder to scan for duplicate results folders.')
    parser.add_argument('--force', '-f',
                        action='store_true',
                        help="Don't ask for confirmation before deleting.")
    return parser.parse_args()


def main():
    args = parse_args()
    duplicate_folders = []
    runs_path: Path = args.raw_data / 'MiSeq' / 'runs'
    for zip_path in runs_path.glob('*/Results/*.zip'):
        results_folder = zip_path.parent
        original_folder = results_folder / zip_path.stem
        rel_path = original_folder.relative_to(runs_path)
        if original_folder.is_dir():
            print(rel_path)
            duplicate_folders.append(original_folder)

    if not duplicate_folders:
        print('No duplicates found.')
        return

    if not args.force:
        confirmation = input(f'Are you sure you want to delete '
                             f'{len(duplicate_folders)} folders? Y/[N] ')
        if confirmation.upper() != 'Y':
            exit('Aborted.')
    duplicate_folders.sort()
    for original_folder in duplicate_folders:
        rel_path = original_folder.relative_to(runs_path)
        print(f'deleting {rel_path}...')
        shutil.rmtree(original_folder)

if __name__ == '__main__':
    main()
