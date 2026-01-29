import os
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
from collections import Counter, defaultdict
from csv import DictWriter
from pathlib import Path
from tempfile import NamedTemporaryFile
from zipfile import ZipFile, ZIP_DEFLATED

import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd

from micall.utils.list_fastq_files import find_fastq_source_folder


def parse_args():
    parser = ArgumentParser(description='Scan run folders and report disk usage.',
                            formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('run_sizes_csv',
                        nargs='?',
                        default='run_sizes.csv',
                        type=Path,
                        help='CSV file to write the size of each run in MB')
    parser.add_argument('--runs',
                        default='/media/raw_data/MiSeq/runs',
                        type=Path,
                        help='path to MiSeq run folders')
    parser.add_argument('--group_size',
                        default=50,
                        type=int,
                        help='scan the first run from each group')
    return parser.parse_args()


def format_size(size):
    return f'{int(size/1024/1024+0.5)}'


def zip_factory():
    return ZipFile(NamedTemporaryFile(), 'w', ZIP_DEFLATED)


def scan_run_folders(runs_folder: Path, run_sizes_csv, group_size: int):
    writer = DictWriter(run_sizes_csv, ['run', 'version', 'fastqs', 'outputs', 'zipped'])
    writer.writeheader()
    run_folders = sorted(run_folder
                         for run_folder in runs_folder.iterdir()
                         if run_folder.is_dir())
    chosen_folders = (run_folder
                      for i, run_folder in enumerate(run_folders)
                      if i % group_size == 0)
    for run_folder in chosen_folders:
        print(run_folder)
        data_size = 0
        version_sizes = Counter()
        version_zips = defaultdict(zip_factory)
        results_path = run_folder / 'Results'
        # Find the actual FASTQ folder (could be BaseCalls or Alignment_*/*/Fastq)
        fastq_folder = find_fastq_source_folder(run_folder, '*_R1_*')
        fastq_folder_str = str(fastq_folder) if fastq_folder else None
        # noinspection PyTypeChecker
        for dirpath, dirnames, filenames in os.walk(run_folder):
            dirpath = Path(dirpath)
            for file_name in filenames:
                if file_name.endswith('.DS_Store'):
                    continue
                file_path = dirpath / file_name
                if (file_name.endswith('.fastq.gz') and
                        fastq_folder_str and
                        str(file_path).startswith(fastq_folder_str)):
                    data_size += file_path.stat().st_size
                else:
                    try:
                        rel_path = file_path.relative_to(results_path)
                        parents = list(rel_path.parents)
                        if len(parents) >= 2:
                            version = parents[-2]
                            entry_path = rel_path.relative_to(version)
                        else:
                            entry_path = version = Path(file_name)

                        version_sizes[version] += file_path.stat().st_size
                        version_zip: ZipFile = version_zips[version]
                        version_zip.write(file_path, entry_path)
                    except ValueError:
                        pass
        if not version_sizes:
            version_sizes[None] = 0
        for version, size in sorted(version_sizes.items()):
            version_zip = version_zips[version]
            raw_file = version_zip.fp
            version_zip.close()  # Writes index.
            zip_size = raw_file.tell()
            raw_file.close()
            writer.writerow(dict(run=run_folder.name,
                                 version=version,
                                 fastqs=format_size(data_size),
                                 outputs=format_size(size),
                                 zipped=format_size(zip_size)))
        run_sizes_csv.flush()


def main():
    args = parse_args()
    run_sizes_path: Path = args.run_sizes_csv
    if not run_sizes_path.exists():
        with run_sizes_path.open('w') as run_sizes_csv:
            scan_run_folders(args.runs, run_sizes_csv, args.group_size)

    runs = pd.read_csv(run_sizes_path)
    runs['version'] = runs['version'].str.replace('version_', 'v')
    runs['version'] = runs['version'].str.replace('-UPDATED-BOWTIE', '')
    runs['version'] = runs['version'].str.replace('RC1', 'r1')
    runs.sort_values('version', inplace=True)
    plain_runs = runs.copy()
    plain_runs['size'] = plain_runs['outputs']
    plain_runs['type'] = 'unzipped'
    zipped_runs = runs.copy()
    zipped_runs['size'] = zipped_runs['zipped']
    zipped_runs['type'] = 'zipped'
    all_runs = pd.concat([plain_runs, zipped_runs])

    sns.boxplot(x='size', y='version', hue='type', data=all_runs)
    plt.xlabel('Output size (MB)')
    plt.title('MiSeq Disk Usage')

    plt.show()


if __name__ == '__main__':
    main()
