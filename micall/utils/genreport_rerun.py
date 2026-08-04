import csv
import os
import shutil
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter, FileType
from csv import DictReader, DictWriter
from glob import glob
from itertools import groupby
from operator import itemgetter

from datetime import datetime

from micall.resistance.genreport import gen_report
from micall.monitor.find_groups import find_groups
from micall.resistance.resistance import report_resistance


def parse_args():
    parser = ArgumentParser(
        description='Rerun resistance interpretations.',
        formatter_class=ArgumentDefaultsHelpFormatter)

    parser.add_argument(
        '-r',
        '--results',
        help='source results folder in one run')
    parser.add_argument(
        '-w',
        '--working',
        help='working folder')
    parser.add_argument(
        '-d',
        '--raw_data',
        help='Main RAW_DATA folder to search for runs.')
    parser.add_argument(
        '-s',
        '--samples_csv',
        type=FileType(),
        help='CSV file with sample and run names.')
    parser.add_argument(
        '-m',
        '--min_run_name',
        help='Select all later runs, (e.g., "161201").')
    parser.add_argument(
        '--pipeline_version',
        default='0-dev',
        help='version number')
    args = parser.parse_args()
    source_names = ('results', 'samples_csv', 'min_run_name')
    source_count = sum(1 for name in source_names if getattr(args, name) is not None)
    if source_count != 1:
        parser.error('Must provide one of --results, --samples_csv, or --min_run_name.')
    if args.results is None and args.raw_data is None:
        parser.error('Must provide --raw_data if not using --results.')
    return args


def rewrite_file(filename):
    backup_csv = filename + '.original.csv'
    os.rename(filename, backup_csv)
    with open(backup_csv) as source, open(filename, 'w') as dest:
        reader = DictReader(source)
        writer = DictWriter(dest, reader.fieldnames)
        writer.writeheader()
        for row in reader:
            yield row  # Caller can modify it.
            if row:  # Skips a row that got cleared.
                writer.writerow(row)


def genreport_rerun(source, working):
    run_path = os.path.dirname(os.path.dirname(source))
    run_name = os.path.basename(run_path)
    publish_path = os.path.join(working,
                                'rerun_results',
                                run_name)
    os.makedirs(publish_path)
    print('##', run_name)
    working_paths = split_files(source, working)
    file_names = (os.path.basename(working_path) for working_path in working_paths)
    sorted_file_names = sorted(file_names)
    sample_sheet_path = os.path.join(source, '../../SampleSheet.csv')
    groups = list(find_groups(sorted_file_names,
                              sample_sheet_path,
                              included_projects={'HCV'}))
    for group in groups:
        working_path = os.path.join(working, group.names[0])
        if group.names[1] is None:
            midi_name = ''
            midi_path = working_path
        else:
            midi_name = group.names[1]
            midi_path = os.path.join(
                working,
                group.names[1])
            if not os.path.isdir(midi_path):
                midi_name = 'failed MidHCV'
                midi_path = working_path
        print(working_path, midi_name)
        with open(os.path.join(working_path, 'amino.csv')) as amino_csv, \
                open(os.path.join(midi_path, 'amino.csv')) as midi_amino_csv, \
                open(os.path.join(working_path, 'nuc.csv')) as nuc_csv, \
                open(os.path.join(working_path, 'resistance.csv'), 'w') as resistance_csv, \
                open(os.path.join(working_path, 'mutations.csv'), 'w') as mutations_csv, \
                open(os.path.join(working_path, 'nuc_mutations.csv'), 'w') as nuc_mutations_csv, \
                open(os.path.join(working_path, 'resistance_fail.csv'), 'w') as resistance_fail_csv:
            report_resistance(amino_csv,
                              midi_amino_csv,
                              nuc_csv,
                              resistance_csv,
                              mutations_csv,
                              nuc_mutations_csv,
                              resistance_fail_csv)
        sample_name = os.path.basename(working_path)
        with open(os.path.join(working_path, 'resistance.csv')) as resistance_csv, \
                open(os.path.join(working_path, 'mutations.csv')) as mutations_csv, \
                open(os.path.join(working_path, 'resistance_report.pdf'), 'wb') as resistance_report_csv:
            gen_report(resistance_csv,
                       mutations_csv,
                       resistance_report_csv,
                       sample_name=sample_name)
    for file_name in ('resistance.csv', 'mutations.csv', 'resistance_fail.csv'):
        with open(os.path.join(publish_path, file_name), 'w') as dest:
            dest_writer = csv.writer(dest)
            for i, group in enumerate(groups):
                working_path = os.path.join(working, group.names[0])
                sample_name = os.path.basename(working_path)
                with open(os.path.join(working_path, file_name), 'r') as source:
                    source_reader = csv.reader(source)
                    for j, row in enumerate(source_reader):
                        if j != 0:
                            row.insert(0, sample_name)
                        elif i == 0:
                            row.insert(0, 'sample')
                        else:
                            continue
                        dest_writer.writerow(row)


def split_files(source, working):
    working_paths = set()
    file_name = 'amino.csv'
    file_path = os.path.join(source, file_name)
    with open(file_path) as f:
        reader = DictReader(f)
        for sample, rows in groupby(reader, itemgetter('sample')):
            working_path = os.path.join(working, sample)
            working_paths.add(working_path)
            os.makedirs(working_path, exist_ok=True)
            target_path = os.path.join(working_path, file_name)
            with open(target_path, 'w') as target_csv:
                writer = DictWriter(target_csv,
                                    reader.fieldnames[1:])
                writer.writeheader()
                for row in rows:
                    del row['sample']
                    writer.writerow(row)
    return working_paths


def find_recent_results(min_run_name, source_folder, pipeline_version):
    pattern = os.path.join(source_folder, 'MiSeq', 'runs', '*')
    matches = glob(pattern)
    folder_names = [os.path.basename(match) for match in matches]
    recent_folders = sorted(os.path.join(source_folder,
                                         'MiSeq',
                                         'runs',
                                         folder,
                                         'Results',
                                         'version_' + pipeline_version)
                            for folder in folder_names
                            if folder >= min_run_name and folder[0].isnumeric())
    folders_with_data = [
        folder
        for folder in recent_folders
        if os.path.exists(os.path.join(folder, 'doneprocessing'))]
    return folders_with_data


def find_sample_results(samples_csv, source_folder, pipeline_version):
    run_names = {row['run'] for row in DictReader(samples_csv)}
    run_paths = []
    for run_name in run_names:
        run_path = os.path.join(source_folder, 'MiSeq', 'runs', run_name)
        if os.path.exists(run_path):
            run_paths.append(run_path)
        else:
            run_parts = run_name.split('.')
            if len(run_parts) < 2:
                run_parts.append('M01841')
            format_string = '%d-%b-%y' if len(run_parts[0]) == 9 else '%d-%b-%Y'
            run_date = datetime.strptime(run_parts[0], format_string)
            base_run_name = run_date.strftime('%y%m%d') + '_' + run_parts[1]
            pattern = os.path.join(source_folder,
                                   'MiSeq',
                                   'runs',
                                   base_run_name+'*',
                                   'needsprocessing')
            matches = glob(pattern)
            if len(matches) != 1:
                raise RuntimeError('Expected one match for {}, but found: {}'.format(
                    pattern,
                    matches))
            run_paths.append(os.path.dirname(matches[0]))
    run_paths.sort()
    results_folders = []
    for run_path in run_paths:
        done_processing_file = os.path.join(source_folder,
                                            'MiSeq',
                                            'runs',
                                            run_path,
                                            'Results',
                                            'version_' + pipeline_version,
                                            'doneprocessing')
        if os.path.exists(done_processing_file):
            results_folders.append(os.path.dirname(done_processing_file))
        else:
            print('!! Results not found:', done_processing_file)
    return results_folders


def clear_working_folder(working):
    for f in os.listdir(working):
        if f != 'rerun_results':
            file_path = os.path.join(working, f)
            if os.path.isdir(file_path):
                shutil.rmtree(file_path)
            else:
                os.remove(file_path)


def main():
    args = parse_args()
    if args.results is not None:
        results_folders = [args.results]
    elif args.min_run_name:
        results_folders = find_recent_results(args.min_run_name,
                                              args.raw_data,
                                              args.pipeline_version)
    else:
        results_folders = find_sample_results(args.samples_csv,
                                              args.raw_data,
                                              args.pipeline_version)
    rerun_results_path = os.path.join(args.working, 'rerun_results')
    shutil.rmtree(rerun_results_path, ignore_errors=True)
    for results_folder in results_folders:
        clear_working_folder(args.working)
        genreport_rerun(results_folder, args.working)
    clear_working_folder(args.working)

if __name__ == '__main__':
    main()
