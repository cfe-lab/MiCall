#!/usr/bin/env python

import argparse
from datetime import datetime
from glob import glob
import os.path
import random
import time

from micall import settings
from sample_sheet_parser import sample_sheet_parser


def parseOptions():
    parser = argparse.ArgumentParser(
        description='Assemble a mixed set of sample files in a new run folder.')

    parser.add_argument('source_runs',
                        help='path to raw data run folders')
    parser.add_argument('--sample_count',
                        '-n',
                        help='number of samples to assemble',
                        type=int,
                        default=96)
    parser.add_argument('--seed',
                        '-s',
                        help='random number seed, defaults to timer',
                        type=int,
                        default=long(time.time() * 256))

    return parser.parse_args()


def choose_samples(args):
    """ Randomly choose samples from all runs under the source runs folder.

    @param args: the command line arguments
    @return: [(run_folder, sample_name, sample_line)] a list of tuples where
        sample_line is the original line from its sample sheet.
    """
    runs = glob(os.path.join(args.source_runs, '*', settings.NEEDS_PROCESSING))
    runs = map(os.path.dirname, runs)
    sample_names = []
    for run in runs:
        if os.path.exists(os.path.join(run, settings.ERROR_PROCESSING)):
            continue
        sample_sheet_file = os.path.join(run, 'SampleSheet.csv')
        with open(sample_sheet_file, 'rU') as f:
            sample_sheet = sample_sheet_parser(f)
            f.seek(0)
            for line in f:
                if line == '[Data]\n':
                    break
            sample_lines = list(f)

        for sample_name, sample_data in sample_sheet['Data'].iteritems():
            files_exist = (
                os.path.isfile(calculate_data_file_path(run, sample_name, 1)) and
                os.path.isfile(calculate_data_file_path(run, sample_name, 2)))
            if files_exist:
                sample_number = int(sample_name.split('_')[-1][1:])
                if sample_number < len(sample_lines):
                    sample_line = sample_lines[sample_number]
                    sample_names.append((run,
                                         sample_name,
                                         sample_line,
                                         sample_data['tags']))

    sample_names = random.sample(sample_names, args.sample_count)
    sample_names.sort()
    chosen_samples = []
    for i, (run, sample_name, sample_line, tags) in enumerate(sample_names):
        row = i / 12
        column = i % 12
        new_tags = "N%d-N%d" % (row+501, column+701)
        retagged_line = sample_line.replace(tags, new_tags)
        chosen_samples.append((run, sample_name, retagged_line))
    return chosen_samples


def calculate_data_file_path(run, sample_name, read_number):
    read_file = '{}_L001_R{}_001.fastq.gz'.format(sample_name, read_number)
    read_path = os.path.join(run,
                             'Data',
                             'Intensities',
                             'BaseCalls',
                             read_file)
    return read_path


def create_target_run():
    name_count = 1
    while True:
        run_name = datetime.now().strftime('%y%m%d_potpourri_test')
        if name_count > 1:
            run_name = '{}_{}'.format(run_name, name_count)
        target_run = os.path.join(settings.rawdata_mount,
                                  'MiSeq',
                                  'runs',
                                  run_name)
        if not os.path.exists(target_run):
            break
        name_count += 1

    os.makedirs(os.path.join(target_run, 'Data', 'Intensities', 'BaseCalls'))
    return target_run


def create_links(target_run, sample_names):
    for sample_number, (run, sample_name, _) in enumerate(sample_names, start=1):
        for read_number in 1, 2:
            read_path = calculate_data_file_path(run, sample_name, read_number)
            name_parts = sample_name.split('_')
            name_parts[-1] = 'S{}'.format(sample_number)
            new_sample_name = '_'.join(name_parts)
            target_path = calculate_data_file_path(target_run,
                                                   new_sample_name,
                                                   read_number)
            os.symlink(read_path, target_path)


def main():
    args = parseOptions()
    random.seed(args.seed)

    target_run = create_target_run()
    print 'Creating run {} with seed {}'.format(target_run, args.seed)

    sample_names = choose_samples(args)
    create_links(target_run, sample_names)
    header = """\
[Header]
IEMFileVersion,3
Investigator Name,TST
Project Name,{date:%d-%b-%Y}
Experiment Name,{date:%d-%b-%Y}
Date,{date:%m/%d/%Y}
Workflow,GenerateFASTQ
Assay,Nextera
Description,Nextera
Chemistry,Nextera
[Reads]
251
251
[Settings]
[Data]
Sample_ID,Sample_Name,Sample_Plate,Sample_Well,index,index2,Sample_Project,Description,GenomeFolder
""".format(date=datetime.today())

    with open(os.path.join(target_run, 'SampleSheet.csv',), 'w') as f:
        f.write(header)
        for _run, _sample_name, sample_line in sample_names:
            f.write(sample_line)

    with open(os.path.join(target_run, settings.NEEDS_PROCESSING), 'w') as f:
        pass  # empty file
    print "Done. Don't forget to import the run in QAI."
main()
