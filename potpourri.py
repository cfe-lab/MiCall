#!/usr/bin/env python

import argparse
from glob import glob
import os.path
import random
import time

import settings
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
    @return: [(run_folder, sample_name)] a list of tuples
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
        for sample_name in sample_sheet['Data'].iterkeys():
            sample_names.append((run, sample_name))
    
    sample_names = random.sample(sample_names, args.sample_count)
    sample_names.sort()
    return sample_names

def main():
    args = parseOptions()
    random.seed(args.seed)
    sample_names = choose_samples(args)
    for run, sample_name in sample_names:
        print run, sample_name
main()
