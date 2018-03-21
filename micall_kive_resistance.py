import logging
import shutil
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter

import os

from micall.drivers.run_info import RunInfo
from micall.drivers.sample import Sample
from micall.drivers.sample_group import SampleGroup

logger = logging.getLogger('foo')


def parse_args():
    parser = ArgumentParser(description='Map FASTQ files to references.',
                            formatter_class=ArgumentDefaultsHelpFormatter)
    # inputs
    parser.add_argument('main_amino_csv',
                        help='CSV containing amino frequencies from main sample')
    parser.add_argument('midi_amino_csv',
                        help='CSV containing amino frequencies from MIDI sample')

    # outputs
    parser.add_argument('resistance_csv',
                        help='CSV containing resistance calls.')
    parser.add_argument('mutations_csv',
                        help='CSV containing resistance mutations.')
    parser.add_argument('resistance_fail_csv',
                        help='CSV containing failure reasons.')
    parser.add_argument('resistance_pdf',
                        help='resistance report')

    return parser.parse_args()


def load_sample(args):
    """ Load the data from Kive's command-line arguments. """
    scratch_path = os.path.join(os.path.dirname(args.main_amino_csv), 'scratch')
    shutil.rmtree(scratch_path, ignore_errors=True)

    sample1 = Sample(amino_csv=args.main_amino_csv,
                     resistance_csv=args.resistance_csv,
                     mutations_csv=args.mutations_csv,
                     resistance_fail_csv=args.resistance_fail_csv,
                     resistance_pdf=args.resistance_pdf,
                     scratch_path=scratch_path)
    sample2 = Sample(amino_csv=args.midi_amino_csv)
    return SampleGroup(sample1, sample2)


def main():
    logging.basicConfig(level=logging.WARN)
    args = parse_args()
    sample_group = load_sample(args)

    sample_group.process_resistance(RunInfo([sample_group]))


main()
