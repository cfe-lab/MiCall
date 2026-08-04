import sys
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
from collections import defaultdict
import logging
import os
import re
from csv import DictReader
from typing import Union

from micall.core.trim_fastqs import trim
from micall.utils.dd import DD
from micall.core.aln2counts import aln2counts

logger = logging.getLogger(__name__)

ALIGNED_CSV_HEADER = 'refname,qcut,rank,count,offset,seq'
SUBSEQ_ENV_VARNAME = 'MICALL_DD_SUBSEQ'

def parse_args(argv):
    parser = ArgumentParser(
        description='Find the simplest list of aligned reads that reproduces a chosen problem.',
        formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('filename',
                        help='Input file with the initial aligned reads')
    parser.add_argument('simple',
                        help='Output file with the simplified aligned reads')
    parser.add_argument('--test',
                        help='name of the test to run',
                        choices=MicallDD.test_names,
                        default=MicallDD.test_names[0])

    return parser.parse_args(argv)


class MicallDD(DD):
    test_names = ['subseq']

    def __init__(self,
                 filename,
                 simple,
                 test_name):
        super(MicallDD, self).__init__()
        self.filename = filename
        self.simple = simple
        self.get_result = getattr(self, 'check_' + test_name)
        self.reads = read_aligned(self.filename)

        expected_subsequence = os.environ.get(SUBSEQ_ENV_VARNAME, None)
        if expected_subsequence is None:
            raise RuntimeError(f"Expected ${SUBSEQ_ENV_VARNAME!r} environment variable to be set for the 'subseq' test")

        self.expected_subsequence_re = re.compile(expected_subsequence)

    def _test(self, read_indexes):
        read_count = len(read_indexes)
        self.write_simple_aligned(self.simple, read_indexes)
        workdir = os.path.dirname(self.simple)

        def writer(filename):
            return open(os.path.join(workdir, filename), 'w+')

        with open(self.simple, 'r') as aligned_csv, \
             writer('nuc.csv') as nuc_csv, \
             writer('amino.csv') as amino_csv, \
             writer('insertions.csv') as insertions_csv, \
             writer('conseq.csv') as conseq_csv, \
             writer('failed_align.csv') as failed_align_csv, \
             writer('nuc_detail.csv') as nuc_detail_csv, \
             writer('stitched.csv') as stitched_csv:

            # noinspection PyBroadException
            try:
                aln2counts(# Inputs #
                           aligned_csv,
                           nuc_csv,
                           amino_csv,
                           insertions_csv,
                           conseq_csv,
                           failed_align_csv,
                           # Outputs #
                           nuc_detail_csv=nuc_detail_csv,
                           )

                exception = None
            except Exception as ex:
                logger.warning(f'Read counting failed: {ex!r}.', exc_info=True)
                exception = ex

            stitched_csv.seek(0)
            result = self.get_result(stitched_csv, read_count, exception)
            if result == DD.FAIL:
                save_best(aligned_csv)
                save_best(nuc_csv)
                save_best(amino_csv)
                save_best(insertions_csv)
                save_best(conseq_csv)
                save_best(failed_align_csv)
                save_best(stitched_csv)

            return result

    def check_subseq(self, stitched_csv, read_count, exception):
        if exception is not None:
            return DD.UNRESOLVED

        simple_count = len(stitched_csv.readlines()) - 1

        logger.debug('Result: %d stitched sequences from %d selected reads.',
                     simple_count, read_count)

        stitched_csv.seek(0)
        success = self.expected_subsequence_re.search(stitched_csv.read())

        return DD.FAIL if success else DD.PASS

    def write_simple_aligned(self, filename, read_indexes):
        selected_reads = (self.reads[i] for i in read_indexes)
        with open(filename, 'w') as f:
            f.write(ALIGNED_CSV_HEADER)
            f.write('\n')
            for line in selected_reads:
                f.write(line)

    def coerce(self, c):
        if c is None:
            return 'None'
        blocks = []  # [[first, last]] indexes for all contiguous blocks
        for i in c:
            if (not blocks) or blocks[-1][-1] != i-1:
                blocks.append([i, i])
            else:
                blocks[-1][-1] = i
        return '[' + ', '.join(str(block[0]) if block[0] == block[1]
                               else '{}-{}'.format(*block)
                               for block in blocks) + ']'


def save_best(file: Union[str, '_io.TextIOWrapper']):
    """ Save the current best version of a file.
    """

    filename = file if type(file) is str else file.name
    base, ext = os.path.splitext(filename)
    best = base + '_best' + ext

    os.rename(filename, best)


def read_aligned(filename):
    """ Load all the reads from an aligned reads file into a dictionary.

    @param filename: the aligned.csv file to open
    @param reads: dict({index: line})
    """

    with open(filename, 'r') as f:
        header = next(f)

        # Sanity check that may detect instances where an incorrect file has been passed as input.
        if header.strip() != ALIGNED_CSV_HEADER.strip():
            raise ValueError(f'Aligned reads file {filename!r} does not start with a known header')

        return f.readlines()


def main(argv):
    logging.basicConfig(
        level=logging.DEBUG,
        format='%(asctime)s[%(levelname)s]%(module)s:%(lineno)d - %(message)s',
        stream=sys.stdout)
    args = parse_args(argv)
    try:
        logger.info('Starting.')
        dd = MicallDD(args.filename,
                      args.simple,
                      args.test)
        read_indexes = list(range(len(dd.reads)))
        min_indexes = dd.ddmin(read_indexes)
        dd.test(min_indexes)
        logger.info('Done.')
    except Exception as ex:
        logger.error('Failed.', exc_info=ex)


if __name__ == '__main__':
    main(sys.argv[1:])
