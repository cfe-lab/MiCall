import sys
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
from collections import defaultdict
import logging
import os
from csv import DictReader

from micall.core.trim_fastqs import trim
from micall.utils.dd import DD
from micall.core.aln2counts import aln2counts

logger = logging.getLogger(__name__)


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
        base, ext = os.path.splitext(simple)
        self.best = base + '_best' + ext
        self.get_result = getattr(self, 'check_' + test_name)
        self.reads = read_aligned(self.filename)

    def _test(self, read_indexes):
        read_count = len(read_indexes)
        self.write_simple_aligned(self.simple, read_indexes)
        workdir = os.path.dirname(self.simple)

        def writer(filename):
            return open(os.path.join(workdir, filename), 'w+')

        with open(self.simple, 'r') as aligned_csv, \
             writer('stitched.csv') as output:
            # noinspection PyBroadException
            try:
                aln2counts(aligned_csv,
                           # TODO: maybe redirect to os.devnull instead.
                           writer('nuc.csv'),
                           writer('amino.csv'),
                           writer('insertions.csv'),
                           writer('conseq.csv'),
                           writer('failed_align.csv'),
                           # --- #
                           conseq_stitched_csv=output,
                           nuc_detail_csv=writer("nuc_detail.csv"),
                           )

                exception = None
            except Exception as ex:
                logger.warning(f'Assembly failed: {ex!r}.', exc_info=True)
                print(f'Assembly failed: {ex!r}.')
                exception = ex

            output.seek(0)
            result = self.get_result(output, read_count, exception)
            if result == DD.FAIL:
                os.rename(self.simple, self.best)
            return result

    @staticmethod
    def check_subseq(output, read_count, exception):
        if exception is not None:
            return DD.UNRESOLVED

        simple_count = len(output.readlines()) - 1

        logger.debug('Result: %d simplified reads from %d selected reads.',
                     simple_count,
                     read_count)

        expected_substring = os.environ.get('MICALL_DD_SUBSEQ', None)
        if expected_substring is None:
            raise RuntimeError(f"Expected ${'MICALL_DD_SUBSEQ'!r} environment value to be set for the {'subseq'!r} test")
        output.seek(0)
        success = any((expected_substring in line) for line in output)

        return DD.FAIL if success else DD.PASS

    def write_simple_aligned(self, filename, read_indexes):
        selected_reads = (self.reads[i] for i in read_indexes)
        with open(filename, 'w') as f:
            f.write('refname,qcut,rank,count,offset,seq\n')
            for line in selected_reads:
                f.write(line)


def read_aligned(filename):
    """ Load all the reads from an aligned reads file into a dictionary.

    @param filename: the aligned.csv file to open
    @param reads: dict({index: line})
    """

    with open(filename, 'r') as f:
        header = next(f)

        # Sanity check that may detect instances where an incorrect file has been passed as input.
        if header.strip() != 'refname,qcut,rank,count,offset,seq':
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
