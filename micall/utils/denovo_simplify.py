from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
from collections import defaultdict
import logging
import os

from micall.core.trim_fastqs import trim
from micall.utils.dd import DD
from micall.utils.denovo import main as denovo_main

logger = logging.getLogger(__name__)


def parse_args():
    parser = ArgumentParser(
        description='Find the simplest FASTQ that can still be assembled.',
        formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('fastq1',
                        help='FASTQ read 1 file')
    parser.add_argument('fastq2',
                        help='FASTQ read 2 file')
    parser.add_argument('bad_cycles_csv',
                        help='CSV file with bad cycles')
    parser.add_argument('simple1',
                        help='FASTQ read 1 simplified file to write')
    parser.add_argument('simple2',
                        help='FASTQ read 2 simplified file to write')

    return parser.parse_args()


class MicallDD(DD):
    def __init__(self,
                 filename1,
                 filename2,
                 bad_cycles_filename,
                 simple1,
                 simple2):
        super(MicallDD, self).__init__()
        self.filename1 = filename1
        self.filename2 = filename2
        self.bad_cycles_filename = bad_cycles_filename
        self.simple1 = simple1
        self.simple2 = simple2
        reads = defaultdict(list)
        read_fastq(self.filename1, reads)
        read_count = len(reads)
        read_fastq(self.filename2, reads)
        added_count = len(reads) - read_count
        if added_count > 0:
            raise RuntimeError('Found {} new reads.'.format(added_count))
        self.reads = list(reads.values())

    def _test(self, read_indexes):
        read_indexes = reversed(read_indexes)
        self.write_simple_fastq(self.simple1, self.simple2, read_indexes)
        workdir = os.path.dirname(self.simple1)
        os.chdir(workdir)
        trimmed_filename1 = os.path.join(workdir, 'rerun.trimmed1.fastq')
        trimmed_filename2 = os.path.join(workdir, 'rerun.trimmed2.fastq')
        contigs_filename = os.path.join(workdir, 'rerun.contigs.csv')
        with open(contigs_filename, 'w+') as contigs_csv:
            trim((self.simple1, self.simple2),
                 self.bad_cycles_filename,
                 (trimmed_filename1, trimmed_filename2),
                 use_gzip=False)
            # noinspection PyBroadException
            try:
                denovo_main(trimmed_filename1, trimmed_filename2, contigs_csv)
            except Exception:
                logger.warning('Assembly failed.', exc_info=True)
                return DD.UNRESOLVED
            contigs_csv.seek(0)
            contig_count = len(contigs_csv.readlines())

        return self.get_result(contig_count, expected_count=0)

    @staticmethod
    def get_result(contig_count, expected_count):
        diff = expected_count - contig_count
        print('Result: {} contigs, expected {} ({}).'.format(contig_count,
                                                             expected_count,
                                                             diff))
        return DD.FAIL if diff else DD.PASS

    def write_simple_fastq(self, filename1, filename2, read_indexes):
        selected_reads = (self.reads[i] for i in read_indexes)
        with open(filename1, 'w') as f1, open(filename2, 'w') as f2:
            for lines in selected_reads:
                for line in lines[:4]:
                    f1.write(line)
                for line in lines[4:]:
                    f2.write(line)

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


def read_fastq(filename, reads):
    """ Load all the reads from a FASTQ file into a list.

    @param filename: the FASTQ file to open
    @param reads: defaultdict({qname: [line1, line2, line3, line4, line1, line2, line3, line4]}
    """
    with open(filename, 'rU') as f:
        for line1, line2, line3, line4 in zip(f, f, f, f):
            qname = line1.split()[0]
            lines = reads[qname]
            lines.append(line1)
            lines.append(line2)
            lines.append(line3)
            lines.append(line4)


def main():
    args = parse_args()
    try:
        logger.info('Starting.')
        dd = MicallDD(args.fastq1,
                      args.fastq2,
                      args.bad_cycles_csv,
                      args.simple1,
                      args.simple2)
        read_indexes = list(range(len(dd.reads)))
        min_indexes = dd.ddmin(read_indexes)
        dd._test(min_indexes)
        logger.info('Done.')
    except Exception as ex:
        logger.error('Failed.', exc_info=ex)


if __name__ == '__main__':
    main()
