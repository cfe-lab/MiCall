import sys
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
from collections import defaultdict
import logging
import os
from csv import DictReader

from micall.core.trim_fastqs import trim
from micall.utils.dd import DD
from micall.core.denovo import denovo
from micall.utils.work_dir import WorkDir
from pathlib import Path

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
    parser.add_argument('--test',
                        help='name of the test to run',
                        choices=MicallDD.test_names,
                        default=MicallDD.test_names[0])
    parser.add_argument('--project_code',
                        help='Project code for trimming primers')

    return parser.parse_args()


class MicallDD(DD):
    test_names = ('one_contig',
                  'multiple_genotypes',
                  'type_error',
                  'one_long_contig',
                  'one_medium_contig',
                  'two_long_contigs')

    def __init__(self,
                 filename1,
                 filename2,
                 bad_cycles_filename,
                 simple1,
                 simple2,
                 test_name,
                 project_code=None):
        super(MicallDD, self).__init__()
        self.filename1 = filename1
        self.filename2 = filename2
        self.bad_cycles_filename = bad_cycles_filename
        self.simple1 = simple1
        self.simple2 = simple2
        self.project_code = project_code
        base, ext = os.path.splitext(simple1)
        self.best1 = base + '_best' + ext
        base, ext = os.path.splitext(simple2)
        self.best2 = base + '_best' + ext
        self.get_result = getattr(self, 'check_' + test_name)
        reads = defaultdict(list)
        read_fastq(self.filename1, reads)
        read_count = len(reads)
        read_fastq(self.filename2, reads)
        added_count = len(reads) - read_count
        if added_count > 0:
            raise RuntimeError('Found {} new reads.'.format(added_count))
        self.reads = list(reads.values())

    def _test(self, read_indexes):
        read_count = len(read_indexes)
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
                 use_gzip=False,
                 project_code=self.project_code)
            exception = None
            # noinspection PyBroadException
            try:
                # Use WorkDir for dynamic scoping of work_dir
                with WorkDir.using(Path(workdir)):
                    denovo(Path(trimmed_filename1),
                           Path(trimmed_filename2),
                           Path(contigs_csv.name))
            except Exception as ex:
                logger.warning('Assembly failed.', exc_info=True)
                exception = ex
            contigs_csv.seek(0)

            result = self.get_result(contigs_csv, read_count, exception)
            if result == DD.FAIL:
                os.rename(self.simple1, self.best1)
                os.rename(self.simple2, self.best2)
            return result

    @staticmethod
    def check_one_contig(contigs_csv, read_count, exception):
        if exception is not None:
            return DD.UNRESOLVED
        contig_count = len(contigs_csv.readlines()) - 1
        logger.debug('Result: %d contigs from %d reads.',
                     contig_count,
                     read_count)
        return DD.FAIL if contig_count == 1 else DD.PASS

    @staticmethod
    def check_multiple_genotypes(contigs_csv, read_count, exception):
        if exception is not None:
            return DD.UNRESOLVED
        reader = DictReader(contigs_csv)
        genotypes = sorted({row['genotype'] for row in reader})
        genotype_count = len(genotypes)
        logger.debug('Result: %d genotypes from %d reads: %s.',
                     genotype_count,
                     read_count,
                     genotypes)
        return DD.FAIL if genotype_count > 2 else DD.PASS

    @staticmethod
    def check_type_error(_contigs_csv, read_count, exception):
        logger.debug('Result: %s exception from %d reads.',
                     exception,
                     read_count)
        return DD.FAIL if isinstance(exception, TypeError) else DD.PASS

    @staticmethod
    def check_one_long_contig(contigs_csv, read_count, exception):
        if exception is not None:
            return DD.UNRESOLVED
        contig_sizes = MicallDD.get_contig_sizes(contigs_csv, read_count)
        return (DD.FAIL
                if len(contig_sizes) == 1 and min(contig_sizes.values()) >= 1000
                else DD.PASS)

    @staticmethod
    def check_one_medium_contig(contigs_csv, read_count, exception):
        if exception is not None:
            return DD.UNRESOLVED
        contig_sizes = MicallDD.get_contig_sizes(contigs_csv, read_count)
        return (DD.FAIL
                if len(contig_sizes) == 1 and min(contig_sizes.values()) >= 585
                else DD.PASS)

    @staticmethod
    def check_two_long_contigs(contigs_csv, read_count, exception):
        if exception is not None:
            return DD.UNRESOLVED
        contig_sizes = MicallDD.get_contig_sizes(contigs_csv, read_count)
        return (DD.FAIL
                if len(contig_sizes) == 2 and min(contig_sizes.values()) >= 3000
                else DD.PASS)

    @staticmethod
    def get_contig_sizes(contigs_csv, read_count):
        reader = DictReader(contigs_csv)
        contig_sizes = {f"{row['ref']}[{i}]": len(row['contig'])
                        for i, row in enumerate(reader)}
        genotype_summary = ', '.join(
            f'{genotype} ({size})'
            for genotype, size in sorted(contig_sizes.items()))
        logger.debug('Result: contigs %s from %d reads.',
                     genotype_summary or 'not found',
                     read_count)
        return contig_sizes

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
    with open(filename, 'r') as f:
        for line1, line2, line3, line4 in zip(f, f, f, f):
            qname = line1.split()[0]
            lines = reads[qname]
            lines.append(line1)
            lines.append(line2)
            lines.append(line3)
            lines.append(line4)


def main():
    logging.basicConfig(
        level=logging.DEBUG,
        format='%(asctime)s[%(levelname)s]%(module)s:%(lineno)d - %(message)s',
        stream=sys.stdout)
    args = parse_args()
    try:
        logger.info('Starting.')
        dd = MicallDD(args.fastq1,
                      args.fastq2,
                      args.bad_cycles_csv,
                      args.simple1,
                      args.simple2,
                      args.test,
                      args.project_code)
        read_indexes = list(range(len(dd.reads)))
        min_indexes = dd.ddmin(read_indexes)
        dd.test(min_indexes)
        logger.info('Done.')
    except Exception as ex:
        logger.error('Failed.', exc_info=ex)


if __name__ == '__main__':
    main()
