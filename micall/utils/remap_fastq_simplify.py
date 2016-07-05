from collections import defaultdict
from itertools import imap
import logging
import os

from micall.core.miseq_logging import init_logging_console_only
from micall.core.remap import remap
from micall.core.prelim_map import prelim_map
from csv import DictReader
from micall.core.sam2aln import sam2aln
from micall.core.aln2counts import aln2counts, AMINO_ALPHABET
from micall.utils.dd import DD


class MicallDD(DD):
    def __init__(self, filename1):
        super(MicallDD, self).__init__()
        self.filename1 = filename1
        reads = defaultdict(list)
        read_fastq(filename1, reads)
        read_count = len(reads)
        read_fastq(get_reverse_filename(filename1), reads)
        added_count = len(reads) - read_count
        if added_count > 0:
            raise RuntimeError('Found {} new reads.'.format(added_count))
        self.reads = reads.values()

    def _test(self, read_indexes, debug_file_prefix=None):
        BOWTIE_THREADS = 11
        simple_filename1 = self.filename1 + '_simple.fastq'
        self.write_simple_fastq(simple_filename1, read_indexes)
        workdir = os.path.dirname(self.filename1)
        simple_filename2 = get_reverse_filename(simple_filename1)
        prelim_filename = os.path.join(workdir, 'temp70.prelim.csv')
        remap_filename = os.path.join(workdir, 'temp70.remap.csv')
        remap_counts_filename = os.path.join(workdir, 'temp70.remap_counts.csv')
        aligned_filename = os.path.join(workdir, 'temp70.aligned.csv')
        nuc_filename = os.path.join(workdir, 'temp70.nuc.csv')
        amino_filename = os.path.join(workdir, 'temp70.amino.csv')
        failed_align_filename = os.path.join(workdir, 'temp70.failed_align.csv')
        with open(prelim_filename, 'w+') as prelim_csv, \
                open(remap_filename, 'w+') as remap_csv, \
                open(remap_counts_filename, 'w+') as remap_counts_csv, \
                open(aligned_filename, 'w+') as aligned_csv, \
                open(nuc_filename, 'w+') as nuc_csv, \
                open(amino_filename, 'w+') as amino_csv, \
                open(failed_align_filename, 'w+') as failed_align_csv, \
                open(os.devnull, 'w+') as real_devnull:
            devnull = DevNullWrapper(real_devnull)
            prelim_map(simple_filename1,
                       simple_filename2,
                       prelim_csv,
                       nthreads=BOWTIE_THREADS)
            prelim_csv.seek(0)
            remap(simple_filename1,
                  simple_filename2,
                  prelim_csv,
                  remap_csv,
                  remap_counts_csv,
                  devnull,
                  devnull,
                  devnull,
                  nthreads=BOWTIE_THREADS,
                  debug_file_prefix=debug_file_prefix)
            remap_csv.seek(0)
            sam2aln(remap_csv,
                    aligned_csv,
                    devnull,
                    failed_align_csv,
                    nthreads=BOWTIE_THREADS)
            aligned_csv.seek(0)
            aln2counts(aligned_csv,
                       nuc_csv,
                       amino_csv,
                       devnull,
                       devnull,
                       devnull,
                       devnull)

        with open(amino_filename, 'rU') as amino_csv:
            amino_reader = DictReader(amino_csv)
            region_counts = [row
                             for row in amino_reader
                             if row['region'] == 'HCV3-S52-E2']
        region_coverage = [sum(int(row[aa]) for aa in AMINO_ALPHABET if aa != '*')
                           for row in region_counts]
        if not region_coverage:
            return DD.PASS
        max_coverage = max(region_coverage)
        if max_coverage < 5 or region_coverage[1] < max_coverage * 0.7:
            # print('{} max {}'.format(region_coverage[1], max_coverage))
            return DD.PASS
        if region_coverage[0] == 0:
            return DD.FAIL
        # print('{}, {} max {}'.format(region_coverage[0], region_coverage[1], max_coverage))
        return DD.PASS

    def write_simple_fastq(self, filename1, read_indexes):
        selected_reads = imap(self.reads.__getitem__, read_indexes)
        filename2 = get_reverse_filename(filename1)
        with open(filename1, 'w') as f1, open(filename2, 'w') as f2:
            for lines in selected_reads:
                for line in lines[:4]:
                    f1.write(line)
                for line in lines[4:]:
                    f2.write(line)


class DevNullWrapper(object):
    def __init__(self, devnull):
        self.devnull = devnull

    def __getattr__(self, name):
        return getattr(self.devnull, name)

    def truncate(self):
        pass


def get_reverse_filename(fastq1_filename):
    return fastq1_filename.replace('censored1.fastq', 'censored2.fastq')


def read_fastq(filename, reads):
    """ Load all the reads from a FASTQ file into a list.

    @param filename: the FASTQ file to open
    @param reads: defaultdict({line1: [line2, line3, line4, line2, line3, line4]}
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
    logger = init_logging_console_only(logging.INFO)
    try:
        logger.info('Starting.')
        fname = ('censored1.fastq')
        dd = MicallDD(fname)
        read_indexes = range(len(dd.reads))
        min_indexes = dd.ddmin(read_indexes[50:85])
        dd._test(min_indexes, debug_file_prefix='micall_debug')
        # dd.write_simple_fastq(fname + '_min.fastq', min_indexes)
        logger.info('Done.')
    except:
        logger.error('Failed.', exc_info=True)

if __name__ == '__main__':
    main()
