from collections import defaultdict
from itertools import imap
import logging
import os

from micall.core.miseq_logging import init_logging_console_only
from micall.core.remap import remap
from micall.core.prelim_map import prelim_map
from csv import DictReader
from micall.core.sam2aln import sam2aln
from micall.core.aln2counts import aln2counts
from micall.utils.dd import DD
from micall.g2p.sam_g2p import sam_g2p, DEFAULT_MIN_COUNT
from micall.g2p.pssm_lib import Pssm

BOWTIE_THREADS = 11


class MicallDD(DD):
    def __init__(self, filename1):
        super(MicallDD, self).__init__()
        self.filename1 = filename1  # self.filter_fastqs(filename1)
        self.pssm = Pssm()
        reads = defaultdict(list)
        read_fastq(self.filename1, reads)
        read_count = len(reads)
        read_fastq(get_reverse_filename(self.filename1), reads)
        added_count = len(reads) - read_count
        if added_count > 0:
            raise RuntimeError('Found {} new reads.'.format(added_count))
        self.reads = reads.values()

    def filter_fastqs(self, filename1):
        filter_name1 = filename1 + '_filter.fastq'
        if os.path.exists(filter_name1):
            logging.info('Already filtered.')
            return filter_name1
        filename2 = get_reverse_filename(filename1)
        filter_name2 = filename2 + '_filter.fastq'
        workdir = os.path.dirname(filename1)
        prelim_filename = os.path.join(workdir, 'temp70.prelim.csv')
        with open(prelim_filename, 'w+') as prelim_csv:
            prelim_map(filename1,
                       filename2,
                       prelim_csv,
                       nthreads=BOWTIE_THREADS)
            prelim_csv.seek(0)
            reader = DictReader(prelim_csv)
            mapped_qnames = {row['qname']
                             for row in reader
                             if row['rname'].startswith('HCV') and 2400 < int(row['pos']) < 3700}
        self.filter_reads(filename1, filter_name1, mapped_qnames)
        self.filter_reads(filename2, filter_name2, mapped_qnames)
        logging.info('Finished filtering with %d reads.', len(mapped_qnames))
        return filter_name1

    def filter_reads(self, filename, filter_name, qnames):
        with open(filename, 'rU') as fin, open(filter_name, 'w') as fout:
            for read in zip(fin, fin, fin, fin):
                qname = read[0].split()[0][1:]
                if qname in qnames:
                    for line in read:
                        fout.write(line)

    def _test(self, read_indexes, debug_file_prefix=None):
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
        g2p_filename = os.path.join(workdir, 'temp70.g2p.csv')
        with open(prelim_filename, 'w+') as prelim_csv, \
                open(remap_filename, 'w+') as remap_csv, \
                open(remap_counts_filename, 'w+') as remap_counts_csv, \
                open(aligned_filename, 'w+') as aligned_csv, \
                open(nuc_filename, 'w+') as nuc_csv, \
                open(amino_filename, 'w+') as amino_csv, \
                open(g2p_filename, 'w+') as g2p_csv, \
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
                    devnull,
                    nthreads=BOWTIE_THREADS)
            aligned_csv.seek(0)
            aln2counts(aligned_csv,
                       nuc_csv,
                       amino_csv,
                       devnull,
                       devnull,
                       devnull,
                       devnull)
            remap_csv.seek(0)
            nuc_csv.seek(0)
            sam_g2p(pssm=self.pssm,
                    remap_csv=remap_csv,
                    nuc_csv=nuc_csv,
                    g2p_csv=g2p_csv,
                    g2p_summary_csv=devnull,
                    min_count=DEFAULT_MIN_COUNT)

        return self.get_result(g2p_filename, len(read_indexes))

    def disabled_resolve(self, csub, c, direction):
        sub_size = len(csub)
        if direction == DD.REMOVE:
            # result = csub[:sub_size/2]
            result = None
        else:
            # ADD
            add_count = (len(c) - sub_size + 1) / 2
            result = []
            for i in c:
                if i in csub:
                    result.append(i)
                elif add_count > 0:
                    result.append(i)
                    add_count -= 1
        return result

    def get_result(self, g2p_filename, read_count):
        with open(g2p_filename, 'rU') as g2p_csv:
            g2p_reader = DictReader(g2p_csv)
            rows = list(g2p_reader)
        good_counts = [int(row['count']) for row in rows if not row['error']]
        good_sum = sum(good_counts)
        good_max = good_counts and max(good_counts) or 0
        bad_count = sum(int(row['count']) for row in rows if row['error'] == 'zerolength')
        if good_max < 10 and (read_count - good_sum - bad_count != 0):
            # Other errors can be distracting.
            logging.info('PASS: %d good and %d bad out of %d.',
                         good_sum,
                         bad_count,
                         read_count)
            return DD.PASS
        return DD.FAIL if 0.4 * good_sum < bad_count <= good_sum else DD.PASS

    def write_simple_fastq(self, filename1, read_indexes):
        selected_reads = imap(self.reads.__getitem__, read_indexes)
        filename2 = get_reverse_filename(filename1)
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
    logger = init_logging_console_only(logging.INFO)
    try:
        logger.info('Starting.')
        fname = ('censored1.fastq')
        dd = MicallDD(fname)
        read_indexes = range(len(dd.reads))[800:1600]
        run_test = True
        if run_test:
            min_indexes = dd.ddmin(read_indexes)
        else:
            min_indexes = read_indexes
        dd._test(min_indexes, debug_file_prefix='micall_debug')
        # dd.write_simple_fastq(fname + '_min.fastq', min_indexes)
        logger.info('Done.')
    except:
        logger.error('Failed.', exc_info=True)

if __name__ == '__main__':
    main()
