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
from micall.g2p.pssm_lib import Pssm
from operator import itemgetter

BOWTIE_THREADS = 11


class MicallDD(DD):
    def __init__(self, filename1):
        super(MicallDD, self).__init__()
        if 'filter' in filename1:
            self.filename1 = filename1
        else:
            self.filename1 = self.filter_fastqs(filename1)
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
        prelim_filename = os.path.join(workdir, 'filter.prelim.csv')
        remap_filename = os.path.join(workdir, 'filter.remap.csv')
        with open(prelim_filename, 'w+') as prelim_csv, \
                open(remap_filename, 'w') as remap_csv, \
                open(os.devnull, 'w+') as real_devnull:
            print('Preliminary mapping to filter.')
            devnull = DevNullWrapper(real_devnull)
            prelim_map(filename1,
                       filename2,
                       prelim_csv,
                       nthreads=BOWTIE_THREADS)
            prelim_csv.seek(0)
            print('Remapping to filter.')
            remap(filename1,
                  filename2,
                  prelim_csv,
                  remap_csv,
                  devnull,
                  devnull,
                  devnull,
                  devnull,
                  nthreads=BOWTIE_THREADS)
        with open(remap_filename, 'rU') as remap_csv:
            print('Filtering.')
            reader = DictReader(remap_csv)
            mapped_qnames = {row['qname']
                             for row in reader
                             if row['rname'].startswith('HCV-1') and 2500 < int(row['pos']) < 3500}
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
        prelim_filename = os.path.join(workdir, 'rerun.prelim.csv')
        remap_filename = os.path.join(workdir, 'rerun.remap.csv')
        remap_counts_filename = os.path.join(workdir, 'rerun.remap_counts.csv')
        aligned_filename = os.path.join(workdir, 'rerun.aligned.csv')
        clipping_filename = os.path.join(workdir, 'rerun.clipping.csv')
        nuc_filename = os.path.join(workdir, 'rerun.nuc.csv')
        amino_filename = os.path.join(workdir, 'rerun.amino.csv')
        conseq_filename = os.path.join(workdir, 'rerun.conseq.csv')
        with open(prelim_filename, 'w+') as prelim_csv, \
                open(remap_filename, 'w+') as remap_csv, \
                open(remap_counts_filename, 'w+') as remap_counts_csv, \
                open(aligned_filename, 'w+') as aligned_csv, \
                open(clipping_filename, 'w+') as clipping_csv, \
                open(nuc_filename, 'w+') as nuc_csv, \
                open(amino_filename, 'w+') as amino_csv, \
                open(conseq_filename, 'w+') as conseq_csv, \
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
                    clipping_csv=clipping_csv,
                    nthreads=None)
            aligned_csv.seek(0)
            clipping_csv.seek(0)
            aln2counts(aligned_csv,
                       nuc_csv,
                       amino_csv,
                       devnull,
                       conseq_csv,
                       devnull,
                       devnull,
                       clipping_csv=clipping_csv)

        return self.get_result(amino_filename, len(read_indexes))

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

    def get_result(self, amino_filename, read_count):
        rows = []
        with open(amino_filename, 'rU') as amino_csv:
            amino_reader = DictReader(amino_csv)
            for row in amino_reader:
                row['refseq.aa.pos'] = int(row['refseq.aa.pos'])
                if 71 <= row['refseq.aa.pos'] <= 80 and row['region'] == 'HCV1A-H77-NS2':
                    row['coverage'] = sum(map(int,
                                              (row[amino]
                                               for amino in AMINO_ALPHABET)))
                    row['coverage'] += int(row['del'])
                    row['clip'] = int(row['clip'])
                    rows.append(row)
        clipped_rows = [row
                        for row in rows
                        if 71 <= row['refseq.aa.pos'] <= 75]
        unclipped_rows = [row
                          for row in rows
                          if 76 <= row['refseq.aa.pos'] <= 80]
        if not (clipped_rows and unclipped_rows):
            print len(clipped_rows), len(unclipped_rows)
            return DD.PASS
        min_clips = min(row['clip'] for row in clipped_rows)
        max_coverage = max(row['coverage'] for row in unclipped_rows)
        print map(itemgetter('coverage'), rows)
        print map(itemgetter('clip'), rows)
        return DD.FAIL if max_coverage * 2 < min_clips else DD.PASS

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
        read_indexes = range(len(dd.reads))
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
