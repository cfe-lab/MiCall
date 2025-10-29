import csv
from collections import defaultdict
import logging
import os
from pathlib import Path

from micall.core.remap import remap
from micall.core.prelim_map import prelim_map
from csv import DictReader
from micall.core.trim_fastqs import trim, censor
from micall.utils.dd import DD
from micall.utils.work_dir import WorkDir
from micall.g2p.pssm_lib import Pssm

BOWTIE_THREADS = 11


class MicallDD(DD):
    def __init__(self, filename1, bad_cycles_filename):
        super(MicallDD, self).__init__()
        if True or 'filter' in filename1:
            self.filename1 = filename1
        else:
            self.filename1 = self.filter_fastqs(filename1)
        self.bad_cycles_filename = bad_cycles_filename
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
        remap_counts_filename = os.path.join(workdir, 'filter.remap_counts.csv')
        remap_conseq_filename = os.path.join(workdir, 'filter.remap_conseq.csv')
        unmapped1_filename = os.path.join(workdir, 'filter.unmapped1.fastq')
        unmapped2_filename = os.path.join(workdir, 'filter.unmapped2.fastq')

        print('Preliminary mapping to filter.')
        with WorkDir.using(Path(workdir)):
            prelim_map(Path(filename1),
                       Path(filename2),
                       Path(prelim_filename),
                       nthreads=BOWTIE_THREADS)
        print('Remapping to filter.')
        with WorkDir.using(Path(workdir)):
            remap(Path(filename1),
                  Path(filename2),
                  Path(prelim_filename),
                  Path(remap_filename),
                  Path(remap_counts_filename),
                  Path(remap_conseq_filename),
                  Path(unmapped1_filename),
                  Path(unmapped2_filename))
        with open(remap_filename, 'r') as remap_csv:
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
        with open(filename, 'r') as fin, open(filter_name, 'w') as fout:
            for read in zip(fin, fin, fin, fin):
                qname = read[0].split()[0][1:]
                if qname in qnames:
                    for line in read:
                        fout.write(line)

    def _test(self, read_indexes, debug_file_prefix=None):
        read_indexes = reversed(read_indexes)
        simple_filename1 = self.filename1 + '_simple.fastq'
        self.write_simple_fastq(simple_filename1, read_indexes)
        workdir = os.path.dirname(self.filename1)
        os.chdir(workdir)
        simple_filename2 = get_reverse_filename(simple_filename1)
        censored_filename1 = os.path.join(workdir, 'rerun.censored1.fastq')
        censored_filename2 = os.path.join(workdir, 'rerun.censored2.fastq')
        trimmed_filename1 = os.path.join(workdir, 'rerun.trimmed1.fastq')
        trimmed_filename2 = os.path.join(workdir, 'rerun.trimmed2.fastq')
        prelim_censored_filename = os.path.join(workdir, 'rerun_censored.prelim.csv')
        prelim_trimmed_filename = os.path.join(workdir, 'rerun_trimmed.prelim.csv')
        with open(self.bad_cycles_filename, 'r') as bad_cycles:
            bad_cycles = list(csv.DictReader(bad_cycles))
        with open(simple_filename1, 'r') as simple1, \
                open(censored_filename1, 'w') as censored1:
            censor(simple1, bad_cycles, censored1, use_gzip=False)
        with open(simple_filename2, 'r') as simple2, \
                open(censored_filename2, 'w') as censored2:
            censor(simple2, bad_cycles, censored2, use_gzip=False)
        with WorkDir.using(Path(workdir)):
            prelim_map(Path(censored_filename1),
                       Path(censored_filename2),
                       Path(prelim_censored_filename),
                       nthreads=BOWTIE_THREADS)
            trim((simple_filename1, simple_filename2),
                 self.bad_cycles_filename,
                 (trimmed_filename1, trimmed_filename2),
                 use_gzip=False)
            prelim_map(Path(trimmed_filename1),
                       Path(trimmed_filename2),
                       Path(prelim_trimmed_filename),
                       nthreads=BOWTIE_THREADS)
        with open(prelim_censored_filename, 'r') as prelim_censored_csv:
            censored_map_count = self.count_mapped(prelim_censored_csv)
        with open(prelim_trimmed_filename, 'r') as prelim_trimmed_csv:
            trimmed_map_count = self.count_mapped(prelim_trimmed_csv)

        return self.get_result(censored_map_count, trimmed_map_count)

    def count_mapped(self, prelim_csv):
        # return sum(1
        #            for row in csv.DictReader(prelim_csv)
        #            if row['rname'] != '*')
        n = 0
        for row in csv.DictReader(prelim_csv):
            if int(row['flag']) & 4 == 0:
                n += 1
        return n

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

    def get_result(self, censored_map_count, trimmed_map_count):
        diff = trimmed_map_count - censored_map_count
        print('Result: censored {}, trimmed {} ({}).'.format(censored_map_count,
                                                             trimmed_map_count,
                                                             diff))
        return DD.FAIL if diff >= 20 else DD.PASS

    def write_simple_fastq(self, filename1, read_indexes):
        selected_reads = map(self.reads.__getitem__, read_indexes)
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
    if 'censored1.fastq' in fastq1_filename:
        reverse_filename = fastq1_filename.replace('censored1.fastq', 'censored2.fastq')
    else:
        reverse_filename = fastq1_filename.replace('_R1_', '_R2_')
    assert reverse_filename != fastq1_filename
    return reverse_filename


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
    logger = logging.getLogger(__name__)
    try:
        logger.info('Starting.')
        fname = 'censored1.fastq'
        bad_cycles_filename = 'bad_cycles.csv'
        dd = MicallDD(fname, bad_cycles_filename)
        read_indexes = range(len(dd.reads))
        run_test = True
        if run_test:
            min_indexes = dd.ddmin(read_indexes)
        else:
            min_indexes = read_indexes
        dd._test(min_indexes, debug_file_prefix='micall_debug')
        # dd.write_simple_fastq(fname + '_min.fastq', min_indexes)
        logger.info('Done.')
    except Exception as ex:
        logger.error('Failed.', exc_info=ex)


if __name__ == '__main__':
    main()
