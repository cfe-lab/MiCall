import argparse
from collections import defaultdict, Counter
import glob
import logging
import os
import subprocess

from micall.core.miseq_logging import init_logging_console_only
from micall.core.remap import remap
from micall.core.prelim_map import prelim_map
from csv import DictReader
from micall.core.sam2aln import sam2aln
from micall.core.aln2counts import aln2counts


def write_simple_fastq(filename1, reads):
    filename2 = get_reverse_filename(filename1)
    with open(filename1, 'w') as f1, open(filename2, 'w') as f2:
        for lines in reads:
            for line in lines[:4]:
                f1.write(line)
            for line in lines[4:]:
                f2.write(line)


def test(reads, simple_filename):
    """ Map with 6.8 and 7.0, then compare average coverage in NS3.

    @return: 'PASS' if 7.0 coverage is as good as or better than 6.8 coverage,
        'FAIL' otherwise
    """
    workdir = os.path.dirname(simple_filename)
    for bamfile in glob.glob(os.path.join(workdir, 'temp.bam*')):
        os.remove(bamfile)

    write_simple_fastq(simple_filename, reads)
    return test_file(simple_filename)


def test_file(simple_filename):
    ns3_coverage68 = remap68(simple_filename, do_counts=True)
    ns3_coverage70 = remap70(simple_filename, do_counts=True)
    if ns3_coverage70 >= ns3_coverage68 * 0.75:
        print '6.8: {}, 7.0: {}'.format(ns3_coverage68, ns3_coverage70)
        return 'PASS'
    print '6.8: {}, 7.0: {}'.format(ns3_coverage68, ns3_coverage70)
    return 'FAIL'


def remap68(fastq1_filename, do_counts=False):
    workdir = os.path.dirname(fastq1_filename)
    rundir = '/mnt/data/don/git/MiCall6_8'
    fastq2_filename = get_reverse_filename(fastq1_filename)
    prelim_filename = os.path.join(workdir, 'temp68.prelim.csv')
    remap_filename = os.path.join(workdir, 'temp68.remap.csv')
    remap_counts_filename = os.path.join(workdir, 'temp68.remap_counts.csv')
    aligned_filename = os.path.join(workdir, 'temp68.aligned.csv')
    nuc_filename = os.path.join(workdir, 'temp68.nuc.csv')
    subprocess.check_call([os.path.join(rundir, 'prelim_map.py'),
                           fastq1_filename,
                           fastq2_filename,
                           prelim_filename],
                          cwd=rundir)
    subprocess.check_call([os.path.join(rundir, 'remap.py'),
                           fastq1_filename,
                           fastq2_filename,
                           prelim_filename,
                           remap_filename,
                           remap_counts_filename,
                           os.devnull,
                           os.devnull,
                           os.devnull],
                          cwd=rundir)
    if not do_counts:
        return get_max_mapped_counts(remap_counts_filename)
    subprocess.check_call([os.path.join(rundir, 'sam2aln.py'),
                           remap_filename,
                           aligned_filename,
                           os.devnull,
                           os.devnull],
                          cwd=rundir)
    subprocess.check_call([os.path.join(rundir, 'aln2counts.py'),
                           aligned_filename,
                           nuc_filename,
                           os.devnull,
                           os.devnull,
                           os.devnull,
                           os.devnull,
                           os.devnull],
                          cwd=rundir)
    return get_gene_coverage(nuc_filename, 'NS3')


def get_gene_coverage(nuc_filename, gene_name):
    gene_suffix = '-' + gene_name
    counts = Counter()  # {pos: count}
    with open(nuc_filename, 'rU') as nuc_csv:
        reader = DictReader(nuc_csv)
        for row in reader:
            if row['region'].endswith(gene_suffix):
                nuc_total = sum(map(int, (row[nuc] for nuc in 'ACGT')))
                nuc_pos = int(row['refseq.nuc.pos'])
                counts[nuc_pos] += nuc_total
    if not counts:
        return 0
    return sum(counts.itervalues()) / float(len(counts))


class DevNullWrapper(object):
    def __init__(self, devnull):
        self.devnull = devnull

    def __getattr__(self, name):
        return getattr(self.devnull, name)

    def truncate(self):
        pass


def get_max_mapped_counts(remap_counts_filename):
    prelim_count = remap_count = 0
    with open(remap_counts_filename, 'rU') as remap_counts_csv:
        reader = DictReader(remap_counts_csv)
        for row in reader:
            row_count = int(row['count'])
            row_type = row['type']
            if row_type.startswith('prelim'):
                prelim_count = max(row_count, prelim_count)
            elif row_type.startswith('remap'):
                remap_count = max(row_count, remap_count)
    if remap_count == prelim_count:
        return 0
    return remap_count


def remap70(fastq1_filename, do_counts=False):
    workdir = os.path.dirname(fastq1_filename)
    fastq2_filename = get_reverse_filename(fastq1_filename)
    prelim_filename = os.path.join(workdir, 'temp70.prelim.csv')
    remap_filename = os.path.join(workdir, 'temp70.remap.csv')
    remap_counts_filename = os.path.join(workdir, 'temp70.remap_counts.csv')
    aligned_filename = os.path.join(workdir, 'temp70.aligned.csv')
    nuc_filename = os.path.join(workdir, 'temp70.nuc.csv')
    failed_align_filename = os.path.join(workdir, 'temp70.failed_align.csv')
    with open(prelim_filename, 'w+') as prelim_csv, \
            open(remap_filename, 'w+') as remap_csv, \
            open(remap_counts_filename, 'w+') as remap_counts_csv, \
            open(aligned_filename, 'w+') as aligned_csv, \
            open(nuc_filename, 'w+') as nuc_csv, \
            open(failed_align_filename, 'w+') as failed_align_csv, \
            open(os.devnull, 'w+') as real_devnull:
        devnull = DevNullWrapper(real_devnull)
        prelim_map(fastq1_filename, fastq2_filename, prelim_csv)
        prelim_csv.seek(0)
        remap(fastq1_filename,
              fastq2_filename,
              prelim_csv,
              remap_csv,
              remap_counts_csv,
              devnull,
              devnull,
              devnull)
        if not do_counts:
            remap_counts_csv.close()
            return get_max_mapped_counts(remap_counts_filename)
        remap_csv.seek(0)
        sam2aln(remap_csv, aligned_csv, devnull, failed_align_csv)
        aligned_csv.seek(0)
        aln2counts(aligned_csv,
                   nuc_csv,
                   devnull,
                   devnull,
                   devnull,
                   devnull,
                   devnull)
    return get_gene_coverage(nuc_filename, 'NS3')


def ddmin(reads, simple_filename):
    # assert test(sam_lines) == "FAIL"

    n = 2     # Initial granularity
    while len(reads) >= 2:
        start = 0
        subset_length = len(reads) / n
        some_complement_is_failing = False

        while start < len(reads):
            complement = reads[:start] + reads[start + subset_length:]

            result = test(complement, simple_filename)
            complement_description = '0-{} + {}-{}'.format(start,
                                                           start+subset_length,
                                                           len(reads))
            print len(complement), result, complement_description, simple_filename

            if result == "FAIL":
                reads = complement
                n = max(n - 1, 2)
                some_complement_is_failing = True
                break

            start += subset_length

        if not some_complement_is_failing:
            if n == len(reads):
                break
            n = min(n*2, len(reads))

    return reads


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


def compare_remap(txtfilename, logger):
    result = test_file(txtfilename)

    if result == 'FAIL':
        print 'Simplifying sample {}'.format(txtfilename)
        reads = defaultdict(list)
        read_fastq(txtfilename, reads)
        read_count = len(reads)
        read_fastq(get_reverse_filename(txtfilename), reads)
        added_count = len(reads) - read_count
        if added_count > 0:
            raise RuntimeError('Found {} new reads.'.format(added_count))
        reads = reads.values()
        simple_filename = txtfilename.replace('censored1.fastq',
                                              'simple_censored1.fastq')
        simple_fastq_lines = ddmin(reads, simple_filename)
        write_simple_fastq(simple_filename, simple_fastq_lines)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Find the simplest test failure by trimming FASTQ files.')

    parser.add_argument('workdir', help='path to folder holding FASTQ files')
    parser.add_argument('--pattern',
                        default='*censored1.fastq',
                        help='File name pattern to match FASTQ files')

    args = parser.parse_args()

    logger = init_logging_console_only(logging.INFO)
    filenames = glob.glob(os.path.join(args.workdir, args.pattern))
    filenames.sort()
    for txtfilename in filenames:
        if txtfilename.find('simple') < 0:
            logger.info(os.path.basename(txtfilename))
            compare_remap(txtfilename, logger)
    logger.info('Done.')
