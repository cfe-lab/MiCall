import argparse
from collections import defaultdict
import glob
import logging
import os
import subprocess

from micall.core.miseq_logging import init_logging_console_only
from micall.core.remap import remap
from micall.core.prelim_map import prelim_map


def write_simple_fastq(filename1, reads):
    filename2 = get_reverse_filename(filename1)
    with open(filename1, 'w') as f1, open(filename2, 'w') as f2:
        for lines in reads:
            for line in lines[:4]:
                f1.write(line)
            for line in lines[4:]:
                f2.write(line)


def test(reads, simple_filename):
    """ Map with 6.8 and 7.0, then compare unmapped read counts.

    @return: 'PASS' if the results match or it's a difference we're not
        interested in, 'FAIL' otherwise
    """
    workdir = os.path.dirname(simple_filename)
    for bamfile in glob.glob(os.path.join(workdir, 'temp.bam*')):
        os.remove(bamfile)

    write_simple_fastq(simple_filename, reads)
    unmapped_count68 = remap68(simple_filename)
    unmapped_count70 = remap70(simple_filename)
    if unmapped_count68 == unmapped_count70:
        return 'PASS'
    print '6.8: {}, 7.0: {}'.format(unmapped_count68, unmapped_count70)
    return 'FAIL'


def remap68(fastq1_filename):
    workdir = os.path.dirname(fastq1_filename)
    fastq2_filename = get_reverse_filename(fastq1_filename)
    prelim_filename = os.path.join(workdir, 'temp68.prelim.csv')
    remap_filename = os.path.join(workdir, 'temp.remap.csv')
    remap_counts_filename = os.path.join(workdir, 'temp.remap_counts.csv')
    remap_conseq_filename = os.path.join(workdir, 'temp.remap_conseq.csv')
    unmapped1_filename = os.path.join(workdir, 'temp.unmapped1.fastq')
    unmapped2_filename = os.path.join(workdir, 'temp.unmapped2.fastq')
    subprocess.check_call(['/mnt/data/don/git/MiCall6_8/prelim_map.py',
                           fastq1_filename,
                           fastq2_filename,
                           prelim_filename],
                          cwd='/mnt/data/don/git/MiCall6_8')
    subprocess.check_call(['/mnt/data/don/git/MiCall6_8/remap.py',
                           fastq1_filename,
                           fastq2_filename,
                           prelim_filename,
                           remap_filename,
                           remap_counts_filename,
                           remap_conseq_filename,
                           unmapped1_filename,
                           unmapped2_filename],
                          cwd='/mnt/data/don/git/MiCall6_8')
    return get_unmapped_counts(remap_counts_filename)


def get_unmapped_counts(remap_counts_filename):
    with open(remap_counts_filename, 'rU') as remap_counts_csv:
        unmapped_counts_line = list(remap_counts_csv)[-1]
        unmapped_counts = int(unmapped_counts_line.split(',')[-2])
    return unmapped_counts


def remap70(fastq1_filename):
    workdir = os.path.dirname(fastq1_filename)
    fastq2_filename = get_reverse_filename(fastq1_filename)
    prelim_filename = os.path.join(workdir, 'temp.prelim.csv')
    remap_filename = os.path.join(workdir, 'temp.remap.csv')
    remap_counts_filename = os.path.join(workdir, 'temp.remap_counts.csv')
    remap_conseq_filename = os.path.join(workdir, 'temp.remap_conseq.csv')
    unmapped1_filename = os.path.join(workdir, 'temp.unmapped1.fastq')
    unmapped2_filename = os.path.join(workdir, 'temp.unmapped2.fastq')
    with open(prelim_filename, 'w+') as prelim_csv, \
            open(remap_filename, 'w+') as remap_csv, \
            open(remap_counts_filename, 'w+') as remap_counts_csv, \
            open(remap_conseq_filename, 'w+') as remap_conseq_csv, \
            open(unmapped1_filename, 'w+') as unmapped1, \
            open(unmapped2_filename, 'w+') as unmapped2:
        prelim_map(fastq1_filename, fastq2_filename, prelim_csv)
        prelim_csv.seek(0)
        remap(fastq1_filename,
              fastq2_filename,
              prelim_csv,
              remap_csv,
              remap_counts_csv,
              remap_conseq_csv,
              unmapped1,
              unmapped2)
    return get_unmapped_counts(remap_counts_filename)


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
            n = min(n*2, len(reads))
            if n == len(reads):
                break

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
    unmapped_count68 = remap68(txtfilename)
    unmapped_count70 = remap70(txtfilename)

    if unmapped_count68 != unmapped_count70:
        print 'Mismatch with sample', txtfilename
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
    for txtfilename in glob.glob(os.path.join(args.workdir, args.pattern)):
        if txtfilename.find('simple') < 0:
            logger.info(os.path.basename(txtfilename))
            compare_remap(txtfilename, logger)
    logger.info('Done.')
