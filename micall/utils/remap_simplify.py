import argparse
import glob
import logging
import os
import tempfile

from micall import settings
from micall.core.miseq_logging import init_logging_console_only
from micall.core.remap import build_conseqs_with_python, \
    build_conseqs_with_samtools, line_counter
from micall.utils import externals


def write_simple_sam(samfile, sam_lines, included_positions):
    for line in sam_lines:
        if included_positions is not None and not line.startswith('@'):
            fields = line.rstrip('\n').split('\t')
            seq = list(fields[9])
            qual = list(fields[10])
            for i in range(len(seq)):
                if i not in included_positions:
                    seq[i] = 'A'
                    qual[i] = '#'
                fields[9] = ''.join(seq)
                fields[10] = ''.join(qual)
            
            line = '\t'.join(fields) + '\n'
        samfile.write(line)

def test(sam_lines, temp_prefix, samtools, included_positions=None):
    """ Build a consensus sequence using samtools and Python, then compare.
    
    @return: 'PASS' if the results match or it's a difference we're not
        interested in, 'FAIL' otherwise
    """
    workdir = os.path.dirname(temp_prefix)
    for bamfile in glob.glob(os.path.join(workdir, 'temp.bam*')):
        os.remove(bamfile)
    
    samfd, txtfilename = tempfile.mkstemp(suffix=".sam",
                                            prefix=temp_prefix,
                                            text=True)
    try:
        samfile = os.fdopen(samfd, 'w')
        try:
            write_simple_sam(samfile, sam_lines, included_positions)
        finally:
            samfile.close()
        try:
            samtools_conseqs = build_conseqs_with_samtools(txtfilename,
                                                           samtools,
                                                           raw_count=len(sam_lines))
        except:
            return 'PASS' # We want simplest valid input that gives diff. result
        
        try:
            python2_conseqs = build_conseqs_with_python(txtfilename)
        except:
            return 'FAIL' #It was valid for samtools, want it valid for Python
        
        if samtools_conseqs.keys() != python2_conseqs.keys():
            # samtools only uses reads that map to regions listed in the header,
            # so we don't care if it ignored a whole region.
            return 'PASS'
        result = 'PASS' if samtools_conseqs == python2_conseqs else 'FAIL'
        print len(sam_lines), result, txtfilename
        return result
    finally:
        os.remove(txtfilename)

def ddmin(sam_lines, temp_prefix, samtools):
    #assert test(sam_lines) == "FAIL"

    n = 2     # Initial granularity
    while len(sam_lines) >= 2:
        start = 0
        subset_length = len(sam_lines) / n
        some_complement_is_failing = False

        while start < len(sam_lines):
            complement = sam_lines[:start] + sam_lines[start + subset_length:]

            if test(complement, temp_prefix, samtools) == "FAIL":
                sam_lines = complement
                n = max(n - 1, 2)
                some_complement_is_failing = True
                break
                
            start += subset_length

        if not some_complement_is_failing:
            n = min(n*2, len(sam_lines))
            if n == len(sam_lines):
                break

    return sam_lines

def ddmin_pos(sam_lines, temp_prefix, samtools, included_positions):
    #assert test(sam_lines) == "FAIL"

    n = 2     # Initial granularity
    while len(included_positions) >= 2:
        start = 0
        subset_length = len(included_positions) / n
        some_complement_is_failing = False

        while start < len(included_positions):
            complement = included_positions[:start] + included_positions[start + subset_length:]

            if test(sam_lines, temp_prefix, samtools, complement) == "FAIL":
                included_positions = complement
                n = max(n - 1, 2)
                some_complement_is_failing = True
                break
                
            start += subset_length

        if not some_complement_is_failing:
            n = min(n*2, len(included_positions))
            if n == len(included_positions):
                break

    return included_positions

def compare_conseqs(txtfilename, samtools):
    raw_count = line_counter.count(txtfilename)
    try:
        samtools_conseqs = build_conseqs_with_samtools(txtfilename, samtools, raw_count)
    except:
        samtools_conseqs = "Error"
    
    python2_conseqs = build_conseqs_with_python(txtfilename)

    if samtools_conseqs != python2_conseqs:
        # Some debugging code to compare the two results.
        with open(txtfilename, 'rU') as samfile:
            sam_lines = list(samfile)
        simple_prefix = os.path.splitext(txtfilename)[0] + '_simple'
        simple_sam_lines = ddmin(sam_lines, simple_prefix, samtools)
        max_len = 0
        for line in simple_sam_lines:
            if not line.startswith('@'):
                fields = line.split('\t')
                seq = fields[9]
                max_len = max(max_len, len(seq))
        included_positions = range(max_len)
        simple_included_positions = ddmin_pos(simple_sam_lines,
                                              simple_prefix,
                                              samtools,
                                              included_positions)
        with open(simple_prefix + '.txt', 'w') as simplefile:
            write_simple_sam(simplefile,
                             simple_sam_lines,
                             simple_included_positions)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Find the simplest test failure by trimming SAM files.')

    parser.add_argument('workdir', help='path to folder holding SAM files')
    parser.add_argument('--pattern',
                        default='temp.sam*.sam',
                        help='File name pattern to match SAM files')

    args = parser.parse_args()
    
    logger = init_logging_console_only(logging.INFO)
    samtools = externals.Samtools(settings.samtools_version,
                                  settings.samtools_path,
                                  logger)
    for txtfilename in glob.glob(os.path.join(args.workdir, args.pattern)):
        logger.info(os.path.basename(txtfilename))
        compare_conseqs(txtfilename, samtools)
    logger.info('Done.')
