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
import itertools
from operator import attrgetter
import re


def write_simple_sam(samfile, sam_lines):
    for line in sam_lines:
        samfile.write(line)

def test(sam_lines, temp_prefix, samtools):
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
            write_simple_sam(samfile, sam_lines)
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
        
        result = 'PASS' if samtools_conseqs == python2_conseqs else 'FAIL'
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

            result = test(complement, temp_prefix, samtools)
            complement_description = '0-{} + {}-{}'.format(start,
                                                           start+subset_length,
                                                           len(sam_lines))
            print len(complement), result, complement_description, txtfilename

            if result == "FAIL":
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

class SamHeader(object):
    def __init__(self, line):
        self.line = line
        self.qname = self.flag = self.pos = None
        
    def __repr__(self):
        return 'SamHeader({!r})'.format(self.line)
        
class SamBase(object):
    def __init__(self, nuc, quality, pos, token_type, qname, flag, rname, mapq):
        self.nuc = nuc
        self.quality = quality
        self.pos = pos
        self.token_type = token_type
        self.qname = qname
        self.flag = flag
        self.rname = rname
        self.mapq = mapq
    
    @classmethod
    def split_bases(cls, sam_lines):
        bases = []
        for line in sam_lines:
            if line.startswith('@'):
                bases.append(SamHeader(line.rstrip('\n')))
                continue
            fields = line.split('\t')
            qname = fields[0]
            flag = fields[1]
            rname = fields[2]
            pos = int(fields[3])
            mapq = fields[4]
            cigar = fields[5]
            token_parts = re.split('(\D)', cigar)
            token_types = []
            for token_size, token_type in zip(*[iter(token_parts)]*2):
                token_types.extend([token_type] * int(token_size))
            seq = fields[9]
            quality = fields[10]
            next_pos = pos
            seq_index = 0
            for _i, token_size, token_type in cls._split_cigar(cigar):
                if token_type not in 'MD':
                    base_pos = None
                else:
                    base_pos = next_pos
                    next_pos += 1
                if token_type == 'D':
                    nuc = q = None
                else:
                    nuc = seq[seq_index]
                    q = quality[seq_index]
                    seq_index += 1
                bases.append(cls(nuc,
                                 q,
                                 base_pos,
                                 token_type,
                                 qname,
                                 flag,
                                 rname,
                                 mapq))
                pass
        return bases
    
    @classmethod
    def _fill_positions(cls, sam_bases):
        """ Fill in the gaps in a sequence of bases.
        
        Any gaps in the pos attributes will be filled with SamBase('A', '#').
        @return: a sequence of bases without gaps
        """
        prev_pos = None
        for sam_base in sam_bases:
            next_pos = sam_base.pos
            if next_pos is not None:
                if prev_pos is not None:
                    for pos in range(prev_pos+1, next_pos):
                        yield SamBase('A',
                                      '#',
                                      pos,
                                      'M',
                                      sam_base.qname,
                                      sam_base.flag,
                                      sam_base.rname,
                                      sam_base.mapq)
                prev_pos = next_pos
            yield sam_base
    
    @classmethod
    def join(cls, sam_bases):
        sam_lines = []
        for (qname, flag), line_bases in itertools.groupby(sam_bases,
                                                           attrgetter('qname',
                                                                      'flag')):
            seq = ''
            quality = ''
            cigar = ''
            token_type = None
            token_size = 0
            fields = None
            pos = None
            for sam_base in cls._fill_positions(line_bases):
                try:
                    fields = [sam_base.line]
                except:
                    if pos == None:
                        pos = sam_base.pos
                    if token_type != sam_base.token_type:
                        if token_size > 0:
                            cigar += '{}{}'.format(token_size, token_type)
                        token_size = 1
                        token_type = sam_base.token_type
                    else:
                        token_size += 1
                    if sam_base.nuc:
                        seq += sam_base.nuc
                        quality += sam_base.quality
            if not fields:
                cigar += '{}{}'.format(token_size, token_type)
                len_str = str(len(seq))
                if flag == '147':
                    len_str = '-' + len_str
                fields = [qname,
                          flag,
                          sam_base.rname,
                          str(pos),
                          sam_base.mapq,
                          cigar,
                          '=',
                          str(pos),
                          len_str,
                          seq,
                          quality]
            sam_lines.append('\t'.join(fields) + '\n')
        return sam_lines

    @classmethod
    def _split_cigar(cls, cigar):
        """ Yield a series of items representing steps in the cigar string.
        
        @return: a generated sequence of items (i, token_size, token_type)
        """
        token_size = 0
        for c in cigar:
            if c in '0123456789':
                token_size = 10*token_size + int(c)
            else:
                token_type = c
                for i in range(token_size):
                    yield i, token_size, token_type
                token_size = 0
    
    def __repr__(self):
        return 'SamBase({!r}, {!r})'.format(self.nuc, self.quality)


def ddmin_bases(sam_bases, temp_prefix, samtools):
    #assert test(sam_bases) == "FAIL"

    n = 2     # Initial granularity
    while len(sam_bases) >= 2:
        start = 0
        subset_length = len(sam_bases) / n
        some_complement_is_failing = False

        while start < len(sam_bases):
            complement = sam_bases[:start] + sam_bases[start + subset_length:]
            
            sam_lines = SamBase.join(complement)
            if test(sam_lines, temp_prefix, samtools) == "FAIL":
                sam_bases = complement
                n = max(n - 1, 2)
                some_complement_is_failing = True
                break
                
            start += subset_length

        if not some_complement_is_failing:
            n = min(n*2, len(sam_bases))
            if n == len(sam_bases):
                break

    return sam_bases

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
        sam_bases = SamBase.split_bases(simple_sam_lines)
        simple_sam_bases = ddmin_bases(sam_bases, simple_prefix, samtools)
        really_simple_sam_lines = SamBase.join(simple_sam_bases)
        with open(simple_prefix + '.txt', 'w') as simplefile:
            write_simple_sam(simplefile, really_simple_sam_lines)

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
