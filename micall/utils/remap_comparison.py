import argparse
from collections import defaultdict
import glob
import logging
import os

from micall import settings
from micall.core.miseq_logging import init_logging_console_only
from micall.core.remap import build_conseqs_with_python, \
    build_conseqs_with_samtools, line_counter, sam_to_conseqs
from micall.utils import externals
import shutil
import traceback


def build_comparison(samtools_conseqs, python2_conseqs, samfile):
    report = ""
    debug_reports = {}
    for refname, samtools_conseq in samtools_conseqs.iteritems():
        python2_conseq = python2_conseqs[refname]
        for i, (s, p) in enumerate(map(None, samtools_conseq, python2_conseq)):
            if s != p:
                debug_reports[(refname, i+1)] = None
    
    if debug_reports:
        sam_to_conseqs(samfile, debug_reports=debug_reports)
        ref_positions = defaultdict(list)
        for refname, pos in debug_reports.iterkeys():
            ref_positions[refname].append(pos)
        for refname, positions in ref_positions.iteritems():
            python2_conseq = python2_conseqs[refname]
            samtools_conseq = samtools_conseqs[refname]
            positions.sort()
            for pos in positions:
                pos_report = debug_reports[(refname, pos)]
                python2_nuc = (pos <= len(python2_conseq) and
                               python2_conseq[pos-1] or
                               'None')
                samtools_nuc = (pos <= len(samtools_conseq) and
                               samtools_conseq[pos-1] or
                               'None')
                report += '{}: {}=>{} {}\n'.format((refname, pos),
                                                   samtools_nuc,
                                                   python2_nuc,
                                                   pos_report)
    return report

def compare_conseqs(samfilename, samtools):
    raw_count = line_counter.count(samfilename)
    try:
        samtools_conseqs = build_conseqs_with_samtools(samfilename, samtools, raw_count)
    except:
        traceback.print_exc()
        return
    
    python2_conseqs = build_conseqs_with_python(samfilename)
    
    
    if python2_conseqs != samtools_conseqs:
        # Some debugging code to compare the two results.
        with open(samfilename, 'rU') as samfile:
            first_qname = ''
            while True:
                line = samfile.readline()
                if not line.startswith('@'):
                    first_qname = line.split('\t')[0]
                    samfile.seek(0)
                    break
                
            print first_qname
            print build_comparison(samtools_conseqs,
                                   python2_conseqs,
                                   samfile)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Report the differences when processing simple SAM files.')

    parser.add_argument('workdir', help='path to folder holding SAM files')
    parser.add_argument('--pattern',
                        default='temp.sam*_simple.txt',
                        help='File name pattern to match simplified files')

    args = parser.parse_args()

    logger = init_logging_console_only(logging.INFO)
    samtools = externals.Samtools(settings.samtools_version,
                                  settings.samtools_path,
                                  logger)
    for txtfilename in sorted(glob.glob(os.path.join(args.workdir, args.pattern))):
        logger.info(os.path.basename(txtfilename))
        samfilename = os.path.splitext(txtfilename)[0] + '.sam'
        shutil.copy(txtfilename, samfilename)
        try:
            compare_conseqs(samfilename, samtools)
        finally:
            os.remove(samfilename)
    logger.info('Done.')
