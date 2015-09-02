import argparse
import glob
import logging
import os
from subprocess import check_call, call
from tempfile import NamedTemporaryFile

from micall.core.miseq_logging import init_logging_console_only
from micall.g2p.pssm_lib import Pssm
from micall.g2p.sam_g2p import sam_g2p


def write_simple_sam(samfile, sam_lines):
    for line in sam_lines:
        samfile.write(line)

def test(remap_lines, temp_prefix, pssm, ruby_script, delete_results=True):
    """ Calculate G2P scores using ruby_script and Python, then compare.
    
    @return: 'PASS' if the results match or it's a difference we're not
        interested in, 'FAIL' otherwise
    """
    with NamedTemporaryFile(suffix=".csv", prefix=temp_prefix, delete=True) as remap_file:
        for line in remap_lines:
            remap_file.write(line)
        remap_file.flush()
        remap_file.seek(0)
        
        filename_root = os.path.splitext(os.path.splitext(remap_file.name)[0])[0]
        nuc_filename = filename_root + ".nuc.csv"
        ruby_out_filename = filename_root + "_rbg2p.csv"
        python_out_filename = filename_root + "_pyg2p.csv"
        ruby_path = os.path.dirname(ruby_script)
        
        try:
            check_call([ruby_script, remap_file.name, nuc_filename, ruby_out_filename],
                       cwd=ruby_path)
            with open(nuc_filename, 'rU') as nuc_csv, \
                 open(python_out_filename, 'wb') as g2p_csv:
                
                sam_g2p(pssm, remap_file, nuc_csv, g2p_csv)
            
            with open(os.devnull, 'w') as devnull:
                is_diff = call(['diff', '-q', ruby_out_filename, python_out_filename],
                               stdout=devnull)
            result = 'FAIL' if is_diff else 'PASS'
            logger.info('{} lines: {}'.format(len(remap_lines), result))
            return result
        finally:
            if delete_results:
                if os.path.exists(ruby_out_filename):
                    os.remove(ruby_out_filename)
                if os.path.exists(python_out_filename):
                    os.remove(python_out_filename)

def ddmin(remap_lines, temp_prefix, pssm, ruby_script):
    #assert test(remap_lines) == "FAIL"

    header = remap_lines[:1]
    remap_lines = remap_lines[1:]
    n = 2     # Initial granularity
    while len(remap_lines) >= 2:
        start = 0
        subset_length = len(remap_lines) / n
        some_complement_is_failing = False

        while start < len(remap_lines):
            complement = remap_lines[:start] + remap_lines[start + subset_length:]

            if test(header + complement, temp_prefix, pssm, ruby_script) == "FAIL":
                remap_lines = complement
                n = max(n - 1, 2)
                some_complement_is_failing = True
                break
                
            start += subset_length

        if not some_complement_is_failing:
            n = min(n*2, len(remap_lines))
            if n == len(remap_lines):
                break

    return header + remap_lines

def compare_conseqs(txtfilename, ruby_script, pssm):
    with open(txtfilename, 'rU') as remap_file:
        remap_lines = remap_file.readlines()
    simple_prefix = os.path.splitext(txtfilename)[0] + '_simple'
    if test(remap_lines, simple_prefix, pssm, ruby_script) != 'PASS':
        simple_remap_lines = ddmin(remap_lines, simple_prefix, pssm, ruby_script)
        test(simple_remap_lines,
             simple_prefix,
             pssm,
             ruby_script,
             delete_results=False)
        if not txtfilename.endswith('.txt'):
            with open(simple_prefix + '.txt', 'w') as simplefile:
                for line in simple_remap_lines:
                    simplefile.write(line)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Find the simplest test failure by trimming SAM files.')

    parser.add_argument('workdir', help='path to folder holding SAM files')
    parser.add_argument('ruby_script', help='path to Ruby version of G2P')
    parser.add_argument('--pattern',
                        default='*.remap_head.csv',
                        help='File name pattern to match SAM files')

    args = parser.parse_args()
    
    logger = init_logging_console_only(logging.INFO)
    pssm = Pssm(path_to_lookup='../g2p/g2p_fpr.txt',
                path_to_matrix='../g2p/g2p.matrix')
    for txtfilename in sorted(glob.glob(os.path.join(args.workdir, args.pattern))):
        logger.info(os.path.basename(txtfilename))
        compare_conseqs(txtfilename, args.ruby_script, pssm)
    logger.info('Done.')
