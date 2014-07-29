#! /usr/bin/env python

"""
Shipyard-style MiSeq pipeline, step 2
Takes preliminary SAM as CSV input.  Iterative re-mapping of reads from
original FASTQ files.
Also report the number of reads mapped before and after processing.
Dependencies:
    bowtie2-build
    bowtie2-align
    samtools (with mpileup modified to take higher max per-file depth)
    settings.py
"""

import argparse
import itertools
import logging
import os
import re
import subprocess
import sys

# These are both CodeResourceDependencies
import miseq_logging
from settings import bowtie_threads, consensus_q_cutoff, max_remaps, \
    min_mapping_efficiency

logger = miseq_logging.init_logging_console_only(logging.DEBUG)
indel_re = re.compile('[+-][0-9]+')

def main():
    parser = argparse.ArgumentParser('Iterative remapping of bowtie2 by reference.')
    
    parser.add_argument('fastq1', help='<input> FASTQ containing forward reads')
    parser.add_argument('fastq2', help='<input> FASTQ containing reverse reads')
    parser.add_argument('sam_csv', help='<input> SAM output of bowtie2 in CSV format')
    parser.add_argument('output_csv', help='<output> CSV containing remap output (modified SAM)')
    parser.add_argument('stats_csv', help='<output> CSV containing numbers of mapped reads')
    parser.add_argument('conseq_csv', help='<output> CSV containing mapping consensus sequences')
    
    args = parser.parse_args()
    
    max_pileup_depth = str(2**16)
    
    # check that the inputs exist
    if not os.path.exists(args.fastq1):
        logger.error('No FASTQ found at', args.fastq1)
        sys.exit(1)

    if not os.path.exists(args.fastq2):
        logger.error('No FASTQ found at', args.fastq2)
        sys.exit(1)

    # check that we have access to bowtie2
    try:
        redirect_call(['bowtie2', '-h'], os.devnull)
    except OSError:
        logger.error('bowtie2 not found; check if it is installed and in $PATH\n')
        sys.exit(1)

    # check that the output paths are valid
    for path in [args.output_csv, args.stats_csv, args.conseq_csv]:
        output_path = os.path.split(path)[0]
        if not os.path.exists(output_path) and output_path != '':
            logger.error('Output path does not exist: %s', output_path)
            sys.exit(1)

    # generate initial *.faidx file
    is_ref_found = False
    possible_refs = ('cfe.fasta', '../reference_sequences/cfe.fasta')
    for ref in possible_refs:
        if not os.path.isfile(ref):
            continue
        is_ref_found = True
        log_call(['samtools', 'faidx', ref])
        break
    if not is_ref_found:
        raise RuntimeError('No reference sequences found in {!r}'.format(
            possible_refs))

    # get the raw read count
    raw_count = count_file_lines(args.fastq1) / 2  # 4 lines per record in FASTQ, paired

    stat_file = open(args.stats_csv, 'w')
    stat_file.write('raw,%d\n' % raw_count)

    # group CSV stream by first item
    with open(args.sam_csv, 'rU') as handle:
        prelim_count = 0
        map_counts = {}
        refnames = []
        for refname, group in itertools.groupby(handle, lambda x: x.split(',')[2]):
            refnames.append(refname)
            # reconstitute region-specific SAM files
            tmpfile = open('%s.sam' % refname, 'w')
            count = 0
            for line in group:
                tmpfile.write('\t'.join(line.split(',')))
                prelim_count += 1
                count += 1
            stat_file.write('prelim %s,%d\n' % (refname, count))
            map_counts.update({refname: count})
            tmpfile.close()

    # settings for iterative remapping
    n_remaps = 0
    mapping_efficiency = float(prelim_count) / raw_count
    frozen = []  # which regions to stop re-mapping
    tmpfile = 'temp.sam'  # temporary bowtie2-align output

    conseqs = {}

    while n_remaps < max_remaps:
        if len(frozen) == len(refnames):
            # every region is frozen
            break

        for refname in refnames:
            if refname in frozen:
                # don't attempt to re-map this region
                continue

            samfile = refname+'.sam'
            bamfile = refname+'.bam'
            confile = refname+'.conseq'

            # convert SAM to BAM
            redirect_call(['samtools', 'view', '-b', '-T', confile if refname in conseqs else ref, samfile], bamfile)

            log_call(['samtools', 'sort', bamfile, refname])  # overwrite

            # BAM to pileup
            pileup_path = bamfile+'.pileup'
            redirect_call(['samtools', 'mpileup', '-d', max_pileup_depth, '-A', bamfile], pileup_path)

            # pileup to consensus sequence
            with open(pileup_path, 'rU') as f:
                conseqs[refname] = pileup_to_conseq(f, consensus_q_cutoff)

            # generate *.faidx for later calls to samtools-view
            handle = open(refname+'.conseq', 'w')
            handle.write('>%s\n%s\n' % (refname, conseqs[refname]))
            handle.close()
            log_call(['samtools', 'faidx', confile])

            # consensus to *.bt2
            log_call(['bowtie2-build', '-c', '-q', conseqs[refname], refname])
            log_call(['bowtie2',
                      '--quiet',
                      '-p',
                      str(bowtie_threads),
                      '--local',
                      '-x',
                      refname,
                      '-1',
                      args.fastq1,
                      '-2',
                      args.fastq2,
                      '--no-unal',
                      '-S',
                      tmpfile])

            # how many reads did we map?
            count = count_file_lines(tmpfile) - 3  # ignore SAM header

            if count <= map_counts[refname]:
                # failed to improve the number of mapped reads
                frozen.append(refname)
                continue

            # overwrite previous SAM file
            os.rename(tmpfile, samfile)
            map_counts[refname] = count

        n_remaps += 1
        mapping_efficiency = sum(map_counts.values()) / float(raw_count)
        if mapping_efficiency > min_mapping_efficiency:
            break  # a sufficient fraction of raw data has been mapped


    seqfile = open(args.conseq_csv, 'w')  # record consensus sequences for later use
    outfile = open(args.output_csv, 'w')  # combine SAM files into single CSV output

    for refname in refnames:
        stat_file.write('remap %s,%d\n' % (refname, map_counts[refname]))
        seqfile.write('%s,%s\n' % (refname, conseqs[refname]))
        handle = open(refname+'.sam', 'rU')
        for line in handle:
            if line.startswith('@'):
                continue  # omit SAM header lines
            items = line.strip('\n').split('\t')[:11]
            items[2] = refname  # replace '0' due to passing conseq to bowtie2-build on cmd line
            outfile.write(','.join(items) + '\n')
        handle.close()

    outfile.close()
    stat_file.close()
    seqfile.close()

def pileup_to_conseq (handle, qCutoff):
    conseq = ''
    to_skip = 0
    last_pos = 0
    for line in handle:
        if to_skip > 0:
            to_skip -= 1
            continue

        #label, pos, en, depth, astr, qstr
        _,      pos, _,  _,     astr, qstr = line.strip('\n').split('\t')
        pos = int(pos)  # position in the pileup, 1-index
        if (pos - last_pos) > 1:
            conseq += 'N' * (pos - last_pos - 1)
        last_pos = pos
        alist = []  # alist stores all bases at a given coordinate
        i = 0       # Current index for astr
        j = 0       # Current index for qstr

        while i < len(astr):
            if astr[i] == '^':
                q = ord(qstr[j])-33
                base = astr[i+2] if q >= qCutoff else 'N'
                alist.append(base.upper())
                i += 3
                j += 1
            elif astr[i] in '*':
                alist.append('-')
                i += 1
            elif astr[i] == '$':
                i += 1
            elif i < len(astr)-1 and astr[i+1] in '+-':
                m = indel_re.match(astr[i+1:])
                indel_len = int(m.group().strip('+-'))
                left = i+1 + len(m.group())
                insertion = astr[left:(left+indel_len)]
                q = ord(qstr[j])-33
                base = astr[i].upper() if q >= qCutoff else 'N'
                token = base + m.group() + insertion.upper()
                if astr[i+1] == '+':
                    alist.append(token)
                else:
                    alist.append(base)
                i += len(token)
                j += 1
            else:
                # Operative case: sequence matches reference (And no indel ahead)
                q = ord(qstr[j])-33
                base = astr[i].upper() if q >= qCutoff else 'N'
                alist.append(base)
                i += 1
                j += 1

        atypes = set(alist)
        intermed = []
        for atype in atypes:
            intermed.append((alist.count(atype), atype))
        intermed.sort(reverse=True)
        if intermed:
            token = intermed[0][1]
        else:
            token = 'N'
        if '+' in token:
            m = indel_re.findall(token)[0] # \+[0-9]+
            conseq += token[0] + token[1+len(m):]
        elif token == '-':
            pass
        else:
            conseq += token
    handle.close()
    return conseq


def redirect_call(args, outpath, format_string='%s'):
    """ Launch a subprocess, and redirect the output to a file.
    
    Raise an exception if the return code is not zero.
    Standard error is logged to the debug logger.
    @param args: A list of arguments to pass to subprocess.Popen().
    @param outpath: a filename that stdout should be redirected to. If you 
    don't need to redirect the output, then just use subprocess.check_call().
    @param format_string: A template for the debug message that will have each
    line of standard error formatted with it.
    """
    with open(outpath, 'w') as outfile:
        p = subprocess.Popen(args, stdout=outfile, stderr=subprocess.PIPE)
        for line in p.stderr:
            logger.debug(format_string, line.rstrip())
        if p.returncode:
            raise subprocess.CalledProcessError(p.returncode, args)

def log_call(args, format_string='%s'):
    """ Launch a subprocess, and log any output to the debug logger.
    
    Raise an exception if the return code is not zero. This assumes only a
    small amount of output, and holds it all in memory before logging it.
    @param args: A list of arguments to pass to subprocess.Popen().
    @param format_string: A template for the debug message that will have each
    line of output formatted with it.
    """
    output = subprocess.check_output(args, stderr=subprocess.STDOUT)
    for line in output.splitlines():
        logger.debug(format_string, line)

def count_file_lines(path):
    """ Run the wc command to count lines in a file, as shown here:
    https://gist.github.com/zed/0ac760859e614cd03652
    """
    wc_output = subprocess.check_output(['wc', '-l', path])
    return int(wc_output.split()[0])



if __name__ == '__main__':
    main()
