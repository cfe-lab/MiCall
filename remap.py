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
from glob import glob
import itertools
import logging
import os
import re
import subprocess
import sys

# These are both CodeResourceDependencies
import miseq_logging
from settings import bowtie_threads, consensus_q_cutoff, mapping_ref_path,\
    max_remaps, min_mapping_efficiency

logger = miseq_logging.init_logging_console_only(logging.DEBUG)
indel_re = re.compile('[+-][0-9]+')

def calculate_sample_name(fastq_filepath):
    filename = os.path.basename(fastq_filepath)
    return '_'.join(filename.split('_')[:2])

def parse_bitflag (flag):
    """
    Interpret bitwise flag in SAM field as follows:

    Flag    Chr Description
    =============================================================
    0x0001  p   the read is paired in sequencing
    0x0002  P   the read is mapped in a proper pair
    0x0004  u   the query sequence itself is unmapped
    0x0008  U   the mate is unmapped
    0x0010  r   strand of the query (1 for reverse)
    0x0020  R   strand of the mate
    0x0040  1   the read is the first read in a pair
    0x0080  2   the read is the second read in a pair
    0x0100  s   the alignment is not primary
    0x0200  f   the read fails platform/vendor quality checks
    0x0400  d   the read is either a PCR or an optical duplicate
    """
    labels = ['is_paired', 'is_mapped_in_proper_pair', 'is_unmapped', 'mate_is_unmapped',
            'is_reverse', 'mate_is_reverse', 'is_first', 'is_second', 'is_secondary_alignment',
            'is_failed', 'is_duplicate']

    binstr = bin(int(flag)).replace('0b', '')
    # flip the string
    binstr = binstr[::-1]
    # if binstr length is shorter than 11, pad the right with zeroes
    for i in range(len(binstr), 11):
        binstr += '0'

    bitflags = list(binstr)
    res = {}
    for i, bit in enumerate(bitflags):
        res.update({labels[i]: bool(int(bit))})

    return (res)

def main():
    parser = argparse.ArgumentParser(
        description='Iterative remapping of bowtie2 by reference.')
    
    parser.add_argument('fastq1', help='<input> FASTQ containing forward reads')
    parser.add_argument('fastq2', help='<input> FASTQ containing reverse reads')
    parser.add_argument('prelim_csv', help='<input> CSV containing preliminary map output (modified SAM)')
    parser.add_argument('remap_csv', help='<output> CSV containing remap output (modified SAM)')
    parser.add_argument('remap_counts_csv', help='<output> CSV containing numbers of mapped reads')
    parser.add_argument('remap_conseq_csv', help='<output> CSV containing mapping consensus sequences')
    parser.add_argument('unmapped1', help='<output> FASTQ R1 of reads that failed to map to any region')
    parser.add_argument('unmapped2', help='<output> FASTQ R2 of reads that failed to map to any region')
    
    args = parser.parse_args()
    
    max_pileup_depth = str(2**16)
    
    # check that the inputs exist
    if not os.path.exists(args.fastq1):
        logger.error('No FASTQ found at %s', args.fastq1)
        sys.exit(1)

    if not os.path.exists(args.fastq2):
        logger.error('No FASTQ found at %s', args.fastq2)
        sys.exit(1)

    # check that we have access to bowtie2
    try:
        redirect_call(['bowtie2', '-h'], os.devnull)
    except OSError:
        logger.error('bowtie2 not found; check if it is installed and in $PATH\n')
        sys.exit(1)

    # check that the output paths are valid
    for path in [args.remap_csv, args.remap_counts_csv, args.remap_conseq_csv]:
        output_path = os.path.split(path)[0]
        if not os.path.exists(output_path) and output_path != '':
            logger.error('Output path does not exist: %s', output_path)
            sys.exit(1)

    # generate initial *.faidx file
    is_ref_found = False
    possible_refs = glob('*.fasta')
    if not possible_refs:
        possible_refs = [mapping_ref_path + '.fasta']
    for ref in possible_refs:
        if not os.path.isfile(ref):
            continue
        is_ref_found = True
        log_call(['samtools', 'faidx', ref])
        break
    if not is_ref_found:
        possible_refs.insert(0, '*.fasta')
        raise RuntimeError('No reference sequences found in {!r}'.format(
            possible_refs))

    # get the raw read count
    raw_count = count_file_lines(args.fastq1) / 2  # 4 lines per record in FASTQ, paired

    sample_name = calculate_sample_name(args.fastq1)
    stat_file = open(args.remap_counts_csv, 'w')
    stat_file.write('sample_name,type,count\n')
    stat_file.write('%s,raw,%d\n' % (sample_name, raw_count))

    # group CSV stream by first item
    with open(args.prelim_csv, 'rU') as handle:
        handle.readline()  # skip header
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
            stat_file.write('%s,prelim %s,%d\n' % (sample_name, refname, count))
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

    mapped = {}  # track which reads have been mapped to a region

    # generate outputs
    seqfile = open(args.remap_conseq_csv, 'w')  # record consensus sequences for later use
    outfile = open(args.remap_csv, 'w')  # combine SAM files into single CSV output
    
    seqfile.write('region,sequence\n')
    outfile.write('sample_name,qname,flag,rname,pos,mapq,cigar,rnext,pnext,tlen,seq,qual\n')

    for refname in refnames:
        stat_file.write('%s,remap %s,%d\n' % (sample_name,
                                              refname,
                                              map_counts[refname]))
        seqfile.write('%s,%s\n' % (refname, conseqs[refname]))
        # transfer contents of last SAM file to CSV
        handle = open(refname+'.sam', 'rU')
        for line in handle:
            if line.startswith('@'):
                continue  # omit SAM header lines
            items = line.strip('\n').split('\t')[:11]
            qname = items[0]
            bits = parse_bitflag(items[1])

            if qname not in mapped:
                mapped.update({qname: 0})

            # track how many times this read has mapped to a region with integer value
            # 0(00) = neither; 2(10) = forward only; 1(01) = reverse only; 3(11) both
            mapped[qname] += (2 if bits['is_first'] else 1)

            items[2] = refname  # replace '0' due to passing conseq to bowtie2-build on cmd line
            items.insert(0, sample_name)
            outfile.write(','.join(items) + '\n')
        handle.close()

    outfile.close()
    seqfile.close()

    # screen raw data for reads that have not mapped to any region
    outfile = open(args.unmapped1, 'w')
    n_unmapped = 0
    with open(args.fastq1, 'rU') as f:
        # http://stackoverflow.com/questions/1657299/how-do-i-read-two-lines-from-a-file-at-a-time-using-python
        for ident, seq, opt, qual in itertools.izip_longest(*[f]*4):
            qname = ident.lstrip('@').rstrip('\n').split()[0]
            if qname not in mapped or mapped[qname] < 2:
                # forward read not mapped
                outfile.write(''.join([ident, seq, opt, qual]))
                n_unmapped += 1
    outfile.close()

    # write out the other pair
    outfile = open(args.unmapped2, 'w')
    with open(args.fastq2, 'rU') as f:
        for ident, seq, opt, qual in itertools.izip_longest(*[f]*4):
            qname = ident.lstrip('@').rstrip('\n').split()[0]
            if qname not in mapped or mapped[qname] % 2 == 0:
                # reverse read not mapped
                outfile.write(''.join([ident, seq, opt, qual]))
                n_unmapped += 1
    outfile.close()

    # report number of unmapped reads
    stat_file.write('%s,unmapped,%d\n' % (sample_name, n_unmapped))
    stat_file.close()


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
