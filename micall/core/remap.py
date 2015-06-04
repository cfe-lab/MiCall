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
import csv
import itertools
import logging
import os
import re
import subprocess
import sys

# These are both CodeResourceDependencies
import miseq_logging
import project_config
from micall.settings import bowtie_threads, consensus_q_cutoff,\
    max_remaps, min_mapping_efficiency
from micall.core.sam2aln import cigar_re, is_first_read
from operator import itemgetter

# SAM file format
fieldnames = [
    'qname',
    'flag',
    'rname',
    'pos',
    'mapq',
    'cigar',
    'rnext',
    'pnext',
    'tlen',
    'seq',
    'qual'
]


def resource_path(target):
    return os.path.join('' if not hasattr(sys, '_MEIPASS') else sys._MEIPASS, target)

logger = miseq_logging.init_logging_console_only(logging.DEBUG)
indel_re = re.compile('[+-][0-9]+')

def calculate_sample_name(fastq_filepath):
    filename = os.path.basename(fastq_filepath)
    return '_'.join(filename.split('_')[:2])

def is_first_read(flag):
    """
    Interpret bitwise flag from SAM field.
    Returns True or False indicating whether the read is the first read in a pair.
    """
    IS_FIRST_SEGMENT = 0x40
    return (int(flag) & IS_FIRST_SEGMENT) != 0

def is_primer(read_row, max_primer_length):
    cigar = read_row['cigar']
    sizes = re.findall(r'(\d+)M', cigar)
    match_length = max(map(int, sizes))
    return match_length <= max_primer_length


def remap(fastq1, fastq2, prelim_csv, remap_csv, remap_counts_csv, remap_conseq_csv, unmapped1, unmapped2,
          cwd=None, nthreads=None, callback=None):
    if cwd is not None:
        os.chdir(cwd)
    max_pileup_depth = str(2**16)
    
    # check that the inputs exist
    if not os.path.exists(fastq1):
        logger.error('No FASTQ found at %s', fastq1)
        sys.exit(1)

    if not os.path.exists(fastq2):
        logger.error('No FASTQ found at %s', fastq2)
        sys.exit(1)

    # check that we have access to bowtie2
    try:
        redirect_call([resource_path('bowtie2'), '-h'], os.devnull)
    except OSError:
        logger.error('bowtie2 not found; check if it is installed and in $PATH\n')
        sys.exit(1)

    # generate initial *.faidx file
    projects = project_config.ProjectConfig.loadDefault()
    ref_path = 'micall.fasta'  # TODO: this should not be hard-coded, move to settings
    with open(ref_path, 'w') as ref:
        projects.writeSeedFasta(ref)
    log_call([resource_path('samtools'), 'faidx', ref_path])

    # get the raw read count
    raw_count = count_file_lines(fastq1) / 2  # 4 lines per record in FASTQ, paired

    sample_name = calculate_sample_name(fastq1)
    remap_counts_writer = csv.DictWriter(remap_counts_csv, ['sample_name',
                                                                 'type',
                                                                 'count',
                                                                 'filtered_count'])
    remap_counts_writer.writeheader()
    remap_counts_writer.writerow(dict(sample_name=sample_name,
                                      type='raw',
                                      count=raw_count))

    # group CSV stream by first item
    with prelim_csv:
        reader = csv.DictReader(prelim_csv)
        map_counts = {}
        refgroups = {} # { group_name: (refname, count) }
        # We had a problem where a large number of primer reads would map to
        # a different seed and overwhelm the seed selection. Now we ignore
        # short reads when selecting the seed, because they might be primers.
        max_primer_length = 50
        count_threshold = 10 # must have more than this to get remapped at all
        for refname, group in itertools.groupby(reader, itemgetter('rname')):
            # reconstitute region-specific SAM files
            #tmpfile = open('%s.sam' % refname, 'w')
            count = 0
            filtered_count = 0
            for row in group:
                #tmpfile.write('\t'.join([row[field]
                #                         for field in reader.fieldnames]))
                #tmpfile.write('\n')
                count += 1
                if not is_primer(row, max_primer_length):
                    filtered_count += 1
            remap_counts_writer.writerow(dict(sample_name=sample_name,
                                              type='prelim %s' % refname,
                                              count=count,
                                              filtered_count=filtered_count))
            map_counts[refname] = count
            #tmpfile.close()
            refgroup = projects.getSeedGroup(refname)
            _best_ref, best_count = refgroups.get(refgroup,
                                                  (None, count_threshold))
            if filtered_count > best_count:
                refgroups[refgroup] = (refname, filtered_count)
    
    refnames = [refname for refname, _count in refgroups.itervalues()]

    # settings for iterative remapping
    n_remaps = 0
    frozen = []  # which regions to stop re-mapping

    samfile = open(prelim_csv.name, 'rU')  # re-open file
    conseqs = {}

    while n_remaps < max_remaps:
        if callback:
            callback(0)  # reset progress bar (standalone app only)

        # update references by generating consensus sequences from current SAM CSV
        tmpfile = open('remap-aligned.csv', 'w')
        nullfile = open(os.devnull, 'w')
        sam2aln(remap_csv=samfile, aligned_csv=tmpfile, insert_csv=nullfile, failed_csv=nullfile)

        tmpfile = open('remap-aligned.csv', 'rU')
        confile = open('remap-conseq.csv', 'w')
        aln2counts(aligned_csv=tmpfile, nuc_csv=nullfile, amino_csv=nullfile, coord_ins_csv=nullfile,
                   conseq_csv=confile, failed_align_csv=nullfile, nuc_variants_csv=nullfile)

        handle = open(confile.name, 'rU')
        outfile = open('remap-conseq.fa', 'w')
        _ = handle.next()  # header
        for line in handle:
            region, qcut, pcut, sequence = line.strip('\n').split(',')
            if pcut == 'MAX':
                conseqs.update({region: sequence})
                outfile.write('>%s\n%s\n' % (region, sequence))
        handle.close()
        outfile.close()

        # regenerate bowtie2 index files
        log_call([resource_path('bowtie2-build'), '-f', '-q', outfile.name, outfile.name])

        # stream output from bowtie2
        bowtie_args = [resource_path('bowtie2'),
                       '--quiet',
                       '-x', outfile.name,
                       '-1', fastq1,
                       '-2', fastq2,
                       '--no-unal', # don't report reads that failed to align
                       '--no-hd', # no header lines (start with @)
                       '--local',
                       '--rdg 12,3',  # increase gap open penalties
                       '--rfg 12,3',
                       '-p', str(bowtie_threads) if nthreads is None else str(nthreads)]

        # regenerate SAM CSV output
        remap_csv.close()
        remap_csv = open(remap_csv.name, 'w')  # reset file
        remap_writer = csv.DictWriter(remap_csv, fieldnames)
        remap_writer.writeheader()

        mapped = {}
        new_counts = {}
        p = subprocess.Popen(bowtie_args, stdout=subprocess.PIPE)
        with p.stdout:
            for i, line in enumerate(p.stdout):
                if callback and i%1000 == 0:
                    callback(i)
                items = line.split('\t')
                qname, bitflag, refname = items[:3]
                if not refname in new_counts:
                    new_counts.update({refname: 0})
                new_counts[refname] += 1

                if qname not in mapped:
                    mapped.update({qname: 0})
                # track how many times this read has mapped to a region with integer value
                # 0(00) = neither; 2(10) = forward only; 1(01) = reverse only; 3(11) both
                mapped[qname] += (2 if is_first_read(bitflag) else 1)

                # write SAM line to remap CSV
                remap_writer.writerow(dict(zip(fieldnames, items)))

        samfile = remap_csv  # now target this file for regenerating consensus sequences
        n_remaps += 1
        if callback:
            callback('... remap iteration %d' % n_remaps)

        # stopping criterion 1
        if all([(count <= map_counts[refname]) for refname, count in new_counts.iteritems()]):
            break

        map_counts = new_counts  # update counts

        # stopping criterion 2
        mapping_efficiency = sum(map_counts.values()) / float(raw_count)
        if mapping_efficiency > min_mapping_efficiency:
            break  # a sufficient fraction of raw data has been mapped

    ## finished iterative phase

    # write consensus sequences and counts
    remap_conseq_csv.write('region,sequence\n')  # record consensus sequences for later use
    for refname in map_counts.iterkeys():
        remap_counts_writer.writerow(dict(sample_name=sample_name,
                                          type='remap %s' % refname,
                                          count=map_counts[refname]))
        conseq = conseqs.get(refname) or projects.getReference(refname)
        remap_conseq_csv.write('%s,%s\n' % (refname, conseq))

    remap_conseq_csv.close()

    # screen raw data for reads that have not mapped to any region
    n_unmapped = 0
    with open(fastq1, 'rU') as f:
        # http://stackoverflow.com/a/1657385/4794
        # izip_longest will call f.next() four times for each iteration, so
        # the four variables are assigned the next four lines.
        for ident, seq, opt, qual in itertools.izip_longest(f, f, f, f):
            qname = ident.lstrip('@').rstrip('\n').split()[0]
            if qname not in mapped or mapped[qname] < 2:
                # forward read not mapped
                unmapped1.write(''.join([ident, seq, opt, qual]))
                n_unmapped += 1
    unmapped1.close()

    # write out the other pair
    with open(fastq2, 'rU') as f:
        for ident, seq, opt, qual in itertools.izip_longest(f, f, f, f):
            qname = ident.lstrip('@').rstrip('\n').split()[0]
            if qname not in mapped or mapped[qname] % 2 == 0:
                # reverse read not mapped
                unmapped2.write(''.join([ident, seq, opt, qual]))
                n_unmapped += 1
    unmapped2.close()

    # report number of unmapped reads
    remap_counts_writer.writerow(dict(sample_name=sample_name,
                                      type='unmapped',
                                      count=n_unmapped))


def csv_to_pileup (sam_csv):
    """
    Convert a SAM-style CSV file into samtools pileup-like format in memory.
    Process inline by region.  Also return numbers of reads mapped to each region.
    :param sam_csv: Product of prelim_map.py or iterative remapping.
    :return: dictionary by region, keying dictionaries of concatenated bases by position
    """
    pileup = {}
    counts = {}
    reader = csv.DictReader(sam_csv)
    for rcount, row in enumerate(reader):
        refname = row['rname']
        if refname not in pileup:
            pileup.update({refname: {}})
        if refname not in counts:
            counts.update({refname: 0})

        # update mapped read counts
        counts[refname] += 1

        is_first = is_first_read(row['flag'])
        seq = row['seq']
        qual = row['qual']

        pos = 0  # position in sequence
        refpos = int(row['pos']) # position in reference

        tokens = cigar_re.findall(row['cigar'])

        if tokens[0].endswith('S'):
            # skip left soft clip
            pos = int(tokens[0][:-1])
            tokens.pop(0)  # remove this first token
            print 'soft clip to query position', pos

        if not tokens[0].endswith('M'):
            # the leftmost token must end with M
            print 'ERROR: CIGAR token after soft clip must be match interval'
            sys.exit()

        # record start of read
        if refpos not in pileup[refname]:
            # keys 's' = sequence, 'q' = quality string
            pileup[refname].update({refpos: {'s': '', 'q': ''}})
        pileup[refname][refpos]['s'] += '^' + chr(int(row['mapq'])+33)

        for token in tokens:
            length = int(token[:-1])
            if token.endswith('M'):
                # match
                for i in range(length):
                    if refpos not in pileup[refname]: pileup[refname].update({refpos: {'s': '', 'q': ''}})
                    pileup[refname][refpos]['s'] += seq[pos] if is_first else seq[pos].lower()
                    pileup[refname][refpos]['q'] += qual[pos]
                    pos += 1
                    refpos += 1

            elif token.endswith('D'):
                # deletion relative to reference
                pileup[refname][refpos-1]['s'] += '-' + str(length) + ('N' if is_first else 'n')*length

                # append deletion placeholders downstream
                for i in range(refpos, refpos+length):
                    if i not in pileup[refname]: pileup[refname].update({i: {'s': '', 'q': ''}})
                    pileup[refname][i]['s'] += '*'

                refpos += length

            elif token.endswith('I'):
                # insertion relative to reference
                # FIXME: pileup does not record quality scores of inserted bases
                insert = seq[pos:(pos+length)]
                pileup[refname][refpos-1]['s'] += '+' + str(length) + (insert if is_first else insert.lower())
                pos += length

            elif token.endswith('S'):
                # soft clip
                break

            else:
                print 'ERROR: Unknown token in CIGAR string', token
                sys.exit()

        # record end of read
        pileup[refname][refpos-1]['s'] += '$'


    return pileup, counts



def pileup_to_conseq (handle, qCutoff):
    """
    Generate a consensus sequence from a samtools pileup file.
    Each line in a pileup file corresponds to a nucleotide position in the
     reference.
    Tokens are interpreted as follows:
    ^               start of read
    $               end of read
    +[1-9]+[ACGT]+  insertion relative to ref of length \1 and substring \2
    -[1-9]+N+  deletion relative to ref of length \1 and substring \2
    *               placeholder for deleted base

    FIXME: this cannot handle combinations of insertions (e.g., 1I3M2I)
    because a pileup loses all linkage information.  For now we have to
    restrict all insertions to those divisible by 3 to enforce a reading
    frame.
    """
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
                # i-th base is followed by an indel indicator
                m = indel_re.match(astr[i+1:])
                indel_len = int(m.group().strip('+-'))
                left = i+1 + len(m.group())
                insertion = astr[left:(left+indel_len)]
                q = ord(qstr[j])-33
                base = astr[i].upper() if q >= qCutoff else 'N'
                token = base + m.group() + insertion.upper()  # e.g., A+3ACG
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
            conseq += token[0]
            if int(m) % 3 == 0:
                # only add insertions that retain reading frame
                conseq += token[1+len(m):]
        elif token == '-':
            conseq += '-'
        else:
            conseq += token
    handle.close()

    # remove in-frame deletions (multiples of 3), if any
    pat = re.compile('([ACGT])(---)+([ACGT])')
    conseq = re.sub(pat, r'\g<1>\g<3>', conseq)
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


def main():
    parser = argparse.ArgumentParser(
        description='Iterative remapping of bowtie2 by reference.')

    parser.add_argument('fastq1', help='<input> FASTQ containing forward reads')
    parser.add_argument('fastq2', help='<input> FASTQ containing reverse reads')
    parser.add_argument('prelim_csv',
                        type=argparse.FileType('rU'),
                        help='<input> CSV containing preliminary map output (modified SAM)')
    parser.add_argument('remap_csv',
                        type=argparse.FileType('w'),
                        help='<output> CSV containing remap output (modified SAM)')
    parser.add_argument('remap_counts_csv',
                        type=argparse.FileType('w'),
                        help='<output> CSV containing numbers of mapped reads')
    parser.add_argument('remap_conseq_csv',
                        type=argparse.FileType('w'),
                        help='<output> CSV containing mapping consensus sequences')
    parser.add_argument('unmapped1',
                        type=argparse.FileType('w'),
                        help='<output> FASTQ R1 of reads that failed to map to any region')
    parser.add_argument('unmapped2',
                        type=argparse.FileType('w'),
                        help='<output> FASTQ R2 of reads that failed to map to any region')

    args = parser.parse_args()
    remap(fastq1=args.fastq1, fastq2=args.fastq2, prelim_csv=args.prelim_csv, remap_csv=args.remap_csv,
          remap_counts_csv=args.remap_counts_csv, remap_conseq_csv=args.remap_conseq_csv,
          unmapped1=args.unmapped1, unmapped2=args.unmapped2)

if __name__ == '__main__':
    main()
