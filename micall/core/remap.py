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
from tempfile import gettempdir

tmp = gettempdir()

# These are both CodeResourceDependencies
import miseq_logging
import project_config

from micall import settings
from micall.core.sam2aln import cigar_re, is_first_read
from micall.utils.externals import Bowtie2, Bowtie2Build, Samtools

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

max_pileup_depth = str(2**16)


logger = miseq_logging.init_logging_console_only(logging.DEBUG)
indel_re = re.compile('[+-][0-9]+')

def is_short_read(read_row, max_primer_length):
    """
    Return length of matched intervals in read
    :param read_row:
    :param max_primer_length:
    :return:
    """
    cigar = read_row['cigar']
    sizes = re.findall(r'(\d+)M', cigar)
    match_length = max(map(int, sizes))
    return match_length <= max_primer_length




def build_conseqs(use_samtools, samfile, samtools, conseqs):
    """ Build the new consensus sequences from the mapping results.
    
    @param use_samtools:  user has the option to use native Python code exclusively of samtools,
                            but this is slower
    @param samfile: the mapping results in SAM format
    @param samtools: a command wrapper for samtools
    @param conseqs: the consensus sequences from the previous iteration
    """
    if use_samtools: # convert SAM to BAM
        bamfile = samfile.replace('.sam', '.bam')
        pileup_path = bamfile + '.pileup'
        samtools.redirect_call(['view', '-b', samfile], bamfile)
        samtools.log_call(['sort', bamfile, bamfile.replace('.bam', '')]) # overwrite
        
        # BAM to pileup
        samtools.redirect_call(['mpileup', '-d', max_pileup_depth, bamfile], pileup_path)
        with open(pileup_path, 'rU') as f2:
            conseqs = pileup_to_conseq(f2, settings.consensus_q_cutoff)
    else:
        pileup, _ = sam_to_pileup(samfile, max_primer_length=50, max_count=1E5)
        # use slower internal code
        conseqs = make_consensus(pileup=pileup, 
                                 last_conseqs=conseqs, 
                                 qCutoff=settings.consensus_q_cutoff)
    return conseqs

def remap(fastq1,
          fastq2,
          prelim_csv,
          remap_csv,
          remap_counts_csv,
          remap_conseq_csv,
          unmapped1,
          unmapped2,
          cwd=None,
          nthreads=None,
          callback=None,
          count_threshold=10,
          use_samtools=True,
          rdgopen=None,
          rfgopen=None):
    """
    Iterative re-map reads from raw paired FASTQ files to a reference sequence set that
    is being updated as the consensus of the reads that were mapped to the last set.
    :param fastq1: input R1 FASTQ
    :param fastq2: input R2 FASTQ
    :param prelim_csv: input CSV output from prelim_csv()
    :param remap_csv:  output CSV, contents of bowtie2 SAM output
    :param remap_counts_csv:  output CSV, counts of reads mapped to regions
    :param remap_conseq_csv:  output CSV, sample- and region-specific consensus sequences
                                generated while remapping reads
    :param unmapped1:  output FASTQ containing R1 reads that did not map to any region
    :param unmapped2:  output FASTQ containing R2 reads that did not map to any region
    :param cwd:  optional setting to change current working directory
    :param nthreads:  optional setting to modify the number of threads used by bowtie2
    :param callback:  optional setting to pass a callback function, used for progress
                        monitoring in GUI
    :param count_threshold:  minimum number of reads that map to a region for it to be remapped
    :param use_samtools:  user has the option to use native Python code exclusively of samtools,
                            but this is slower
    :param rdgopen: read gap open penalty
    :param rfgopen: reference gap open penalty
    :return:
    """

    reffile = os.path.join(tmp, 'temp.fa')
    samfile = os.path.join(tmp, 'temp.sam')

    if cwd is not None:
        os.chdir(cwd)
    nthreads = nthreads or settings.bowtie_threads
    bowtie2 = Bowtie2(settings.bowtie_version, settings.bowtie_path)
    bowtie2_build = Bowtie2Build(settings.bowtie_version,
                                 settings.bowtie_build_path,
                                 logger)
    samtools = Samtools(settings.samtools_version, settings.samtools_path, logger)

    # check that the inputs exist
    if not os.path.exists(fastq1):
        logger.error('No FASTQ found at %s', fastq1)
        sys.exit(1)

    if not os.path.exists(fastq2):
        logger.error('No FASTQ found at %s', fastq2)
        sys.exit(1)

    # retrieve reference sequences used for preliminary mapping
    projects = project_config.ProjectConfig.loadDefault()
    conseqs = {}
    for seed, vals in projects.config['regions'].iteritems():
        seqs = vals['reference']
        conseqs.update({str(seed): ''.join(seqs)})

    # record the raw read count
    raw_count = count_file_lines(fastq1) / 2  # 4 lines per record in FASTQ, paired

    remap_counts_writer = csv.DictWriter(remap_counts_csv, ['type',
                                                            'count',
                                                            'filtered_count'])
    remap_counts_writer.writeheader()
    remap_counts_writer.writerow(dict(type='raw', count=raw_count))

    # convert preliminary CSV to SAM, count reads
    if callback:
        callback('... processing preliminary map')
        callback(0)

    with open(samfile, 'w') as f:
        # write SAM header
        f.write('@HD\tVN:1.0\tSO:unsorted\n')
        for rname, refseq in conseqs.iteritems():
            f.write('@SQ\tSN:%s\tLN:%d\n' % (rname, len(refseq)))
        f.write('@PG\tID:bowtie2\tPN:bowtie2\tVN:2.2.3\tCL:""\n')

        # iterate through prelim CSV and record counts, transfer rows to SAM
        map_counts = {}
        filtered_counts = {}
        reader = csv.DictReader(prelim_csv)
        for i, row in enumerate(reader):
            if callback and i%1000 == 0:
                callback(i)

            rname = row['rname']
            if rname not in map_counts:
                map_counts.update({rname: 0})
            map_counts[rname] += 1

            if is_short_read(row, max_primer_length=50):
                # exclude short reads
                continue

            if rname not in filtered_counts:
                filtered_counts.update({rname: 0})
            filtered_counts[rname] += 1

            # write SAM row
            f.write('\t'.join([row[field] for field in fieldnames]) + '\n')

    # report preliminary counts to file
    for rname, count in map_counts.iteritems():
        remap_counts_writer.writerow(dict(type='prelim %s' % rname,
                                          count=count,
                                          filtered_count=filtered_counts.get(rname, 0)))

    # regenerate consensus sequences based on preliminary map
    conseqs = build_conseqs(use_samtools, samfile, samtools, conseqs)

    # exclude references with low counts (post filtering)
    new_conseqs = {}
    for rname, conseq in conseqs.iteritems():
        count = filtered_counts.get(rname, 0)
        map_counts[rname] = count  # transfer filtered counts to map counts for remap loop
        if count < count_threshold:
            continue
        new_conseqs.update({rname: conseq})
    conseqs = dict([(k, v) for k, v in new_conseqs.iteritems()])  # deep copy


    # start remapping loop
    n_remaps = 0
    new_counts = {}
    mapped = {}
    while n_remaps < settings.max_remaps and conseqs:
        if callback:
            callback('... remap iteration %d' % n_remaps)
            callback(0)  # reset progress bar (standalone app only)

        # generate reference file from current set of consensus sequences
        outfile = open(reffile, 'w')
        for region, conseq in conseqs.iteritems():
            outfile.write('>%s\n%s\n' % (region, conseq))
        outfile.close()

        # regenerate bowtie2 index files
        bowtie2_build.log_call(['-f', '-q', reffile, reffile])

        read_gap_open_penalty = rdgopen or settings.read_gap_open_remap
        ref_gap_open_penalty = rfgopen or settings.ref_gap_open_remap

        # stream output from bowtie2
        bowtie_args = ['--quiet',
                       '-x', reffile,
                       '--rdg', "{},{}".format(read_gap_open_penalty,
                                               settings.read_gap_extend_remap),
                       '--rfg', "{},{}".format(ref_gap_open_penalty,
                                               settings.ref_gap_extend_remap),
                       '-1', fastq1,
                       '-2', fastq2,
                       '--no-unal', # don't report reads that failed to align
                       '--no-hd', # no header lines (start with @)
                       '--local',
                       '-p', str(nthreads)]

        p = bowtie2.create_process(bowtie_args, stdout=subprocess.PIPE)

        # capture stdout stream to count reads before writing to file
        mapped.clear()  # track which reads have mapped to something
        new_counts.clear()
        with open(samfile, 'w') as f:
            # write SAM header
            f.write('@HD\tVN:1.0\tSO:unsorted\n')
            for rname, refseq in conseqs.iteritems():
                f.write('@SQ\tSN:%s\tLN:%d\n' % (rname, len(refseq)))
            f.write('@PG\tID:bowtie2\tPN:bowtie2\tVN:2.2.3\tCL:""\n')

            for i, line in enumerate(p.stdout):
                if callback and i%1000 == 0:
                    callback(i)  # progress monitoring in GUI

                items = line.split('\t')
                qname, bitflag, rname = items[:3]

                if rname not in new_counts:
                    new_counts.update({rname: 0})
                new_counts[rname] += 1

                if qname not in mapped:
                    mapped.update({qname: 0})
                # track how many times this read has mapped to a region with integer value
                # 0(00) = neither; 2(10) = forward only; 1(01) = reverse only; 3(11) both
                mapped[qname] += (2 if is_first_read(bitflag) else 1)

                f.write(line)

        # stopping criterion 1 - none of the regions gained reads
        if all([(count <= map_counts[refname]) for refname, count in new_counts.iteritems()]):
            break

        # stopping criterion 2 - a sufficient fraction of raw data has been mapped
        mapping_efficiency = sum(new_counts.values()) / float(raw_count)
        if mapping_efficiency > settings.min_mapping_efficiency:
            break

        # deep copy of mapping counts
        map_counts = dict([(k, v) for k, v in new_counts.iteritems()])

        # regenerate consensus sequences
        conseqs = build_conseqs(use_samtools, samfile, samtools, conseqs)
        n_remaps += 1


    ## finished iterative phase

    # generate SAM CSV output
    remap_writer = csv.DictWriter(remap_csv, fieldnames)
    remap_writer.writeheader()
    with open(samfile, 'rU') as f:
        for line in f:
            if line.startswith('@'):
                continue  # this shouldn't happen because we set --no-hd
            items = line.strip('\n').split('\t')[:11]
            remap_writer.writerow(dict(zip(fieldnames, items)))
    remap_csv.close()

    # write consensus sequences and counts
    remap_conseq_csv.write('region,sequence\n')  # record consensus sequences for later use
    for refname in new_counts.iterkeys():
        remap_counts_writer.writerow(dict(type='remap %s' % refname,
                                          count=new_counts[refname]))
        # NOTE this is the consensus sequence to which the reads were mapped, NOT the
        # current consensus!
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
    remap_counts_writer.writerow(dict(type='unmapped',
                                      count=n_unmapped))
    remap_counts_csv.close()


def sam_to_pileup (samfile, max_primer_length, max_count=None, delimiter='\t'):
    """
    Convert SAM file into samtools pileup-like format in memory.
    Process inline by region.  Also return numbers of reads mapped to each region.
    :param sam_csv: Product of prelim_map.py or iterative remapping.
    :param max_primer_length: Argument to pass to is_short_read, exclude possible
        primer-derived reads.
    :param max_count: Process no more than these many reads.  If None, then do
        complete read set.
    :return: dictionary by region, keying dictionaries of concatenated bases by position
    """
    pileup = {}
    counts = {}
    rcount = 0
    for line in samfile:
        if line.startswith('@'):
            continue

        _qname, flag, rname, refpos, mapq, cigar, _rnext, _pnext, _tlen, seq, qual = line.strip('\n').split(delimiter)
        if rname not in pileup:
            pileup.update({rname: {}})
        if rname not in counts:
            counts.update({rname: 0})
        counts[rname] += 1
        rcount += 1

        if max_count and rcount > max_count:
            # skip sequence parsing, just get counts
            continue

        is_first = is_first_read(flag)

        pos = 0  # position in sequence
        refpos = int(refpos) # position in reference

        tokens = cigar_re.findall(cigar)

        if tokens[0].endswith('S'):
            # skip left soft clip
            pos = int(tokens[0][:-1])
            tokens.pop(0)  # remove this first token

        if not tokens[0].endswith('M'):
            # the leftmost token must end with M
            print 'ERROR: CIGAR token after soft clip must be match interval'
            sys.exit()

        # record start of read
        if refpos not in pileup[rname]:
            # keys 's' = sequence, 'q' = quality string
            pileup[rname].update({refpos: {'s': '', 'q': ''}})
        pileup[rname][refpos]['s'] += '^' + chr(int(mapq)+33)

        for token in tokens:
            length = int(token[:-1])
            if token.endswith('M'):
                # match
                for i in range(length):
                    if refpos not in pileup[rname]: pileup[rname].update({refpos: {'s': '', 'q': ''}})
                    pileup[rname][refpos]['s'] += seq[pos] if is_first else seq[pos].lower()
                    pileup[rname][refpos]['q'] += qual[pos]
                    pos += 1
                    refpos += 1

            elif token.endswith('D'):
                # deletion relative to reference
                pileup[rname][refpos-1]['s'] += '-' + str(length) + ('N' if is_first else 'n')*length

                # append deletion placeholders downstream
                for i in range(refpos, refpos+length):
                    if i not in pileup[rname]: pileup[rname].update({i: {'s': '', 'q': ''}})
                    pileup[rname][i]['s'] += '*'

                refpos += length

            elif token.endswith('I'):
                # insertion relative to reference
                # FIXME: pileup does not record quality scores of inserted bases
                insert = seq[pos:(pos+length)]
                pileup[rname][refpos-1]['s'] += '+' + str(length) + (insert if is_first else insert.lower())
                pos += length

            elif token.endswith('S'):
                # soft clip
                break

            else:
                print 'ERROR: Unknown token in CIGAR string', token
                sys.exit()

        # record end of read
        pileup[rname][refpos-1]['s'] += '$'

    return pileup, counts


def make_consensus (pileup, last_conseqs, qCutoff):
    """
    Convert pileup dictionary returned by csv_to_pileup() and extract a consensus sequence,
    taking into account indel polymorphisms.  Iterate over regions in the pileup.
    Intervals without coverage are represented by N's.
    :param pileup: dictionary of pileups by region
    :param last_conseqs: dictionary of consensus sequences used to generate pileups, by region
    :param qCutoff: quality cutoff for bases
    :return:
    """
    conseqs = {}
    for region, pile in pileup.iteritems():
        # in case there are intervals without coverage, iterate over entire range
        last_conseq = last_conseqs[region]
        maxpos = len(last_conseq)
        conseq = ''
        for refpos in xrange(1, maxpos+1):
            # note position is 1-index
            v = pile.get(refpos, None)
            if v is None:
                conseq += 'N'
                continue
            astr = v['s']
            qstr = v['q']
            alist = []  # store all bases at this position
            i = 0  # current index for astr
            j = 0  # current index for qstr

            while i < len(astr):
                if astr[i] == '^':
                    # marks start of read
                    q = ord(qstr[j])-33
                    base = astr[i+2] if q >= qCutoff else 'N'
                    alist.append(base.upper())
                    i += 3
                    j += 1
                elif astr[i] in '*':
                    # marks a deletion
                    alist.append('-')
                    i += 1
                elif astr[i] == '$':
                    # marks the end of a read
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

            # determine the most frequent character state at each position, including insertions
            atypes = set(alist)
            intermed = []
            for atype in atypes:
                intermed.append((alist.count(atype), atype))
            intermed.sort(reverse=True)

            if intermed:
                token = intermed[0][1]
            else:
                # if the list was empty, then there was no coverage
                token = 'N'

            if '+' in token:
                # insertion was most common state
                m = indel_re.findall(token)[0] # \+[0-9]+
                conseq += token[0]
                if int(m) % 3 == 0:
                    # only add insertions that retain reading frame
                    conseq += token[1+len(m):]
            elif token == '-':
                # deletion was most common state
                conseq += '-'
            else:
                conseq += token

        # remove in-frame deletions (multiples of 3), if any
        pat = re.compile('([ACGT])(---)+([ACGT])')
        conseq = re.sub(pat, r'\g<1>\g<3>', conseq)
        conseqs.update({region: conseq})

    return conseqs



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
    conseqs = {}
    #counts = {}
    to_skip = 0
    last_pos = 0
    for line in handle:
        if to_skip > 0:
            to_skip -= 1
            continue

        #label, pos, en, depth, astr, qstr
        region,      pos, _,  _,     astr, qstr = line.strip('\n').split('\t')

        if region not in conseqs:
            conseqs.update({region: ''})
            #counts.update({region: 0})

        pos = int(pos)  # position in the pileup, 1-index
        if (pos - last_pos) > 1:
            conseqs[region] += 'N' * (pos - last_pos - 1)
        last_pos = pos
        alist = []  # alist stores all bases at a given coordinate
        i = 0       # Current index for astr
        j = 0       # Current index for qstr

        while i < len(astr):
            if astr[i] == '^':
                # start of new read
                #counts[region] += 1
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
            conseqs[region] += token[0]
            if int(m) % 3 == 0:
                # only add insertions that retain reading frame
                conseqs[region] += token[1+len(m):]
        elif token == '-':
            conseqs[region] += '-'
        else:
            conseqs[region] += token

    # remove in-frame deletions (multiples of 3), if any
    for region, conseq in conseqs.iteritems():
        pat = re.compile('([ACGT])(---)+([ACGT])')
        conseqs[region] = re.sub(pat, r'\g<1>\g<3>', conseq)

    return conseqs#, counts


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
    parser.add_argument("--rdgopen", default=None, help="<optional> read gap open penalty")
    parser.add_argument("--rfgopen", default=None, help="<optional> reference gap open penalty")

    def callback(msg):
        print msg

    args = parser.parse_args()
    remap(fastq1=args.fastq1, fastq2=args.fastq2, prelim_csv=args.prelim_csv, remap_csv=args.remap_csv,
          remap_counts_csv=args.remap_counts_csv, remap_conseq_csv=args.remap_conseq_csv,
          unmapped1=args.unmapped1, unmapped2=args.unmapped2, callback=callback)

if __name__ == '__main__':
    main()
