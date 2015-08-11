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
from datetime import datetime
import itertools
import logging
import os
import re
import subprocess
import sys

from micall import settings
from micall.core import miseq_logging
from micall.core import project_config
from micall.core.sam2aln import cigar_re, is_first_read
from micall.utils.externals import Bowtie2, Bowtie2Build, Samtools
from operator import itemgetter
from collections import Counter

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





def build_conseqs_with_samtools(samfilename, samtools, raw_count):
    bamfile = samfilename.replace('.sam', '.bam')
    pileup_path = bamfile + '.pileup'
    samtools.redirect_call(['view', '-b', samfilename], bamfile)
    samtools.log_call(['sort', bamfile, bamfile.replace('.bam', '')]) # overwrite
    
    # BAM to pileup
    pileup_depth = max(raw_count, 8000)
    samtools.redirect_call(['mpileup', '-d', str(pileup_depth), bamfile], pileup_path)
    with open(pileup_path, 'rU') as f2:
        conseqs = pileup_to_conseq(f2, settings.consensus_q_cutoff)
    return conseqs


def build_conseqs_with_python(samfilename, raw_count):
    with open(samfilename, 'rU') as samfile:
        pileup, _ = sam_to_pileup(samfile, max_primer_length=50)

    is_pileup_dumped = False
    if is_pileup_dumped:    
        with open('/tmp/temp.python.pileup', 'w') as pfile:
            for region, piles in pileup.iteritems():
                for pos, pile in piles.iteritems():
                    pfile.write('\t'.join([region,
                                           str(pos),
                                           'N',
                                           str(len(pile['q'])),
                                           pile['s'],
                                           pile['q']]) + '\n')
    # use slower internal code
    conseqs = make_consensus(pileup=pileup, qCutoff=settings.consensus_q_cutoff)
    return conseqs

def build_conseqs(samfilename, samtools, raw_count):
    """ Build the new consensus sequences from the mapping results.
    
    @param samfilename: the mapping results in SAM format
    @param samtools: a command wrapper for samtools, or None if it's not
        available
    @param raw_count: the maximum number of reads in the SAM file
    """
    use_samtools = samtools is not None
    use_python = not use_samtools
    if use_samtools:
        start = datetime.now()
        samtools_conseqs = build_conseqs_with_samtools(samfilename, samtools, raw_count)
        samtools_seconds = (datetime.now() - start).total_seconds()
        conseqs = samtools_conseqs
    
    if use_python:
        start = datetime.now()
        python_conseqs = build_conseqs_with_python(samfilename, raw_count)
        python_seconds = (datetime.now() - start).total_seconds()
        conseqs = python_conseqs
    
    if use_samtools and use_python:
        # Some debugging code to compare the two results.
        with open(samfilename, 'rU') as samfile:
            first_qname = ''
            for line in samfile:
                if not line.startswith('@'):
                    first_qname = line.split('\t')[0]
                    break
        timing_file_name = samfilename.replace('.sam', '_conseq_timing.csv')
        new_file = not os.path.exists(timing_file_name)
        with open(timing_file_name, 'ab') as timing_file:
            writer = csv.DictWriter(timing_file,
                                    ['first_qname',
                                     'lines',
                                     'samtools_seconds',
                                     'python_seconds',
                                     'samtools_conseqs',
                                     'python_conseqs'])
            if new_file:
                writer.writeheader()
            row = {'first_qname': first_qname,
                   'lines': count_file_lines(samfilename),
                   'samtools_seconds': samtools_seconds,
                   'python_seconds': python_seconds}
            if python_conseqs != samtools_conseqs:
                row['samtools_conseqs'] = samtools_conseqs
                row['python_conseqs'] = python_conseqs
            writer.writerow(row)
        
    return conseqs

def remap(fastq1,
          fastq2,
          prelim_csv,
          remap_csv,
          remap_counts_csv,
          remap_conseq_csv,
          unmapped1,
          unmapped2,
          work_path='',
          nthreads=None,
          callback=None,
          count_threshold=10,
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
    :param work_path:  optional path to store working files
    :param nthreads:  optional setting to modify the number of threads used by bowtie2
    :param callback:  optional setting to pass a callback function, used for progress
                        monitoring in GUI
    :param count_threshold:  minimum number of reads that map to a region for it to be remapped
    :param rdgopen: read gap open penalty
    :param rfgopen: reference gap open penalty
    :return:
    """

    reffile = os.path.join(work_path, 'temp.fasta')
    samfile = os.path.join(work_path, 'temp.sam')

    nthreads = nthreads or settings.bowtie_threads
    bowtie2 = Bowtie2(settings.bowtie_version, settings.bowtie_path)
    bowtie2_build = Bowtie2Build(settings.bowtie_version,
                                 settings.bowtie_build_path,
                                 logger)
    samtools = settings.samtools_version and Samtools(settings.samtools_version,
                                                      settings.samtools_path,
                                                      logger)

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
        refgroups = {} # { group_name: (refname, count) }
        reader = csv.DictReader(prelim_csv)
        row_count = 0
        for refname, group in itertools.groupby(reader, itemgetter('rname')):
            count = 0
            filtered_count = 0
            for row in group:
                if callback and row_count%1000 == 0:
                    callback(row_count)
                
                count += 1
                row_count += 1
    
                if is_short_read(row, max_primer_length=50):
                    # exclude short reads
                    continue
    
                filtered_count += 1
    
                # write SAM row
                f.write('\t'.join([row[field] for field in fieldnames]) + '\n')
            
            # report preliminary counts to file
            remap_counts_writer.writerow(
                dict(type='prelim %s' % refname,
                     count=count,
                     filtered_count=filtered_count))
            refgroup = projects.getSeedGroup(refname)
            _best_ref, best_count = refgroups.get(refgroup,
                                                  (None, count_threshold))
            if filtered_count > best_count:
                refgroups[refgroup] = (refname, filtered_count)

    seed_counts = {best_ref: best_count
                   for best_ref, best_count in refgroups.itervalues()}
    # regenerate consensus sequences based on preliminary map
    conseqs = build_conseqs(samfile, samtools, raw_count)

    # exclude references with low counts (post filtering)
    new_conseqs = {}
    map_counts = {}
    for rname, conseq in conseqs.iteritems():
        count = seed_counts.get(rname, None)
        if count is not None:
            map_counts[rname] = count  # transfer filtered counts to map counts for remap loop
            new_conseqs[rname] = conseq
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
        FORWARD_FLAG = 2
        REVERSE_FLAG = 1
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
                mapped[qname] |= (FORWARD_FLAG
                                  if is_first_read(bitflag)
                                  else REVERSE_FLAG)

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
        conseqs = build_conseqs(samfile, samtools, raw_count)
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
            if qname not in mapped or (mapped[qname] & FORWARD_FLAG) == 0:
                # forward read not mapped
                unmapped1.write(''.join([ident, seq, opt, qual]))
                n_unmapped += 1
    unmapped1.close()

    # write out the other pair
    with open(fastq2, 'rU') as f:
        for ident, seq, opt, qual in itertools.izip_longest(f, f, f, f):
            qname = ident.lstrip('@').rstrip('\n').split()[0]
            if qname not in mapped or (mapped[qname] & REVERSE_FLAG) == 0:
                # reverse read not mapped
                unmapped2.write(''.join([ident, seq, opt, qual]))
                n_unmapped += 1
    unmapped2.close()

    # report number of unmapped reads
    remap_counts_writer.writerow(dict(type='unmapped',
                                      count=n_unmapped))
    remap_counts_csv.close()


def matchmaker(samfile):
    """
    An iterator that returns pairs of reads sharing a common qname from a SAM file.
    Note that unpaired reads will be left in the cached_rows dictionary and
    discarded.
    :param samfile: open file handle to a SAM file
    :return: yields a tuple for each read pair with fields split by tab chars:
        ([qname, flag, rname, ...], [qname, flag, rname, ...])
    """
    cached_rows = {}
    for line in samfile:
        if line.startswith('@'):
            continue
        
        row = line.strip('\n').split('\t')
        qname = row[0]
        old_row = cached_rows.pop(qname, None)
        if old_row is None:
            cached_rows[qname] = row
        else:
            # current row should be the second read of the pair
            yield old_row, row


def sam_to_pileup(samfile, max_primer_length, max_count=None):
    """
    Convert SAM file into samtools pileup-like format in memory.
    Process inline by region.  Also return numbers of reads mapped to each region.
    :param samfile: Product of prelim_map.py or iterative remapping.
    :param max_primer_length: Argument to pass to is_short_read, exclude possible
        primer-derived reads.
    :param max_count: Process no more than these many reads.  If None, then do
        complete read set.
    :return: dictionary by region, keying dictionaries of concatenated bases by position
    """
    pileup = {}
    counts = {}
    rcount = 0
    for read_pair in matchmaker(samfile):
        read1, read2 = read_pair
        if read1[2] != read2[2]:
            # region mismatch, ignore the read pair.
            continue
        last_pos = -1
        for read in read_pair:
            (_qname,
             flag,
             rname,
             refpos,
             mapq,
             cigar,
             _rnext,
             _pnext,
             _tlen,
             seq,
             qual) = read[0:11] # drop optional fields.
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
            is_in_overlap = refpos <= last_pos
    
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
            if refpos > last_pos:
                pileup[rname][refpos]['s'] += '^' + chr(int(mapq)+33)
    
            for token in tokens:
                length = int(token[:-1])
                if token.endswith('M'):
                    # match
                    for i in range(length):
                        if refpos not in pileup[rname]: pileup[rname].update({refpos: {'s': '', 'q': ''}})
                        if refpos > last_pos:
                            pileup[rname][refpos]['s'] += seq[pos] if is_first else seq[pos].lower()
                            pileup[rname][refpos]['q'] += qual[pos]
                            last_pos = refpos
                            is_in_overlap = False
                        pos += 1
                        refpos += 1
    
                elif token.endswith('D'):
                    if refpos > last_pos:
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
                    if refpos > last_pos:
                        insert = seq[pos:(pos+length)]
                        pileup[rname][refpos-1]['s'] += '+' + str(length) + (insert if is_first else insert.lower())
                    pos += length
    
                elif token.endswith('S'):
                    # soft clip
                    break
    
                else:
                    print 'ERROR: Unknown token in CIGAR string', token
                    sys.exit()
    
            if not is_in_overlap:
                # record end of read
                pileup[rname][refpos-1]['s'] += '$'

    return pileup, counts



def make_consensus(pileup, qCutoff):
    """
    Convert pileup dictionary returned by csv_to_pileup() and extract a consensus sequence,
    taking into account indel polymorphisms.  Iterate over regions in the pileup.
    Intervals without coverage are represented by N's.
    :param pileup: dictionary of pileups by region
    :param qCutoff: quality cutoff for bases
    :return:
    """
    conseqs = {}
    for region, pile in pileup.iteritems():
        # in case there are intervals without coverage, iterate over entire range
        maxpos = max(pile.keys())
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
            atypes -= set('N')
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
            last_pos = 0

        pos = int(pos)  # position in the pileup, 1-index
        if (pos - last_pos) > 1:
            conseqs[region] += 'N' * (pos - last_pos - 1)
        last_pos = pos
        base_counts = Counter()
        i = 0       # Current index for astr
        j = 0       # Current index for qstr

        while i < len(astr):
            if astr[i] == '^':
                # start of new read
                #counts[region] += 1
                q = ord(qstr[j])-33
                base = astr[i+2] if q >= qCutoff else 'N'
                base_counts[base.upper()] += 1
                i += 3
                j += 1
            elif astr[i] in '*':
                base_counts['-'] += 1
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
                    base_counts[token] += 1
                else:
                    base_counts[base] += 1
                i += len(token)
                j += 1
            else:
                # Operative case: sequence matches reference (And no indel ahead)
                q = ord(qstr[j])-33
                base = astr[i].upper() if q >= qCutoff else 'N'
                base_counts[base] += 1
                i += 1
                j += 1

        base_counts['N'] = 0
        token = base_counts.most_common(1)[0][0]

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

    args = parser.parse_args()
    remap(fastq1=args.fastq1, fastq2=args.fastq2, prelim_csv=args.prelim_csv, remap_csv=args.remap_csv,
          remap_counts_csv=args.remap_counts_csv, remap_conseq_csv=args.remap_conseq_csv,
          unmapped1=args.unmapped1, unmapped2=args.unmapped2)

if __name__ == '__main__':
    main()
