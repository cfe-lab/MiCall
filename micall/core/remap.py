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
from collections import Counter, defaultdict
import csv
import itertools
import logging
from operator import itemgetter
import os
import re
import shutil
import sys
import tempfile


import miseq_logging
import project_config
from micall import settings
#from micall.core.sam2aln import cigar_re, is_first_read
from micall.utils.externals import Bowtie2, Bowtie2Build, Samtools, LineCounter

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

cigar_re = re.compile('[0-9]+[MIDNSHPX=]')  # CIGAR token

logger = miseq_logging.init_logging_console_only(logging.DEBUG)
indel_re = re.compile('[+-][0-9]+')
line_counter = LineCounter()

def is_first_read(flag):
    """
    Interpret bitwise flag from SAM field.
    Returns True or False indicating whether the read is the first read in a pair.
    """
    IS_FIRST_SEGMENT = 0x40
    return (int(flag) & IS_FIRST_SEGMENT) != 0

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
    samtools.log_call(['sort', bamfile, os.path.splitext(bamfile)[0]]) # overwrite
    
    # BAM to pileup
    # -x doesn't try to merge forward and reverse reads, because the quality
    #     calculation after merging was too hard to understand.
    # -Q0 doesn't eliminate any bases from the pileup based on quality, the
    #     quality is checked in pileup_to_conseq().
    pileup_depth = max(raw_count, 8000)
    samtools.redirect_call(['mpileup', '-xQ0', '-d', str(pileup_depth), bamfile],
                           pileup_path,
                           ignored='\[mpileup\] 1 samples in 1 input files')
    with open(pileup_path, 'rU') as f2:
        conseqs = pileup_to_conseq(f2, settings.consensus_q_cutoff)
    return conseqs


def build_conseqs_with_python(samfilename):
    with open(samfilename, 'rU') as samfile:
        return sam_to_conseqs(samfile, settings.consensus_q_cutoff)

def sam_to_conseqs(samfile, quality_cutoff=0, debug_reports=None):
    """ Build consensus sequences for each reference from a SAM file.
    
    @param samfile: an open file in the SAM format containing reads with their
        mapped position and quality scores
    @param quality_cutoff: minimum quality score for a base to be counted
    @return: {reference_name: consensus_sequence}
    """
    
    if debug_reports:
        for key in debug_reports.iterkeys():
            debug_reports[key] = Counter()
    
    # refmap structure: {refname: {pos: {nuc: count}}}
    def pos_nucs_factory():
        nuc_count_factory = Counter
        return defaultdict(nuc_count_factory)
    refmap = defaultdict(pos_nucs_factory)
    SAM_QUALITY_BASE = 33 # from SAM file format specification
    PROPERLY_ALIGNED_FLAG = 2
    quality_cutoff_char = chr(SAM_QUALITY_BASE + quality_cutoff)
    
    for read_pair in matchmaker(samfile, include_singles=True):
        read1, read2 = read_pair
        if read2 and read1[2] != read2[2]:
            # region mismatch, ignore the read pair.
            continue
        min_pos = None
        max_pos = -1
        for read in read_pair:
            if not read:
                continue
            (_qname,
             flag_str,
             rname,
             refpos_str,
             _mapq,
             cigar,
             _rnext,
             _pnext,
             _tlen,
             seq,
             qual) = read[:11] # ignore optional fields
            flag = int(flag_str)
            if not (flag & PROPERLY_ALIGNED_FLAG):
                continue
            pos_nucs = refmap[rname]
            pos = int(refpos_str)
            tokens = cigar_re.findall(cigar)
    
            token_end_pos = -1
            token_itr = iter(tokens)
            for i, nuc in enumerate(seq):
                is_in_overlap = min_pos is not None and min_pos <= pos <= max_pos
                while i > token_end_pos:
                    token = next(token_itr)
                    token_size = int(token[:-1])
                    token_type = token[-1]
                    if token_type == 'D':
                        for _ in range(token_size):
                            nuc_counts = pos_nucs[pos]
                            # report dash if all reads have deletions
                            # dash is overridden by low quality reads, so use -1
                            nuc_counts['-'] = -1
                            pos += 1
                    else:
                        token_end_pos += token_size
                    if (token_type == 'I' and
                        token_size % 3 == 0 and
                        not is_in_overlap):
                        
                        # replace previous base with insertion
                        prev_nuc_counts = pos_nucs[pos-1]
                        prev_nuc =           seq[i-1]
                        insertion =          seq[i-1:token_end_pos+1]
                        insertion_quality = qual[i-1:token_end_pos+1]
                        min_quality = min(insertion_quality)
                        if min_quality >= quality_cutoff_char:
                            prev_nuc_counts[prev_nuc] -= 1
                            prev_nuc_counts[insertion] += 1
                        
                if token_type == 'S':
                    pass
                elif token_type == 'M':
                    if not is_in_overlap:
                        nuc_counts = pos_nucs[pos]
                        if qual[i] >= quality_cutoff_char:
                            nuc_counts[nuc] += 1
                        else:
                            nuc_counts['N'] = 0
                        if min_pos is None:
                            min_pos = pos
                        max_pos = max(pos, max_pos)
                        if debug_reports:
                            counts = debug_reports.get((rname, pos))
                            if counts is not None:
                                counts[nuc + qual[i]] += 1
                    pos += 1
    if debug_reports:
        for key, counts in debug_reports.iteritems():
            mixtures = []
            nucs = set()
            qualities = set()
            for nuc, quality in counts.iterkeys():
                nucs.add(nuc)
                qualities.add(quality)
            qualities = sorted(qualities)
            for min_quality in qualities:
                filtered_counts = Counter()
                for (nuc, nuc_qual), count in counts.iteritems():
                    if nuc_qual >= min_quality:
                        filtered_counts[nuc] += count
                mixture = []
                for nuc, count in filtered_counts.iteritems():
                    mixture.append('{}: {}'.format(nuc, count))
                mixtures.append('{}{{{}}}'.format(min_quality,
                                                  ', '.join(mixture)))
            debug_reports[key] = ', '.join(mixtures)
            
    conseqs = {}
    for refname, pos_nucs in refmap.iteritems():
        if not any((any(n > 0 for n in counts.itervalues())
                    for counts in pos_nucs.itervalues())):
            #Nothing mapped, so no consensus.
            continue
        conseq = ''
        deletion = ''
        end = max(pos_nucs.keys())+1
        for pos in range(1, end):
            nuc_counts = pos_nucs[pos]
            most_common = find_top_token(nuc_counts)
            if most_common is None:
                conseq += 'N'
            elif most_common == '-':
                deletion += '-'
            else:
                if deletion:
                    if len(deletion) % 3 != 0:
                        conseq += deletion
                    deletion = ''
                conseq += most_common
        conseqs[refname] = conseq
    return conseqs

def build_conseqs(samfilename, samtools, raw_count):
    """ Build the new consensus sequences from the mapping results.
    
    @param samfilename: the mapping results in SAM format
    @param samtools: a command wrapper for samtools, or None if it's not
        available
    @param raw_count: the maximum number of reads in the SAM file
    """
    if samtools:
        conseqs = build_conseqs_with_samtools(samfilename, samtools, raw_count)
    else:
        conseqs = build_conseqs_with_python(samfilename)
    
    if False:
        # Some debugging code to save the SAM file for testing.
        copy_name = tempfile.mktemp(suffix='.sam', prefix=samfilename)
        shutil.copy(samfilename, copy_name)
        
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
          rfgopen=None,
          stderr=sys.stderr,
          gzip=False):
    """
    Iterative re-map reads from raw paired FASTQ files to a reference sequence set that
    is being updated as the consensus of the reads that were mapped to the last set.
    @param fastq1: input R1 FASTQ
    @param fastq2: input R2 FASTQ
    @param prelim_csv: input CSV output from prelim_csv()
    @param remap_csv:  output CSV, contents of bowtie2 SAM output
    @param remap_counts_csv:  output CSV, counts of reads mapped to regions
    @param remap_conseq_csv:  output CSV, sample- and region-specific consensus sequences
                                generated while remapping reads
    @param unmapped1:  output FASTQ containing R1 reads that did not map to any region
    @param unmapped2:  output FASTQ containing R2 reads that did not map to any region
    @param work_path:  optional path to store working files
    @param nthreads:  optional setting to modify the number of threads used by bowtie2
    @param callback: a function to report progress with three optional
        parameters - callback(message, progress, max_progress)
    @param count_threshold:  minimum number of reads that map to a region for it to be remapped
    @param rdgopen: read gap open penalty
    @param rfgopen: reference gap open penalty
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

    # append .gz extension if necessary
    if gzip:
        if not fastq1.endswith('.gz'):
            try:
                os.symlink(fastq1, fastq1+'.gz')
            except OSError:
                # symbolic link already exists
                pass
            fastq1 += '.gz'

        if not fastq2.endswith('.gz'):
            try:
                os.symlink(fastq2, fastq2+'.gz')
            except OSError:
                # symbolic link already exists
                pass
            fastq2 += '.gz'

    # retrieve reference sequences used for preliminary mapping
    projects = project_config.ProjectConfig.loadDefault()
    conseqs = {}
    for seed, vals in projects.config['regions'].iteritems():
        seqs = vals['reference']
        conseqs.update({str(seed): ''.join(seqs)})

    # record the raw read count
    raw_count = line_counter.count(fastq1, gzip=gzip) / 2  # 4 lines per record in FASTQ, paired

    remap_counts_writer = csv.DictWriter(remap_counts_csv,
                                         ['type', 'count', 'filtered_count'],
                                         lineterminator=os.linesep)
    remap_counts_writer.writeheader()
    remap_counts_writer.writerow(dict(type='raw', count=raw_count))

    # convert preliminary CSV to SAM, count reads
    if callback:
        callback(message='... processing preliminary map',
                 progress=0,
                 max_progress=raw_count)

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
                    callback(progress=row_count)
                
                count += 1
                row_count += 1
    
                if is_short_read(row, max_primer_length=50):
                    # exclude short reads
                    continue
    
                filtered_count += 1
    
                # write SAM row
                f.write('\t'.join([row[field] for field in fieldnames]) + '\n')
            if callback:
                callback(progress=raw_count)
            
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
    conseqs = new_conseqs


    # start remapping loop
    n_remaps = 0
    new_counts = {}
    mapped = {}
    while n_remaps < settings.max_remaps and conseqs:
        if callback:
            callback(message='... remap iteration %d' % n_remaps, progress=0)

        # generate reference file from current set of consensus sequences
        outfile = open(reffile, 'w')
        for region, conseq in conseqs.iteritems():
            outfile.write('>%s\n%s\n' % (region, conseq))
        outfile.close()

        # regenerate bowtie2 index files
        bowtie2_build.build(reffile, reffile)

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
                       #'--no-unal', # don't report reads that failed to align
                       '--no-hd', # no header lines (start with @)
                       '--local',
                       '-p', str(nthreads)]

        # capture stdout stream to count reads before writing to file
        mapped.clear()  # track which reads have mapped to something
        new_counts.clear()
        FORWARD_FLAG = 2
        REVERSE_FLAG = 1

        # reset containers with each iteration
        unmapped = {'R1': {}, 'R2': {}}

        with open(samfile, 'w') as f:
            # write SAM header
            f.write('@HD\tVN:1.0\tSO:unsorted\n')
            for rname, refseq in conseqs.iteritems():
                f.write('@SQ\tSN:%s\tLN:%d\n' % (rname, len(refseq)))
            f.write('@PG\tID:bowtie2\tPN:bowtie2\tVN:2.2.3\tCL:""\n')

            for i, line in enumerate(bowtie2.yield_output(bowtie_args, stderr=stderr)):
                if callback and i%1000 == 0:
                    callback(progress=i)  # progress monitoring in GUI

                items = line.split('\t')
                qname, bitflag, rname, _, _, _, _, _, _, seq, qual = items[:11]

                if rname == '*':
                    # did not map to any reference
                    unmapped['R1' if is_first_read(bitflag) else 'R2'].update({qname: (seq, qual)})
                    continue

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
            if callback:
                callback(progress=raw_count)

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
    remap_writer = csv.DictWriter(remap_csv, fieldnames, lineterminator=os.linesep)
    remap_writer.writeheader()
    if mapped:
        # At least one read was mapped, so samfile has relevant data
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

    # output unmapped reads to new FASTQ files
    n_unmapped = 0
    for qname, (seq, qual) in unmapped['R1'].iteritems():
        unmapped1.write('@%s\n%s\n+\n%s\n' % (qname, seq, qual))
        n_unmapped += 1
    unmapped1.close()

    for qname, (seq, qual) in unmapped['R2'].iteritems():
        unmapped2.write('@%s\n%s\n+\n%s\n' % (qname, seq, qual))
        n_unmapped += 1
    unmapped2.close()

    # report number of unmapped reads
    remap_counts_writer.writerow(dict(type='unmapped',
                                      count=n_unmapped))
    remap_counts_csv.close()


def matchmaker(samfile, include_singles=False):
    """
    An iterator that returns pairs of reads sharing a common qname from a SAM file.
    Note that unpaired reads will be left in the cached_rows dictionary and
    discarded.
    @param samfile: open file handle to a SAM file
    @param include_singles: True if unpaired reads should be returned, paired
        with a None value: ([qname, flag, rname, ...], None)
    @return: yields a tuple for each read pair with fields split by tab chars:
        ([qname, flag, rname, ...], [qname, flag, rname, ...])
    """
    ref_names = set()
    cached_rows = {}
    for line in samfile:
        row = line.strip('\n').split('\t')
        
        if line.startswith('@'):
            if row[0] == '@SQ':
                for field in row[1:]:
                    field_name, value = field.split(':', 1)
                    if field_name == 'SN':
                        ref_names.add(value)
            continue
        
        qname = row[0]
        ref_name = row[2]
        if ref_name in ref_names:
            old_row = cached_rows.pop(qname, None)
            if old_row is None:
                cached_rows[qname] = row
            else:
                # current row should be the second read of the pair
                yield old_row, row
    if include_singles:
        for row in cached_rows.itervalues():
            yield row, None

def find_top_token(base_counts):
    top_count = top_token = None
    for token, count in base_counts.most_common():
        if top_count is None:
            top_token = token
            top_count = count
        elif count < top_count:
            break
        if token < top_token:
            top_token = token
    return top_token

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
    positions_with_deletions = set()
    for line in handle:
        if to_skip > 0:
            to_skip -= 1
            continue

        fields = line.strip('\n').split('\t')
        region, pos, _refbase, _depth = fields[:4]
        if len(fields) != 6:
            astr = qstr = ''
        else:
            astr, qstr = fields[-2:]

        if region not in conseqs:
            conseqs.update({region: ''})
            #counts.update({region: 0})
            last_pos = 0

        pos = int(pos)  # position in the pileup, 1-index
        if (pos - last_pos) > 1:
            conseqs[region] += 'N' * (pos - last_pos - 1)
        last_pos = pos
        base_counts = Counter()
        if pos in positions_with_deletions:
            base_counts['-'] = -1
            positions_with_deletions.remove(pos)
        i = 0       # Current index for astr
        j = 0       # Current index for qstr

        while i < len(astr):
            if astr[i] == '^':
                # start of new read
                i += 2
            elif astr[i] in '*':
                i += 1
                j += 1
            elif astr[i] == '$':
                i += 1
            elif i < len(astr)-1 and astr[i+1] in '+-':
                # i-th base is followed by an indel indicator
                m = indel_re.match(astr[i+1:])
                indel_len = int(m.group().strip('+-'))
                left = i+1 + len(m.group())
                insertion = astr[left:(left+indel_len)]
                if astr[i+1] == '-':
                    for deletion_pos in range(pos+1, pos+indel_len+1):
                        positions_with_deletions.add(deletion_pos)
                base = astr[i].upper()
                token = base + m.group() + insertion.upper()  # e.g., A+3ACG
                q = ord(qstr[j])-33
                if q < qCutoff:
                    base_counts['N'] = 0
                else:
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

        if 'N' in base_counts or not base_counts:
            base_counts['N'] = 0
        token = find_top_token(base_counts)
        if '+' not in token:
            # add counts from all insertions to the single base
            ins_tokens = [token
                          for token in base_counts.iterkeys()
                          if '+' in token]
            for token in ins_tokens:
                base = token[0]
                base_counts[base] += base_counts[token]
                del base_counts[token]
            token = find_top_token(base_counts)

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
    trimmed_conseqs = {}
    for region, conseq in conseqs.iteritems():
        if not re.match('^[N-]*$', conseq):
            trimmed_conseqs[region] = re.sub('([ACGTN])(---)+([ACGTN])',
                                             r'\g<1>\g<3>',
                                             conseq)

    return trimmed_conseqs

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
    parser.add_argument("--gzip", help="<optional> FASTQ files are compressed",
                        action='store_true')

    args = parser.parse_args()
    remap(fastq1=args.fastq1,
          fastq2=args.fastq2,
          prelim_csv=args.prelim_csv,
          remap_csv=args.remap_csv,
          remap_counts_csv=args.remap_counts_csv,
          remap_conseq_csv=args.remap_conseq_csv,
          unmapped1=args.unmapped1,
          unmapped2=args.unmapped2,
          gzip=args.gzip)  # defaults to False

if __name__ == '__main__':
    main()
