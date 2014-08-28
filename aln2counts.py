#! /usr/bin/env python

"""
Shipyard-style MiSeq pipeline, post-processing 1
Takes aligned CSV as input (produced by sam2aln).
Re-aligns sequence variants to lab standard references (e.g., HXB2).
Reports nucleotide and amino acid frequencies by reference coordinates.
Outputs consensus sequences in same coordinate system.
This assumes a consistent reading frame across the entire region.

Outputs nucleotide counts for HLA-B in nucleotide frequencies file.
This does not assume any reading frame (because of a frameshift in HLA-B).

Dependencies:
    settings.py
    hyphyAlign.py
"""

import argparse
import HyPhy
from itertools import groupby
import logging
import os
import re
import sys

from hyphyAlign import change_settings, get_boundaries, pair_align
import miseq_logging
import settings

def parseArgs():
    parser = argparse.ArgumentParser(
        description='Post-processing of short-read alignments.')
    
    parser.add_argument('aligned_csv', help='<input> aligned CSF input')
    parser.add_argument('nuc_csv', help='<output> CSV containing nucleotide frequencies')
    parser.add_argument('amino_csv', help='<output> CSV containing amino frequencies')
    parser.add_argument('indels_csv', help='<output> CSV containing insertions')
    parser.add_argument('conseq', help='<output> CSV containing consensus sequences')
    
    return parser.parse_args()

logger = miseq_logging.init_logging_console_only(logging.DEBUG)

hyphy = HyPhy._THyPhy (os.getcwd(), 1)  # @UndefinedVariable
change_settings(hyphy)


amino_alphabet = 'ACDEFGHIKLMNPQRSTVWY*'

codon_dict = {'TTT':'F', 'TTC':'F', 'TTA':'L', 'TTG':'L',
                'TCT':'S', 'TCC':'S', 'TCA':'S', 'TCG':'S',
                'TAT':'Y', 'TAC':'Y', 'TAA':'*', 'TAG':'*',
                'TGT':'C', 'TGC':'C', 'TGA':'*', 'TGG':'W',
                'CTT':'L', 'CTC':'L', 'CTA':'L', 'CTG':'L',
                'CCT':'P', 'CCC':'P', 'CCA':'P', 'CCG':'P',
                'CAT':'H', 'CAC':'H', 'CAA':'Q', 'CAG':'Q',
                'CGT':'R', 'CGC':'R', 'CGA':'R', 'CGG':'R',
                'ATT':'I', 'ATC':'I', 'ATA':'I', 'ATG':'M',
                'ACT':'T', 'ACC':'T', 'ACA':'T', 'ACG':'T',
                'AAT':'N', 'AAC':'N', 'AAA':'K', 'AAG':'K',
                'AGT':'S', 'AGC':'S', 'AGA':'R', 'AGG':'R',
                'GTT':'V', 'GTC':'V', 'GTA':'V', 'GTG':'V',
                'GCT':'A', 'GCC':'A', 'GCA':'A', 'GCG':'A',
                'GAT':'D', 'GAC':'D', 'GAA':'E', 'GAG':'E',
                'GGT':'G', 'GGC':'G', 'GGA':'G', 'GGG':'G',
                '---':'-', 'XXX':'?'}
mixture_regex = re.compile('[WRKYSMBDHVN-]')
mixture_dict = {'W':'AT', 'R':'AG', 'K':'GT', 'Y':'CT', 'S':'CG',
                'M':'AC', 'V':'AGC', 'H':'ATC', 'D':'ATG',
                'B':'TGC', 'N':'ATGC', '-':'ATGC'}
ambig_dict = dict(("".join(sorted(v)), k) for k, v in mixture_dict.iteritems())


def translate(seq, offset):
    seq = '-'*offset + seq
    aa_seq = ''
    for codon_site in xrange(0, len(seq), 3):
        codon = seq[codon_site:codon_site+3]
        if len(codon) < 3:
            break
        if codon.count('-') > 1 or '?' in codon:
            aa_seq += '-' if (codon == '---') else '?'
            continue

        # look for nucleotide mixtures in codon, resolve to alternative codons if found
        num_mixtures = len(mixture_regex.findall(codon))
        if num_mixtures == 0:
            aa_seq += codon_dict[codon]
        elif num_mixtures == 1:
            resolved_AAs = []
            for pos in range(3):
                if codon[pos] in mixture_dict.keys():
                    for r in mixture_dict[codon[pos]]:
                        rcodon = codon[0:pos] + r + codon[(pos+1):]
                        if codon_dict[rcodon] not in resolved_AAs:
                            resolved_AAs.append(codon_dict[rcodon])
            aa_seq += '?' if (len(resolved_AAs) > 1) else resolved_AAs[0]
        else:
            aa_seq += '?'
    return aa_seq


def coordinate_map(aquery, aref):
    """
    Generate coordinate map of query sequence to a reference.
    
    @param aquery: amino acid sequence from the sample
    @param aref: reference amino acid sequence
    @return: (qindex_to_refcoord, inserts) All indexes in the return values
        are indexes into trimmed versions of aquery and aref. They are 
        calculated by counting all non-dash characters. So the index of 'R' in
        '--AC-RA--' is 2 because we don't count the dashes.
        The first item maps indexes in aquery to the corresponding indexes in
        aref.
        The second item is a list of indexes in aquery that are not mapped to
        anything in aref due to insertions in aquery. Note this is not simply
        implied by lack of a query position as a key in the dictionary, because
        ambiguous positions are also ignored.
    """
    qindex_to_refcoord = {}
    inserts = []
    qindex = 0
    rindex = 0
    left, right = get_boundaries(aref)    # Coords of first/last non-gap character

    # For each coordinate on the reference, create a mapping to the query
    for i in range(len(aref)):
        # Do not consider parts of the query outside of the reference
        if i < left:
            qindex += 1
        elif i >= right:
            break
        # A gap in the reference is an insertion in the query which we want to skip in the mapping
        elif aref[i] == '-':
            inserts.append(qindex)  # Store insert location in query coordinate space
            qindex += 1				# Track along the query
        # If there's a gap in the query we are only effectively tracking along the pre-alignment reference
        elif aquery[i] == '-':
            rindex += 1
        elif aquery[i] == '?':
            # consensus '?' is fully ambiguous position, gap spanned by mated pair
            # this is a position we want to ignore
            qindex += 1
            rindex += 1
        # Normal case: tracking forward on both sequences
        else:
            qindex_to_refcoord[qindex] = rindex
            qindex += 1
            rindex += 1

    return qindex_to_refcoord, inserts

class IndelWriter(object):
    def __init__(self, indelfile):
        """ Initialize a writer object.
        
        @param indelfile: an open file that the data will be written to
        """
        self.indelfile = indelfile
        self.indelfile.write('sample,region,qcut,left,insert,count\n')
    
    def start_group(self, sample_name, region, qcut):
        """ Start a new group of reads.
        
        @param sample_name: the name of the sample these reads came from
        @param region: the name of the region these reads mapped to
        @param qcut: the quality cut off used for these reads
        """
        self.sample_name = sample_name
        self.region = region
        self.qcut = qcut
        self.pcache = {}
    
    def add_read(self, offset_sequence, count):
        """ Add a read to the group.
        
        @param offset_sequence: the sequence of the read that has had dashes
            added to offset it into the reference coordinates
        @param count: the number of times this sequence was read
        """
        p = offset_sequence
        if p not in self.pcache:
            # retain linkage info for reporting insertions
            self.pcache.update({p: 0})
        self.pcache[p] += count
        
    def write(self, inserts, min_offset=0):
        """ Write any insert ranges to the file.
        
        Sequence data comes from the reads that were added to the current group.
        @param inserts: indexes into the non-blank characters of the reads. For
            example, if all the reads had at least two leading dashes, then you
            will need to add two to all the insert indexes to find the correct
            position within the reads.
        @param min_offset: the minimum offset from all the reads. Add this to
            the indexes in inserts to find the actual indexes of the inserted
            characters in the read strings.
        """
        if len(inserts) == 0:
            return

        # convert insertion coordinates into contiguous ranges
        insert_ranges = []
        for insert in inserts:
            adjusted = insert + min_offset
            if not insert_ranges or adjusted != insert_ranges[-1][1]:
                # just starting or we hit a gap
                insert_ranges.append([adjusted, adjusted + 1])
            else:
                insert_ranges[-1][1] += 1

        # enumerate indels by popping out all AA sub-string variants
        indel_counts = {} # {left: {insert_seq: count}}
        for left, right in insert_ranges:
            current_counts = {}
            indel_counts[left] = current_counts
            for p, count in self.pcache.iteritems():
                insert_seq = p[left:right]
                current_count = current_counts.get(insert_seq, 0)
                current_counts[insert_seq] = current_count + count

        # record insertions to CSV
        for left, counts in indel_counts.iteritems():
            for insert_seq, count in counts.iteritems():
                self.indelfile.write('%s,%s,%s,%d,%s,%d\n' % (self.sample_name,
                                                              self.region,
                                                              self.qcut,
                                                              left,
                                                              insert_seq,
                                                              count))

class AminoFrequencyWriter(object):
    def __init__(self, aafile, refseqs):
        """ Initialize a writer object.
        
        @param aafile: an open file that the summary will be written to
        @param refseqs: {region: sequence} maps from region name to amino acid 
        sequence
        """
        self.aafile = aafile
        self.refseqs = refseqs
        self.aafile.write(
            'sample,region,q-cutoff,query.aa.pos,refseq.aa.pos,'
            'A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,*\n')
        
    def write(self,
              sample_name,
              region,
              qcut,
              qindex_to_refcoord,
              amino_counts,
              inserts):
        """ Write a summary of the amino acid distribution at each position in
        the reference sequence.
        
        @param sample_name: the name of the sample to go in the first column.
        @param region: the name of the region being processed
        @param qcut: the quality cutoff score that was used to map the consensus
        sequence
        @param qindex_to_refcoord: {qindex: refcoord} maps coordinates from the
        consensus sequence (query) to the reference sequence
        @param amino_counts: {aa_pos: {aa: count}} a dictionary keyed by the amino
        acid position. To calculate the amino acid position, take the position
        within the consensus sequence and add the offset of where the start of the
        consensus sequence maps to the reference sequence. Each entry is a
        dictionary keyed by amino acid letter and containing the number of times
        each amino acid was read at that position.
        @param inserts: a list of indexes for positions in the consensus sequence
        that are inserts relative to the reference sequence
        """

        if not amino_counts:
            return
    
        qindex = 0
        query_offset = min(amino_counts.keys())
        query_end = max(qindex_to_refcoord.keys())
          
        for refcoord in range(len(self.refseqs[region])):
            while True:
                # Skip over any inserts
                aa_pos = qindex + query_offset
                if aa_pos not in inserts:
                    break
                assert aa_pos not in qindex_to_refcoord, aa_pos
                qindex += 1
            
            while True:
                # Skip over any gaps in the query consensus sequence
                mapped_coord = qindex_to_refcoord.get(qindex)
                if mapped_coord is not None or qindex >= query_end:
                    break
                qindex += 1
     
            if mapped_coord == refcoord:
                counts = [amino_counts[aa_pos].get(aa, 0)
                          for aa in amino_alphabet]
                aa_pos_str = str(aa_pos)
                qindex += 1
            else:
                counts = [0 for aa in amino_alphabet]
                aa_pos_str = ''
            outstr = ','.join(map(str, counts))
            self.aafile.write('%s,%s,%s,%s,%d,%s\n' % (sample_name,
                                                       region,
                                                       qcut,
                                                       aa_pos_str,
                                                       refcoord + 1,
                                                       outstr))

class NucleotideFrequencyWriter(object):
    def __init__(self, nucfile, amino_ref_seqs):
        """ Initialize a writer object.
        
        @param nucfile: an open file that the summary will be written to
        @param refseqs: {region: sequence} maps from region name to amino acid 
        sequence
        """
        self.nucfile = nucfile
        self.amino_ref_seqs = amino_ref_seqs
        self.nucfile.write(
            'sample,region,q-cutoff,query.nuc.pos,refseq.nuc.pos,A,C,G,T\n')
        
    def write(self,
              sample_name,
              region,
              qcut,
              nuc_counts,
              qindex_to_refcoord=None,
              min_offset=0):
        """ Write a summary of the nucleotide distribution at each position in
        the reference sequence.
        
        If there is no reference sequence for the region, as in HLA regions,
        just use the positions in the query sequence.
        
        @param sample_name: the name of the sample to go in the first column.
        @param region: the name of the region being processed
        @param qcut: the quality cutoff score that was used to map the consensus
         sequence

        @param nuc_counts: {pos: {nuc: count}} a dictionary keyed by the
         nucleotide position. To calculate the nucleotide position, take the
         position within the consensus sequence and add the offset of where the
         start of the consensus sequence maps to the reference sequence. Each
         entry is a dictionary keyed by nucleotide letter and containing the
         number of times each nucleotide was read at that position.

        @param qindex_to_refcoord: {qindex: refcoord} maps amino acid
         coordinates from the consensus sequence (query) to the reference
         sequence. May also be None if the region has no amino reference.

        @param min_offset: the minimum nucleotide offset that was seen among
         all the reads when aligning them to the reference sequence.
        """

        if region not in self.amino_ref_seqs:
            # a region like HLA
            nuc_coords = nuc_counts.keys()
            left = min(nuc_coords)
            right = max(nuc_coords)
            for pos in range(left, right+1):
                self.nucfile.write('%s,%s,%s,' % (sample_name, region, qcut))
                if pos in nuc_counts:
                    outstr = ','.join(map(str, [nuc_counts[pos].get(nuc, 0) for nuc in 'ACGT']))
                else:
                    outstr = '0,0,0,0'
                self.nucfile.write(','.join(map(str, [pos+1, '', outstr])))
                self.nucfile.write('\n')
        else:
            # output nucleotide frequencies and consensus sequences
            nuc_coords = nuc_counts.keys()
            nuc_coords.sort()

            for query_nuc_pos in nuc_coords:
                # query AA numbering starts at first complete codon, 0-index
                query_aa_pos = query_nuc_pos/3 - min_offset/3
                try:
                    ref_aa_pos = qindex_to_refcoord[query_aa_pos]  # retrieve AA position from coordinate map
                except:
                    logger.warn(
                        'No coordinate mapping for query nuc %d (amino %d) in %s' % (
                            query_nuc_pos,
                            query_aa_pos,
                            sample_name))
                    continue
                ref_nuc_pos = 3 * ref_aa_pos + (query_nuc_pos % 3)

                outstr = ','.join(map(str, [nuc_counts[query_nuc_pos].get(nuc, 0) for nuc in 'ACGT']))
                self.nucfile.write('%s,%s,%s,' % (sample_name, region, qcut))
                self.nucfile.write(','.join(map(str, [query_nuc_pos+1, ref_nuc_pos+1, outstr])))
                self.nucfile.write('\n')

def make_counts(region, group2, indel_writer):
    """
    Generate nucleotide and amino acid frequencies from CSV file grouping
     by region and quality cutoff.

    @param region: the region these lines are from
    @param group2: the lines from the CSV file that describe the reads
    @param indel_writer: this needs to receive each read so it can report the
        contents of any insertions.
    @return (nuc_counts, amino_counts, min_offset) a tuple with the following:
        nuc_counts is a dictionary of dictionaries {position: {nucleotide: count}}
            Its coordinate system is the sample consensus (nucleotide).
            If the leftmost read is downstream of the start of the reading frame,
            the smallest key of nuc_counts will be non-zero.
        amino_counts is the same, but for amino acids {position: {acid: count}}
            Its coordinate system is the sample consensus (amino).
        min_offset is the minimum offset seen in all the reads. The positions
            in nuc_counts and amino_counts already include the offsets.
    """
    nuc_counts = {}
    amino_counts = {}

    min_offset = 1e6
    total_count = 0  # track the total number of aligned and merged reads, given QCUT

    for line in group2:
        _, _, _, _, count, offset, seq = line.strip('\n').split(',')
        offset = int(offset)
        count = int(count)
        total_count += count
        seq = seq.upper()

        if offset < min_offset:
            min_offset = offset

        for i, nuc in enumerate(seq):
            pos = offset + i  # position relative to start of sample consensus (reading frame)
            if pos not in nuc_counts:
                nuc_counts.update({pos: {}})
            if nuc not in nuc_counts[pos]:
                nuc_counts[pos].update({nuc: 0})
            nuc_counts[pos][nuc] += count

        if region.startswith('HLA-'):
            continue

        p = translate('-'*offset + seq, 0)
        indel_writer.add_read(p, count)

        left = offset / 3
        right = (offset + len(seq)) / 3

        for pos in range(left, right):
            if pos not in amino_counts:
                amino_counts.update({pos: {}})
            aa = p[pos]
            if aa not in amino_counts[pos]:
                amino_counts[pos].update({aa: 0})
            amino_counts[pos][aa] += count

    return nuc_counts, amino_counts, min_offset


def make_consensus(nuc_counts):
    """
    Generate consensus sequences from nucleotide frequency data
    at varying frequency cutoffs.  Returns dict keyed by
    frequency cutoff for determining consensus.
    MAX = plurality-rule consensus.
    Nucleotide mixtures are encoded by IUPAC symbols.
    """
    nuc_coords = nuc_counts.keys()
    nuc_coords.sort()  # iterate over positions in ascending order
    conseqs = dict([(cut, '') for cut in settings.conseq_mixture_cutoffs])
    conseqs.update({'MAX': ''})

    for query_nuc_pos in nuc_coords:
        intermed = [(count, nuc) for nuc, count in nuc_counts[query_nuc_pos].iteritems()]
        intermed.sort(reverse=True)
        conseqs['MAX'] += intermed[0][1]  # plurality consensus

        total_count = sum([count for count, nuc in intermed])
        for cut in settings.conseq_mixture_cutoffs:
            mixture = []
            # filter for nucleotides that pass frequency cutoff
            for count, nuc in intermed:
                if float(count) / total_count > cut:
                    mixture.append(nuc)

            if 'N' in mixture:
                if len(mixture) > 1:
                    # N always overruled by any nucleotide (gap handled below)
                    mixture.remove('N')
                else:
                    conseqs[cut] += 'N'
                    continue

            if '-' in mixture:
                if len(mixture) > 1:
                    # gap always overruled by nucleotide
                    mixture.remove('-')
                else:
                    conseqs[cut] += '-'
                    continue

            if len(mixture) > 1:
                mixture.sort()
                conseqs[cut] += ambig_dict[''.join(mixture)]
            elif len(mixture) == 1:
                # no ambiguity
                conseqs[cut] += mixture[0]
            else:
                # no characters left - I don't think this outcome is possible
                conseqs[cut] += 'N'

    return conseqs


def main():
    args = parseArgs()
    
    # check that the amino acid reference input exists
    is_ref_found = False
    
    possible_refs = (os.path.basename(settings.final_alignment_ref_path),
                     settings.final_alignment_ref_path)
    for amino_ref in possible_refs:
        if not os.path.isfile(amino_ref):
            continue
        is_ref_found = True
        break
    if not is_ref_found:
        raise RuntimeError('No reference sequences found in {!r}'.format(
            possible_refs))

    # read in amino acid reference sequences from file
    refseqs = {}
    with open(amino_ref, 'rU') as handle:
        handle.next() # Skip header
        for line in handle:
            region, aaseq = line.strip('\n').split(',')
            refseqs.update({region: aaseq})

    # check that the inputs exist
    if not os.path.exists(args.aligned_csv):
        logger.error('No input CSF found at ' + args.aligned_csv)
        sys.exit(1)

    # check that the output paths are valid
    for path in [args.nuc_csv, args.amino_csv, args.indels_csv, args.conseq]:
        output_path = os.path.split(path)[0]
        if not os.path.exists(output_path) and output_path != '':
            logger.error('Output path does not exist: ' + output_path)
            sys.exit(1)

    # for each region, quality cutoff
    infile = open(args.aligned_csv, 'rU')
    aafile = open(args.amino_csv, 'w')
    nucfile = open(args.nuc_csv, 'w')
    confile = open(args.conseq, 'w')
    indelfile = open(args.indels_csv, 'w')
    amino_writer = AminoFrequencyWriter(aafile, refseqs)
    nuc_writer = NucleotideFrequencyWriter(nucfile, refseqs)
    indel_writer = IndelWriter(indelfile)
    
    infile.readline() # skip header
    confile.write('sample,region,q-cutoff,s-number,consensus-percent-cutoff,sequence\n')
    
    for key, group in groupby(infile, lambda x: x.split(',')[0:2]):
        sample_name, region = key
        sample_name_base, sample_snum = sample_name.split('_', 1)
        
        if region not in refseqs and not region.startswith('HLA-'):
            continue
        for qcut, group2 in groupby(group, lambda x: x.split(',')[2]):
            # gather nucleotide, amino frequencies
            indel_writer.start_group(sample_name, region, qcut)
            nuc_counts, amino_counts, min_offset = make_counts(region,
                                                               group2,
                                                               indel_writer)

            if region not in refseqs:
                nuc_writer.write(sample_name, region, qcut, nuc_counts)
            else:
                # make amino consensus from frequencies
                aa_coords = amino_counts.keys()
                aa_coords.sort()
                aa_max = ''
                for pos in range(min(aa_coords), max(aa_coords) +1):
                    if pos in aa_coords:
                        intermed = [(count, amino) for amino, count in amino_counts[pos].iteritems()]
                        intermed.sort(reverse=True)
                        aa_max += intermed[0][1]
                        continue
                    aa_max += '?'  # no coverage but not a gap


                # map to reference coordinates by aligning consensus
                aquery, aref, _ = pair_align(hyphy, refseqs[region], aa_max)
                qindex_to_refcoord, inserts = coordinate_map(aquery, aref)

                amino_writer.write(sample_name,
                                   region,
                                   qcut,
                                   qindex_to_refcoord,
                                   amino_counts,
                                   inserts)
    
                indel_writer.write(inserts, min_offset)

                nuc_writer.write(sample_name,
                                 region,
                                 qcut,
                                 nuc_counts,
                                 qindex_to_refcoord,
                                 min_offset)

            # output consensus sequences
            conseqs = make_consensus(nuc_counts)
            for cut, conseq in conseqs.iteritems():
                confile.write('%s,%s,%s,%s,%s,%s\n' % (sample_name_base,
                                                          region,
                                                          qcut,
                                                          sample_snum,
                                                          str(cut),
                                                          conseq))

    infile.close()
    aafile.close()
    indelfile.close()

if __name__ == '__main__':
    main()

