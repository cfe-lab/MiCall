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

from hyphyAlign import change_settings, get_boundaries, pair_align
import miseq_logging
import settings

def parseArgs():
    parser = argparse.ArgumentParser(
        description='Post-processing of short-read alignments.')
    
    parser.add_argument('aligned_csv',
                        type=argparse.FileType('rU'),
                        help='aligned CSF input')
    parser.add_argument('nuc_csv',
                        type=argparse.FileType('w'),
                        help='CSV containing nucleotide frequencies')
    parser.add_argument('amino_csv',
                        type=argparse.FileType('w'),
                        help='CSV containing amino frequencies')
    parser.add_argument('indels_csv',
                        type=argparse.FileType('w'),
                        help='CSV containing insertions')
    parser.add_argument('conseq_csv',
                        type=argparse.FileType('w'),
                        help='CSV containing consensus sequences')
    parser.add_argument('failed_align_csv',
                        type=argparse.FileType('w'),
                        help='CSV containing any sequences that failed to align')
    
    return parser.parse_args()

logger = miseq_logging.init_logging_console_only(logging.DEBUG)

hyphy = HyPhy._THyPhy (os.getcwd(), 1)  # @UndefinedVariable
change_settings(hyphy)

MAX_CUTOFF = 'MAX'

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
ambig_dict = dict(("".join(sorted(v)), k)
                  for k, v in mixture_dict.iteritems()
                  if k != '-')


def translate(seq, offset=0):
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

class SequenceReport(object):
    """ Hold the data for several reports related to a sample's genetic sequence.
    
    To use a report object, read a group of aligned reads that mapped to a
    single region, and then write out all the reports for that region.
    """
    def __init__(self,
                 insert_writer,
                 seed_refs,
                 coordinate_refs,
                 conseq_mixture_cutoffs):
        """ Create an object instance.
        
        @param insert_writer: InsertionWriter object that will track reads and
            write out any insertions relative to the coordinate reference.
        @param seed_refs: {name: sequence} dictionary of seed reference
            sequences. Used to determine length of sequence when no coordinate
            reference is available.
        @param coordinate_refs: {name: sequence} dictionary of coordinate
            reference sequences.
        @param conseq_mixture_cutoffs: a list of cutoff fractions used to
            determine what portion a variant must exceed before it will be
            included as a mixture in the consensus.
        """
        self.insert_writer = insert_writer
        self.seed_refs = seed_refs
        self.coordinate_refs = coordinate_refs
        self.conseq_mixture_cutoffs = list(conseq_mixture_cutoffs)
        self.conseq_mixture_cutoffs.insert(0, MAX_CUTOFF)
        
    def read(self, aligned_reads):
        """ Reset all the counters, and read a new section of aligned reads.
        
        A section must have the same sample name, region, and qcut on all lines.
        """
        self.seed_aminos = []
        self.report_aminos = []
        self.inserts = set()
        self.consensus = ''
        for line in aligned_reads:
            (self.sample_name,
             self.region,
             self.qcut,
             _rank,
             count,
             offset,
             nuc_seq) = line.rstrip().split(',')
            offset = int(offset)
            count = int(count)
            if not self.seed_aminos:
                self.insert_writer.start_group(self.sample_name,
                                              self.region,
                                              self.qcut)

            offset_nuc_seq = '-' * offset + nuc_seq
            # pad to a codon boundary
            offset_nuc_seq += '-' * ((3 - (len(offset_nuc_seq) % 3)) % 3)
            start = offset - (offset % 3)
            for nuc_pos in range(start, len(offset_nuc_seq), 3):
                codon_index = nuc_pos / 3
                while len(self.seed_aminos) <= codon_index:
                    self.seed_aminos.append(SeedAmino(len(self.seed_aminos)))
                codon = offset_nuc_seq[nuc_pos:nuc_pos+3]
                self.seed_aminos[codon_index].count_nucleotides(codon, count)
            #TODO: avoid translating twice
            self.insert_writer.add_read(translate(offset_nuc_seq.upper()),
                                        count)
        
        if not self.seed_aminos:
            self.coordinate_ref = None
        else:
            self.consensus = ''.join([seed_amino.get_consensus()
                                      for seed_amino in self.seed_aminos])
            
            self.coordinate_ref = self.coordinate_refs.get(self.region)
            if self.coordinate_ref is None:
                seed_ref = self.seed_refs[self.region]
                while len(self.seed_aminos)*3 < len(seed_ref):
                    self.seed_aminos.append(SeedAmino(len(self.seed_aminos)))
        
        if self.coordinate_ref is not None:
            # map to reference coordinates by aligning consensus
            aquery, aref, _score = pair_align(hyphy,
                                              self.coordinate_ref,
                                              self.consensus)
            consensus_index = ref_index = 0
            self.inserts = set(range(len(self.consensus)))
            empty_seed_amino = SeedAmino(None)
            is_aligned = False
            for i in range(len(aref)):
                if (consensus_index >= len(self.consensus) or
                    aquery[i] != self.consensus[consensus_index]):
                    
                    seed_amino = empty_seed_amino
                else:
                    seed_amino = self.seed_aminos[consensus_index]
                    consensus_index += 1
                if (ref_index < len(self.coordinate_ref) and
                    aref[i] == self.coordinate_ref[ref_index]):
                    
                    self.report_aminos.append(ReportAmino(seed_amino, ref_index+1))
                    if seed_amino.consensus_index is not None:
                        self.inserts.remove(seed_amino.consensus_index)
                        is_aligned = True
                    ref_index += 1
            if not is_aligned:
                self.report_aminos = []
    
    def write_amino_header(self, amino_file):
        amino_file.write(
            'sample,region,q-cutoff,query.aa.pos,refseq.aa.pos,' +
            'A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,*\n')
    
    def write_amino_counts(self, amino_file):
        for report_amino in self.report_aminos:
            seed_amino = report_amino.seed_amino
            query_pos = (str(seed_amino.consensus_index + 1)
                         if seed_amino.consensus_index is not None
                         else '')
            amino_file.write(','.join((self.sample_name,
                                       self.region,
                                       self.qcut,
                                       query_pos,
                                       str(report_amino.position),
                                       seed_amino.get_report())) + '\n')

    def write_nuc_header(self, nuc_file):
        nuc_file.write(
            'sample,region,q-cutoff,query.nuc.pos,refseq.nuc.pos,A,C,G,T\n')
        
    def write_nuc_counts(self, nuc_file):
        def write_counts(seed_amino, report_amino):
            for i, seed_nuc in enumerate(seed_amino.nucleotides):
                query_pos = (str(i + 3*seed_amino.consensus_index + 1)
                             if seed_amino.consensus_index is not None
                             else '')
                ref_pos = (str(i + 3*report_amino.position - 2)
                           if report_amino is not None
                           else '')
                nuc_file.write(','.join((self.sample_name,
                                         self.region,
                                         self.qcut,
                                         query_pos,
                                         ref_pos,
                                         seed_nuc.get_report())) + '\n')
        if self.coordinate_ref is None:
            for seed_amino in self.seed_aminos:
                write_counts(seed_amino, None)
        else:
            for report_amino in self.report_aminos:
                write_counts(report_amino.seed_amino, report_amino)
                
    def write_consensus_header(self, conseq_file):
        conseq_file.write(
            'sample,region,q-cutoff,s-number,consensus-percent-cutoff,sequence\n')
    
    def write_consensus(self, conseq_file):
        sample_name_base, sample_snum = self.sample_name.split('_', 1)
        for mixture_cutoff in self.conseq_mixture_cutoffs:
            consensus = ''
            for seed_amino in self.seed_aminos:
                for seed_nuc in seed_amino.nucleotides:
                    consensus += seed_nuc.get_consensus(mixture_cutoff)
            conseq_file.write('%s,%s,%s,%s,%s,%s\n' % (sample_name_base,
                                                       self.region,
                                                       self.qcut,
                                                       sample_snum,
                                                       format_cutoff(mixture_cutoff),
                                                       consensus))
    
    def write_failure_header(self, fail_file):
        fail_file.write('sample,region,qcut,queryseq,refseq\n')
    
    def write_failure(self, fail_file):
        if self.coordinate_ref is not None and not self.report_aminos:
            fail_file.write(','.join((self.sample_name,
                                      self.region,
                                      self.qcut,
                                      self.consensus,
                                      self.coordinate_ref)) + '\n')
    
    def write_insertions(self):
        self.insert_writer.write(self.inserts)
    
class SeedAmino(object):
    def __init__(self, consensus_index):
        self.consensus_index = consensus_index
        self.counts = {}
        self.nucleotides = [SeedNucleotide() for _ in range(3)]
        
    def count_nucleotides(self, nuc_seq, count):
        """ Record a set of reads at this position in the seed reference.
        @param nuc_seq: a string of three nucleotides that were read at this
        position
        @param count: the number of times they were read
        """
        amino = translate(nuc_seq.upper())
        prev_count = self.counts.get(amino, 0)
        self.counts[amino] = prev_count + count
        for i in range(3):
            self.nucleotides[i].count_nucleotides(nuc_seq[i], count)
    
    def get_report(self):
        """ Build a report string with the counts of each amino acid.
        
        Report how many times each amino acid was seen in count_nucleotides().
        @return: comma-separated list of counts in the same order as the
        amino_alphabet list
        """
        return ','.join([str(self.counts.get(amino, 0))
                         for amino in amino_alphabet])
        
    def get_consensus(self):
        """ Find the amino acid that was seen most often in count_nucleotides().
        
        If there is a tie, just pick one of the tied amino acids.
        @return: the letter of the most common amino acid
        """
        max_count = 0
        consensus = None
        for amino, count in self.counts.iteritems():
            if count > max_count:
                consensus = amino
                max_count = count
        return '-' if consensus is None else consensus

class SeedNucleotide(object):
    def __init__(self):
        self.counts = {}
        
    def count_nucleotides(self, nuc_seq, count):
        """ Record a set of reads at this position in the seed reference.
        @param nuc_seq: a single nucleotide letter that was read at this
        position
        @param count: the number of times it was read
        """
        if nuc_seq == 'n':
            "Represents gap between forward and reverse read, ignore."
        else:
            prev_count = self.counts.get(nuc_seq, 0)
            self.counts[nuc_seq] = prev_count + count
    
    def get_report(self):
        """ Build a report string with the counts of each nucleotide.
        
        Report how many times each nucleotide was seen in count_nucleotides().
        @return: comma-separated list of counts for A, C, G, and T.
        """
        return ','.join(map(str, [self.counts.get(nuc, 0) for nuc in 'ACGT']))
    
    def get_consensus(self, mixture_cutoff):
        """ Choose consensus nucleotide or mixture from the counts.
        
        @param conseq_mixture_cutoffs: the minimum fraction of reads
            that a nucleotide must be found in for it to be considered,
            or MAX_CUTOFF to consider only the most common nucleotide.
        @return: The letter for the consensus nucleotide or mixture.
            Nucleotide mixtures are encoded by IUPAC symbols, and the most common
            nucleotide can be a mixture if there is a tie.
        """
        if not self.counts:
            return ''
        
        intermed = [(count, nuc) for nuc, count in self.counts.iteritems()]
        intermed.sort(reverse=True)
        
        # Remove gaps and low quality reads if there is anything else.
        for i in reversed(range(len(intermed))):
            _count, nuc = intermed[i]
            if nuc in ('N', '-') and len(intermed) > 1:
                intermed.pop(i)
        
        total_count = sum([count for count, nuc in intermed])
        mixture = []
        min_count = (intermed[0][0]
                     if mixture_cutoff == MAX_CUTOFF
                     else total_count * mixture_cutoff)
        # filter for nucleotides that pass frequency cutoff
        for count, nuc in intermed:
            if count >= min_count:
                mixture.append(nuc)

        if len(mixture) > 1:
            mixture.sort()
            consensus = ambig_dict[''.join(mixture)]
        elif len(mixture) == 1:
            # no ambiguity
            consensus = mixture[0]
        else:
            # all reads were below the cutoff
            consensus = 'N'
        return consensus

class ReportAmino(object):
    def __init__(self, seed_amino, position):
        """ Create a new instance.
        
        @param seed_amino: Counts for the 
        """
        self.seed_amino = seed_amino
        self.position = position
        
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

class InsertionWriter(object):
    def __init__(self, insert_file):
        """ Initialize a writer object.
        
        @param insert_file: an open file that the data will be written to
        """
        self.insert_file = insert_file
        self.insert_file.write('sample,region,qcut,left,insert,count\n')
    
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
            added to offset it into the consensus sequence coordinates
        @param count: the number of times this sequence was read
        """
        p = offset_sequence
        if p not in self.pcache:
            # retain linkage info for reporting insertions
            self.pcache.update({p: 0})
        self.pcache[p] += count
        
    def write(self, inserts):
        """ Write any insert ranges to the file.
        
        Sequence data comes from the reads that were added to the current group.
        @param inserts: indexes of positions in the reads that should be
            reported as insertions.
        """
        if len(inserts) == 0:
            return
        
        inserts = list(inserts)
        inserts.sort()

        # convert insertion coordinates into contiguous ranges
        insert_ranges = []
        for insert in inserts:
            if not insert_ranges or insert != insert_ranges[-1][1]:
                # just starting or we hit a gap
                insert_ranges.append([insert, insert + 1])
            else:
                insert_ranges[-1][1] += 1

        # enumerate insertions by popping out all AA sub-string variants
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
                self.insert_file.write('%s,%s,%s,%d,%s,%d\n' % (self.sample_name,
                                                                self.region,
                                                                self.qcut,
                                                                left+1,
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
    def __init__(self, nucfile, amino_ref_seqs, nuc_ref_seqs):
        """ Initialize a writer object.
        
        @param nucfile: an open file that the summary will be written to
        @param amino_ref_seqs: {region: sequence} maps from region name to 
            amino acid sequence
        @param nuc_ref_seqs: {region: sequence} maps from region name to
            nucleotide sequence
        """
        self.nucfile = nucfile
        self.amino_ref_seqs = amino_ref_seqs
        self.nuc_ref_seqs = nuc_ref_seqs
        self.nucfile.write(
            'sample,region,q-cutoff,query.nuc.pos,refseq.nuc.pos,A,C,G,T\n')
        

    def _write_counts(self, sample_name, region, qcut, counts, query_display_pos, ref_display_pos):
        self.nucfile.write('%s,%s,%s,' % (sample_name, region, qcut))
        outstr = ','.join(map(str, [counts.get(nuc, 0) for nuc in 'ACGT']))
        self.nucfile.write(','.join(map(str, (query_display_pos, ref_display_pos, outstr))))
        self.nucfile.write('\n')

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
         zero-based position within the consensus sequence and add the offset
         of where the start of the consensus sequence maps to the reference
         sequence. Each entry is a dictionary keyed by nucleotide letter and
         containing the number of times each nucleotide was read at that
         position.

        @param qindex_to_refcoord: {qindex: refcoord} maps amino acid
         coordinates from the consensus sequence (query) to the reference
         sequence. May also be None if the region has no amino reference.

        @param min_offset: the minimum nucleotide offset that was seen among
         all the reads when aligning them to the reference sequence.
        """

        no_counts = {}
        amino_ref_seq = self.amino_ref_seqs.get(region)
        if amino_ref_seq is None:
            # a region like HLA
            nuc_ref_seq = self.nuc_ref_seqs[region]
            for pos in range(len(nuc_ref_seq)):
                counts = nuc_counts.get(pos, no_counts)
                self._write_counts(sample_name,
                                   region,
                                   qcut,
                                   counts,
                                   pos + 1,
                                   '')
        else:
            # output nucleotide frequencies and consensus sequences
            last_ref_nuc_pos = -1
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
                while last_ref_nuc_pos+1 < ref_nuc_pos:
                    last_ref_nuc_pos += 1
                    self._write_counts(sample_name,
                                       region,
                                       qcut,
                                       no_counts,
                                       '',
                                       last_ref_nuc_pos+1)
                
                counts = nuc_counts[query_nuc_pos]
                self._write_counts(sample_name, region, qcut, counts, query_nuc_pos+1, ref_nuc_pos+1)
                last_ref_nuc_pos = ref_nuc_pos
                
            while last_ref_nuc_pos+1 < len(amino_ref_seq)*3:
                last_ref_nuc_pos += 1
                self._write_counts(sample_name,
                                   region,
                                   qcut,
                                   no_counts,
                                   '',
                                   last_ref_nuc_pos+1)

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

def format_cutoff(cutoff):
    """ Format the cutoff fraction as a string to use as a name. """
    
    if cutoff == MAX_CUTOFF:
        return cutoff
    return '{:0.3f}'.format(cutoff)    

def make_consensus(nuc_counts, conseq_mixture_cutoffs):
    """ Make consensus sequences from nucleotide frequencies at many cut offs.
    
    @param nuc_counts: {pos: {nuc: count}} a dictionary keyed by the
         nucleotide position. To calculate the nucleotide position, take the
         zero-based position within the consensus sequence and add the offset
         of where the start of the consensus sequence maps to the reference
         sequence. Each entry is a dictionary keyed by nucleotide letter and
         containing the number of times each nucleotide was read at that
         position.
    @param conseq_mixture_cutoffs: the minimum fraction of reads at a position
        that a nucleotide must be found in for it to be considered.
    @return: {cutoff: seq} a dictionary where each entry has a key of the
        mixture cutoff name, and the value is the consensus sequence. The
        mixture cutoff name is either 'MAX' for the consensus sequence where
        each position holds the most common nucleotide, or it is the cutoff
        fraction formatted to three decimal places. For example, '0.020'.
        Nucleotide mixtures are encoded by IUPAC symbols, and the most common
        nucleotide can be a mixture if there is a tie.
    """
    nuc_coords = nuc_counts.keys()
    nuc_coords.sort()  # iterate over positions in ascending order
    all_cutoffs = [MAX_CUTOFF]
    all_cutoffs.extend(conseq_mixture_cutoffs)
    conseqs = dict([(format_cutoff(cut), '') for cut in all_cutoffs])

    for query_nuc_pos in nuc_coords:
        intermed = [(count, nuc) for nuc, count in nuc_counts[query_nuc_pos].iteritems()]
        intermed.sort(reverse=True)
        
        # Remove gaps and low quality reads if there is anything else.
        for i in reversed(range(len(intermed))):
            _count, nuc = intermed[i]
            if nuc in ('N', '-') and len(intermed) > 1:
                intermed.pop(i)
        
        total_count = sum([count for count, nuc in intermed])
        for cut in all_cutoffs:
            mixture = []
            min_count = (intermed[0][0]
                         if cut == MAX_CUTOFF
                         else total_count * cut)
            # filter for nucleotides that pass frequency cutoff
            for count, nuc in intermed:
                if count >= min_count:
                    mixture.append(nuc)

            if len(mixture) > 1:
                mixture.sort()
                consensus = ambig_dict[''.join(mixture)]
            elif len(mixture) == 1:
                # no ambiguity
                consensus = mixture[0]
            else:
                # all reads were below the cutoff
                consensus = 'N'
            conseqs[format_cutoff(cut)] += consensus

    return conseqs

def read_references(possible_ref_files):
    # check that the reference input exists
    is_ref_found = False
    for ref_file in possible_ref_files:
        if not os.path.isfile(ref_file):
            continue
        is_ref_found = True
        break
    if not is_ref_found:
        raise RuntimeError('No reference sequences found in {!r}'.format(
            possible_ref_files))

    # read in amino acid reference sequences from file
    refseqs = {}
    with open(ref_file, 'rU') as handle:
        for line_number, line in enumerate(handle):
            stripped_line = line.rstrip('\n')
            if line_number == 0:
                is_multiline = stripped_line.startswith('>')
                if not is_multiline:
                    continue
            if not is_multiline:
                region, seq = stripped_line.split(',')
                refseqs[region] = seq
            elif line_number % 2 == 0:
                region = stripped_line[1:]
            else:
                # second line of multiline entry holds sequence
                refseqs[region] = stripped_line
    return refseqs

def main():
    args = parseArgs()
    
    
    amino_ref_seqs = read_references((
        os.path.basename(settings.final_alignment_ref_path),
        settings.final_alignment_ref_path))
    nuc_ref_file = settings.mapping_ref_path + '.fasta'
    nuc_ref_seqs = read_references((
        os.path.basename(nuc_ref_file),
        nuc_ref_file))
    
    indel_writer = InsertionWriter(args.indels_csv)
    
    args.aligned_csv.readline() # skip header
    
    report = SequenceReport(indel_writer,
                            nuc_ref_seqs,
                            amino_ref_seqs,
                            settings.conseq_mixture_cutoffs)
    report.write_amino_header(args.amino_csv)
    report.write_consensus_header(args.conseq_csv)
    report.write_failure_header(args.failed_align_csv)
    report.write_nuc_header(args.nuc_csv)
    for _key, aligned_reads in groupby(args.aligned_csv,
                                       lambda x: x.split(',')[0:3]):
        report.read(aligned_reads)
        
        report.write_amino_counts(args.amino_csv)
        report.write_consensus(args.conseq_csv)
        report.write_failure(args.failed_align_csv)
        report.write_insertions()
        report.write_nuc_counts(args.nuc_csv)

    args.aligned_csv.close()
    args.amino_csv.close()
    args.nuc_csv.close()
    args.conseq_csv.close()
    args.indels_csv.close()
    args.failed_align_csv.close()

if __name__ == '__main__':
    main()

