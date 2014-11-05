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

from hyphyAlign import change_settings, pair_align
import miseq_logging
import settings
import project_config

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
    parser.add_argument('coord_ins_csv',
                        type=argparse.FileType('w'),
                        help='CSV containing insertions relative to coordinate reference')
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
                 projects,
                 conseq_mixture_cutoffs):
        """ Create an object instance.
        
        @param insert_writer: InsertionWriter object that will track reads and
            write out any insertions relative to the coordinate reference.
        @param projects: ProjectConfig object
        @param conseq_mixture_cutoffs: a list of cutoff fractions used to
            determine what portion a variant must exceed before it will be
            included as a mixture in the consensus.
        """
        self.insert_writer = insert_writer
        self.projects = projects
        self.conseq_mixture_cutoffs = list(conseq_mixture_cutoffs)
        self.conseq_mixture_cutoffs.insert(0, MAX_CUTOFF)
        
    def read(self, aligned_reads):
        """ Reset all the counters, and read a new section of aligned reads.
        
        A section must have the same sample name, region, and qcut on all lines.
        """
        self.seed_aminos = []
        self.reports = {} # {coord_name: [ReportAmino()]}
        self.inserts = set()
        self.consensus = ''
        for line in aligned_reads:
            (self.sample_name,
             self.seed,
             self.qcut,
             _rank,
             count,
             offset,
             nuc_seq) = line.rstrip().split(',')
            offset = int(offset)
            count = int(count)
            if not self.seed_aminos:
                self.insert_writer.start_group(self.sample_name,
                                              self.seed,
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
            self.coordinate_refs = {}
        else:
            self.consensus = ''.join([seed_amino.get_consensus()
                                      for seed_amino in self.seed_aminos])
            
            self.coordinate_refs = self.projects.getCoordinateReferences(self.seed)
            if not self.coordinate_refs:
                seed_ref = self.projects.getReference(self.seed)
                while len(self.seed_aminos)*3 < len(seed_ref):
                    self.seed_aminos.append(SeedAmino(len(self.seed_aminos)))
        
        for coordinate_name, coordinate_ref in self.coordinate_refs.iteritems():
            # map to reference coordinates by aligning consensus
            aquery, aref, _score = pair_align(hyphy,
                                              coordinate_ref,
                                              self.consensus)
            consensus_index = ref_index = 0
            self.inserts = set(range(len(self.consensus)))
            empty_seed_amino = SeedAmino(None)
            report_aminos = []
            is_aligned = False
            for i in range(len(aref)):
                if (consensus_index >= len(self.consensus) or
                    aquery[i] != self.consensus[consensus_index]):
                    
                    seed_amino = empty_seed_amino
                else:
                    seed_amino = self.seed_aminos[consensus_index]
                    consensus_index += 1
                if (ref_index < len(coordinate_ref) and
                    aref[i] == coordinate_ref[ref_index]):
                    
                    report_aminos.append(ReportAmino(seed_amino, ref_index+1))
                    if seed_amino.consensus_index is not None:
                        self.inserts.remove(seed_amino.consensus_index)
                        is_aligned = True
                    ref_index += 1
            if not is_aligned:
                report_aminos = []
            self.reports[coordinate_name] = report_aminos
    
    def write_amino_header(self, amino_file):
        amino_file.write(
            'sample,seed,region,q-cutoff,query.aa.pos,refseq.aa.pos,' +
            'A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,*\n')
    
    def write_amino_counts(self, amino_file):
        regions = self.reports.keys()
        regions.sort()
        for region in regions:
            for report_amino in self.reports[region]:
                seed_amino = report_amino.seed_amino
                query_pos = (str(seed_amino.consensus_index + 1)
                             if seed_amino.consensus_index is not None
                             else '')
                amino_file.write(','.join((self.sample_name,
                                           self.seed,
                                           region,
                                           self.qcut,
                                           query_pos,
                                           str(report_amino.position),
                                           seed_amino.get_report())) + '\n')

    def write_nuc_header(self, nuc_file):
        nuc_file.write(
            'sample,seed,region,q-cutoff,query.nuc.pos,refseq.nuc.pos,A,C,G,T\n')
        
    def write_nuc_counts(self, nuc_file):
        def write_counts(region, seed_amino, report_amino):
            for i, seed_nuc in enumerate(seed_amino.nucleotides):
                query_pos = (str(i + 3*seed_amino.consensus_index + 1)
                             if seed_amino.consensus_index is not None
                             else '')
                ref_pos = (str(i + 3*report_amino.position - 2)
                           if report_amino is not None
                           else '')
                nuc_file.write(','.join((self.sample_name,
                                         self.seed,
                                         region,
                                         self.qcut,
                                         query_pos,
                                         ref_pos,
                                         seed_nuc.get_report())) + '\n')
        if not self.coordinate_refs:
            for seed_amino in self.seed_aminos:
                write_counts(self.seed, seed_amino, None)
        else:
            for region, report_aminos in self.reports.iteritems():
                for report_amino in report_aminos:
                    write_counts(region, report_amino.seed_amino, report_amino)
                
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
                                                       self.seed,
                                                       self.qcut,
                                                       sample_snum,
                                                       format_cutoff(mixture_cutoff),
                                                       consensus))
    
    def write_failure_header(self, fail_file):
        fail_file.write('sample,seed,region,qcut,queryseq,refseq\n')
    
    def write_failure(self, fail_file):
        for region, report_aminos in self.reports.iteritems():
            if not report_aminos:
                coordinate_ref = self.projects.getReference(region)
                fail_file.write(','.join((self.sample_name,
                                          self.seed,
                                          region,
                                          self.qcut,
                                          self.consensus,
                                          coordinate_ref)) + '\n')
    
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
        insert_counts = {} # {left: {insert_seq: count}}
        for left, right in insert_ranges:
            current_counts = {}
            insert_counts[left] = current_counts
            for p, count in self.pcache.iteritems():
                insert_seq = p[left:right]
                current_count = current_counts.get(insert_seq, 0)
                current_counts[insert_seq] = current_count + count

        # record insertions to CSV
        for left, counts in insert_counts.iteritems():
            for insert_seq, count in counts.iteritems():
                self.insert_file.write('%s,%s,%s,%d,%s,%d\n' % (self.sample_name,
                                                                self.region,
                                                                self.qcut,
                                                                left+1,
                                                                insert_seq,
                                                                count))

def format_cutoff(cutoff):
    """ Format the cutoff fraction as a string to use as a name. """
    
    if cutoff == MAX_CUTOFF:
        return cutoff
    return '{:0.3f}'.format(cutoff)    

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
    
    projects = project_config.ProjectConfig.loadDefault()
    
    insert_writer = InsertionWriter(args.coord_ins_csv)
    
    args.aligned_csv.readline() # skip header
    
    report = SequenceReport(insert_writer,
                            projects,
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
    args.coord_ins_csv.close()
    args.failed_align_csv.close()

if __name__ == '__main__':
    main()

