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
from collections import Counter

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
                        help='CSV containing any consensus that failed to align')
    parser.add_argument('nuc_variants_csv',
                        type=argparse.FileType('w'),
                        help='CSV containing top nucleotide variants')
    
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
        

    def _count_reads(self, aligned_reads):
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
                for reading_frame in range(3):
                    self.seed_aminos[reading_frame] = []
            self.insert_writer.add_nuc_read('-'*offset + nuc_seq, count)
                    
            for reading_frame, frame_seed_aminos in self.seed_aminos.iteritems():
                offset_nuc_seq = '-' * (reading_frame + offset) + nuc_seq
                # pad to a codon boundary
                offset_nuc_seq += '-' * ((3 - (len(offset_nuc_seq) % 3)) % 3)
                start = offset - (offset % 3)
                for nuc_pos in range(start, len(offset_nuc_seq), 3):
                    codon_index = nuc_pos / 3
                    while len(frame_seed_aminos) <= codon_index:
                        frame_seed_aminos.append(SeedAmino(len(frame_seed_aminos)))
                    
                    codon = offset_nuc_seq[nuc_pos:nuc_pos + 3]
                    frame_seed_aminos[codon_index].count_nucleotides(codon, count)
                
    def _pair_align(self, reference, query):
        """ Align a query sequence to a reference sequence using hyphy.
        
        @return: (aligned_query, aligned_ref, score) Note that this is the
            opposite order from the input parameters, but that's the order
            hyphy returns.
        """
        return pair_align(hyphy, reference, query)

    def _map_to_coordinate_ref(self, coordinate_name, coordinate_ref):
        # Start max_score with the minimum score we can consider a valid
        # alignment. Anything worse, we won't bother
        consensus_length = len([amino for amino in self.seed_aminos[0] if amino.counts])
        max_score = min(consensus_length, len(coordinate_ref))
        best_alignment = None # (consensus, aquery, aref)
        for reading_frame, frame_seed_aminos in self.seed_aminos.iteritems():
            consensus = ''.join([seed_amino1.get_consensus()
                                for seed_amino1 in frame_seed_aminos])
            if reading_frame == 0:
                # best guess before aligning
                self.consensus[coordinate_name] = consensus
            
            # map to reference coordinates by aligning consensus
            aquery, aref, score = self._pair_align(coordinate_ref, consensus)
            if score < max_score:
                continue
            max_score = score
            best_alignment = (reading_frame,
                              consensus,
                              aquery,
                              aref,
                              frame_seed_aminos)
        
        report_aminos = []
        if best_alignment is not None:
            (reading_frame,
             consensus,
             aquery,
             aref,
             frame_seed_aminos) = best_alignment
            self.reading_frames[coordinate_name] = reading_frame
            self.consensus[coordinate_name] = consensus
            consensus_index = ref_index = 0
            coordinate_inserts = set(range(len(consensus)))
            self.inserts[coordinate_name] = coordinate_inserts
            empty_seed_amino = SeedAmino(None)
            for i in range(len(aref)):
                if (consensus_index >= len(consensus) or 
                    aquery[i] != consensus[consensus_index]):
                    seed_amino = empty_seed_amino
                else:
                    seed_amino = frame_seed_aminos[consensus_index]
                    consensus_index += 1
                if (ref_index < len(coordinate_ref) and 
                    aref[i] == coordinate_ref[ref_index]):
                    report_aminos.append(ReportAmino(seed_amino, ref_index + 1))
                    if seed_amino.consensus_index is not None:
                        coordinate_inserts.remove(seed_amino.consensus_index)
                    ref_index += 1
        
        self.reports[coordinate_name] = report_aminos

    def read(self, aligned_reads):
        """ Reset all the counters, and read a new section of aligned reads.
        
        A section must have the same sample name, region, and qcut on all lines.
        """
        aligned_reads = list(aligned_reads) # let us run multiple passes
        self.seed_aminos = {} # {reading_frame: [SeedAmino(consensus_index)]}
        self.reports = {} # {coord_name: [ReportAmino()]}
        self.reading_frames = {} # {coord_name: reading_frame}
        self.inserts = {} # {coord_name: set([consensus_index])}
        self.consensus = {} # {coord_name: consensus_amino_seq}
        self.variants = {} # {coord_name: [(count, nuc_seq)]}
        self._count_reads(aligned_reads)
        
        if not self.seed_aminos:
            self.coordinate_refs = {}
        else:
            self.coordinate_refs = self.projects.getCoordinateReferences(self.seed)
            if not self.coordinate_refs:
                seed_ref = self.projects.getReference(self.seed)
                while len(self.seed_aminos[0])*3 < len(seed_ref):
                    self.seed_aminos[0].append(SeedAmino(len(self.seed_aminos[0])))
        
        for coordinate_name, coordinate_ref in self.coordinate_refs.iteritems():
            self._map_to_coordinate_ref(coordinate_name, coordinate_ref)
            report_aminos = self.reports[coordinate_name]
            max_variants = self.projects.getMaxVariants(coordinate_name)
            if report_aminos and max_variants:
                variant_counts = Counter() # {seq: count}
                for report_amino in report_aminos:
                    first_amino_index = report_amino.seed_amino.consensus_index
                    if first_amino_index is not None:
                        break
                first_amino_index = first_amino_index or 0
                start_pos = first_amino_index * 3
                for report_amino in reversed(report_aminos):
                    last_amino_index = report_amino.seed_amino.consensus_index
                    if last_amino_index is not None:
                        break
                last_amino_index = last_amino_index or -1
                end_pos = (last_amino_index+1) * 3
                minimum_variant_length = len(coordinate_ref)/2
                for line in aligned_reads:
                    (_sample_name,
                     _seed,
                     _qcut,
                     _rank,
                     count,
                     offset,
                     nuc_seq) = line.rstrip().split(',')
                    count = int(count)
                    offset = int(offset)
                    padded_seq = offset*'-' + nuc_seq
                    clipped_seq = padded_seq[start_pos:end_pos]
                    stripped_seq = clipped_seq.replace('-', '')
                    if len(stripped_seq) > minimum_variant_length:
                        variant_counts[clipped_seq] += count
                coordinate_variants = [(count, seq)
                                       for seq, count in variant_counts.iteritems()]
                coordinate_variants.sort(reverse=True)
                self.variants[coordinate_name] = coordinate_variants[0:max_variants]
    
    def write_amino_header(self, amino_file):
        amino_file.write(
            'sample,seed,region,q-cutoff,query.aa.pos,refseq.aa.pos,' +
            ','.join(amino_alphabet) + '\n')
    
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
            for seed_amino in self.seed_aminos[0]:
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
            for seed_amino in self.seed_aminos[0]:
                for seed_nuc in seed_amino.nucleotides:
                    consensus += seed_nuc.get_consensus(mixture_cutoff)
            conseq_file.write('%s,%s,%s,%s,%s,%s\n' % (sample_name_base,
                                                       self.seed,
                                                       self.qcut,
                                                       sample_snum,
                                                       format_cutoff(mixture_cutoff),
                                                       consensus))
    
    def write_nuc_variants_header(self, nuc_variants_file):
        nuc_variants_file.write('sample,seed,qcut,region,index,count,seq\n')
    
    def write_nuc_variants(self, nuc_variants_file):
        regions = self.variants.keys()
        regions.sort()
        for coordinate_name in regions:
            for i, variant in enumerate(self.variants[coordinate_name]):
                count, nuc_seq = variant
                nuc_variants_file.write(','.join((self.sample_name,
                                                  self.seed,
                                                  self.qcut,
                                                  coordinate_name,
                                                  str(i),
                                                  str(count),
                                                  nuc_seq)) + '\n')
    
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
                                          self.consensus[region],
                                          coordinate_ref)) + '\n')
    
    def write_insertions(self):
        for coordinate_name, coordinate_inserts in self.inserts.iteritems():
            self.insert_writer.write(coordinate_inserts,
                                     coordinate_name,
                                     self.reading_frames[coordinate_name])
    
class SeedAmino(object):
    def __init__(self, consensus_index):
        self.consensus_index = consensus_index
        self.counts = Counter()
        self.nucleotides = [SeedNucleotide() for _ in range(3)]
        
    def count_nucleotides(self, nuc_seq, count):
        """ Record a set of reads at this position in the seed reference.
        @param nuc_seq: a string of three nucleotides that were read at this
        position
        @param count: the number of times they were read
        """
        amino = translate(nuc_seq.upper())
        if amino in amino_alphabet:
            self.counts[amino] += count
        for i in range(3):
            self.nucleotides[i].count_nucleotides(nuc_seq[i], count)
    
    def get_report(self):
        """ Build a report string with the counts of each amino acid.
        
        Report how many times each amino acid was seen in count_nucleotides().
        @return: comma-separated list of counts in the same order as the
        amino_alphabet list
        """
        return ','.join([str(self.counts[amino])
                         for amino in amino_alphabet])
        
    def get_consensus(self):
        """ Find the amino acid that was seen most often in count_nucleotides().
        
        If there is a tie, just pick one of the tied amino acids.
        @return: the letter of the most common amino acid
        """
        consensus = self.counts.most_common(1)
        return '-' if not consensus else consensus[0][0]

class SeedNucleotide(object):
    def __init__(self):
        self.counts = Counter()
        
    def count_nucleotides(self, nuc_seq, count):
        """ Record a set of reads at this position in the seed reference.
        @param nuc_seq: a single nucleotide letter that was read at this
        position
        @param count: the number of times it was read
        """
        if nuc_seq == 'n':
            "Represents gap between forward and reverse read, ignore."
        else:
            self.counts[nuc_seq] += count
    
    def get_report(self):
        """ Build a report string with the counts of each nucleotide.
        
        Report how many times each nucleotide was seen in count_nucleotides().
        @return: comma-separated list of counts for A, C, G, and T.
        """
        return ','.join(map(str, [self.counts[nuc] for nuc in 'ACGT']))
    
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
        
        intermed = self.counts.most_common()
        
        # Remove gaps and low quality reads if there is anything else.
        for i in reversed(range(len(intermed))):
            nuc, _count = intermed[i]
            if nuc in ('N', '-') and len(intermed) > 1:
                intermed.pop(i)
        
        total_count = sum(self.counts.values())
        mixture = []
        min_count = (intermed[0][1]
                     if mixture_cutoff == MAX_CUTOFF
                     else total_count * mixture_cutoff)
        # filter for nucleotides that pass frequency cutoff
        for nuc, count in intermed:
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
        self.insert_file.write('sample,seed,region,qcut,left,insert,count\n')
    
    def start_group(self, sample_name, seed, qcut):
        """ Start a new group of reads.
        
        @param sample_name: the name of the sample these reads came from
        @param seed: the name of the region these reads mapped to
        @param qcut: the quality cut off used for these reads
        """
        self.sample_name = sample_name
        self.seed = seed
        self.qcut = qcut
        self.nuc_seqs = Counter() # {nuc_seq: count}
    
    def add_nuc_read(self, offset_sequence, count):
        """ Add a read to the group.
        
        @param offset_sequence: the nucleotide sequence of the read that has
            had dashes added to offset it into the consensus sequence
            coordinates
        @param count: the number of times this sequence was read
        """
        self.nuc_seqs[offset_sequence] += count
        
    def write(self, inserts, region, reading_frame=0):
        """ Write any insert ranges to the file.
        
        Sequence data comes from the reads that were added to the current group.
        @param inserts: indexes of positions in the reads that should be
            reported as insertions.
        @param region: the name of the coordinate region the current group was
            mapped to
        @param reading_frame: the reading frame to use when translating
            nucleotide sequences to amino acids - an integer from 0 to 2.
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
            current_counts = Counter()
            insert_counts[left] = current_counts
            for nuc_seq, count in self.nuc_seqs.iteritems():
                framed_nuc_seq = reading_frame * '-' + nuc_seq
                insert_nuc_seq = framed_nuc_seq[left*3:right*3]
                is_valid = (insert_nuc_seq and
                            'n' not in insert_nuc_seq and
                            '-' not in insert_nuc_seq)
                if is_valid:
                    insert_amino_seq = translate(insert_nuc_seq)
                    if insert_amino_seq:
                        current_counts[insert_amino_seq] += count

        # record insertions to CSV
        for left, counts in insert_counts.iteritems():
            for insert_seq, count in counts.iteritems():
                self.insert_file.write('%s,%s,%s,%s,%d,%s,%d\n' % (
                    self.sample_name,
                    self.seed,
                    region,
                    self.qcut,
                    left+1,
                    insert_seq,
                    count))

def format_cutoff(cutoff):
    """ Format the cutoff fraction as a string to use as a name. """
    
    if cutoff == MAX_CUTOFF:
        return cutoff
    return '{:0.3f}'.format(cutoff)    

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
    report.write_nuc_variants_header(args.nuc_variants_csv)
    for _key, aligned_reads in groupby(args.aligned_csv,
                                       lambda x: x.split(',')[0:3]):
        report.read(aligned_reads)
        
        report.write_amino_counts(args.amino_csv)
        report.write_consensus(args.conseq_csv)
        report.write_failure(args.failed_align_csv)
        report.write_insertions()
        report.write_nuc_counts(args.nuc_csv)
        report.write_nuc_variants(args.nuc_variants_csv)

    args.aligned_csv.close()
    args.amino_csv.close()
    args.nuc_csv.close()
    args.conseq_csv.close()
    args.coord_ins_csv.close()
    args.failed_align_csv.close()
    args.nuc_variants_csv.close()

if __name__ == '__main__':
    main()

