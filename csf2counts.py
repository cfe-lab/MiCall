#! /usr/bin/env python

"""
Shipyard-style MiSeq pipeline, post-processing 1
Takes aligned CSV as input (produced by sam2csf).
Re-aligns sequence variants to lab standard references (e.g., HXB2).
Reports nucleotide and amino acid frequencies by reference coordinates.
Outputs consensus sequences in same coordinate system.
This assumes a consistent reading frame across the entire region.

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
    parser.add_argument('remap_conseq', help='<input> consensus sequences from remapping step')
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
    Returns map as dictionary object, and a list of non-reference
    positions due to insertions in query.  Note this is not simply
    implied by lack of query position as key in dictionary, because
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

class AminoFrequencyWriter(object):
    def __init__(self, aafile, sample_name, refseqs):
        """ Initialize a writer object.
        
        @param aafile: an open file that the summary will be written to
        @param sample_name: the name of the sample to go in the first column
        @param refseqs: {region: sequence} maps from region name to amino acid 
        sequence
        """
        self.aafile = aafile
        self.sample_name = sample_name
        self.refseqs = refseqs
        self.is_header_written = False
        
    def write(self,
              region,
              qcut,
              qindex_to_refcoord,
              amino_counts,
              inserts):
        """ Write a summary of the amino acid distribution at each position in
        the reference sequence.
        
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
    
        if not self.is_header_written:
            self.aafile.write(
                'sample,region,q-cutoff,query.aa.pos,refseq.aa.pos,'
                'A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,*\n')
            self.is_header_written = True
            
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
            self.aafile.write('%s,%s,%s,%s,%d,%s\n' % (self.sample_name,
                                                       region,
                                                       qcut,
                                                       aa_pos_str,
                                                       refcoord + 1,
                                                       outstr))

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

    # determine reading frames based on region-specific consensus sequences
    # FIXME: this might not be necessary - all offset 0?
    conseqs = {}
    best_frames = {}
    sample_name = ''
    with open(args.remap_conseq, 'rU') as f:
        for i, line in enumerate(f):
            if i == 0:
                #skip header
                continue
            sample_name, region, conseq = line.strip('\n').split(',')
            if region not in refseqs:
                logger.warn('No reference in {} for {} reading {}'.format(
                    amino_ref, 
                    region,
                    args.remap_conseq))
                continue
            refseq = refseqs[region]  # protein sequence

            # determine reading frame based on sample/region consensus
            best_score = -1e5
            for frame in range(3):
                p = translate(conseq, frame)
                #score = align.localds(refseq, p, mat, -20, -5, penalize_end_gaps=False, score_only=True)
                aquery, aref, score = pair_align(hyphy, refseq, p)
                if score > best_score:
                    best_frames[region] = frame
                    best_score = score


    # for each region, quality cutoff
    infile = open(args.aligned_csv, 'rU')
    aafile = open(args.amino_csv, 'w')
    nucfile = open(args.nuc_csv, 'w')
    confile = open(args.conseq, 'w')
    indelfile = open(args.indels_csv, 'w')
    amino_writer = AminoFrequencyWriter(aafile, sample_name, refseqs)

    for region, group in groupby(infile, lambda x: x.split(',')[0]):
        if region not in refseqs:
            continue
        for qcut, group2 in groupby(group, lambda x: x.split(',')[1]):
            # gather nucleotide, amino frequencies
            nuc_counts = {}
            amino_counts = {}
            pcache = {}
            min_offset = 1e6

            total_count = 0
            for line in group2:
                _, _, _, count, offset, seq = line.strip('\n').split(',')
                offset = int(offset)
                count = int(count)
                total_count += count  # track the total number of aligned and merged reads, given QCUT
                seq = seq.upper()

                if offset < min_offset:
                    min_offset = offset

                for i, nuc in enumerate(seq):
                    pos = offset + i
                    if pos not in nuc_counts:
                        nuc_counts.update({pos: {}})
                    if nuc not in nuc_counts[pos]:
                        nuc_counts[pos].update({nuc: 0})
                    nuc_counts[pos][nuc] += count

                p = translate('-'*offset + seq, best_frames[region])
                if p not in pcache:
                    # retain linkage info for reporting insertions
                    pcache.update({p: 0})
                pcache[p] += count

                left = (offset + best_frames[region]) / 3
                right = (offset + len(seq) + best_frames[region]) / 3

                for pos in range(left, right):
                    if pos not in amino_counts:
                        amino_counts.update({pos: {}})
                    aa = p[pos]
                    if aa not in amino_counts[pos]:
                        amino_counts[pos].update({aa: 0})
                    amino_counts[pos][aa] += count

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

            amino_writer.write(region,
                               qcut,
                               qindex_to_refcoord,
                               amino_counts,
                               inserts)

            # output nucleotide frequencies and consensus sequences
            nuc_coords = nuc_counts.keys()
            nuc_coords.sort()
            maxcon = ''
            conseqs = dict([(cut, '') for cut in settings.conseq_mixture_cutoffs])

            for query_nuc_pos in nuc_coords:
                # FIXME: there MUST be a better way to do this :-P
                try:
                    adjustment = best_frames[region] - (3 - min_offset%3)%3
                    query_aa_pos = (query_nuc_pos - min_offset + adjustment) / 3
                    query_codon_pos = (query_nuc_pos - min_offset + adjustment) % 3
                    ref_aa_pos = qindex_to_refcoord[query_aa_pos]
                    ref_nuc_pos = 3*ref_aa_pos + query_codon_pos
                except:
                    logger.warn(
                        'No coordinate mapping for query nuc %d (amino %d) in %s' % (
                            query_nuc_pos,
                            query_aa_pos,
                            args.aligned_csv))
                    continue

                outstr = ','.join(map(str, [nuc_counts[query_nuc_pos].get(nuc, 0) for nuc in 'ACGT']))
                nucfile.write('%s,%s,' % (region, qcut))
                nucfile.write(','.join(map(str, [query_nuc_pos+1, ref_nuc_pos+1, outstr])))
                nucfile.write('\n')

                # append to consensus sequences
                intermed = [(count, nuc) for nuc, count in nuc_counts[query_nuc_pos].iteritems()]
                intermed.sort(reverse=True)
                maxcon += intermed[0][1]  # plurality consensus

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

            confile.write('%s,%s,MAX,%s\n' % (region, qcut, maxcon))
            for cut, conseq in conseqs.iteritems():
                confile.write('%s,%s,%1.3f,%s\n' % (region, qcut, cut, conseq))

            # convert insertion coordinates into contiguous ranges
            if len(inserts) == 0:
                continue

            left = inserts[0]
            last_index = left
            insert_ranges = []
            for i in range(1, len(inserts)):
                index = inserts[i]
                if index - last_index > 1:
                    # skipped an index
                    insert_ranges.append((left, last_index+1))
                    left = index
                last_index = index
            insert_ranges.append((left, last_index+1))

            # enumerate indels by popping out all AA sub-string variants
            indel_counts = {}
            for left, right in insert_ranges:
                indel_counts.update({left: {}})
                for p, count in pcache.iteritems():
                    insert = p[left:right]
                    if insert not in indel_counts[left]:
                        indel_counts[left].update({insert: 0})
                    indel_counts[left][insert] += count

            # record insertions to CSV
            for left in indel_counts.iterkeys():
                for insert in indel_counts[left].iterkeys():
                    indelfile.write('%s,%s,%d,%s,%d\n' % (region, qcut, left, insert, count))

    infile.close()
    aafile.close()
    indelfile.close()

if __name__ == '__main__':
    main()

