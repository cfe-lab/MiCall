#! /usr/bin/python

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
import os
import sys
from itertools import groupby
from settings import *
#from Bio.pairwise2 import align  # too slow!

from hyphyAlign import *
hyphy = HyPhy._THyPhy (os.getcwd(), 1)
change_settings(hyphy)

import re

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


parser = argparse.ArgumentParser('Post-processing of short-read alignments.')

parser.add_argument('input_csf', help='<input> aligned CSF input')
parser.add_argument('input_conseq', help='<input> consensus sequences from remapping step')
parser.add_argument('input_amino_ref', help='<input> amino acid reference sequences')
parser.add_argument('output_nuc', help='<output> CSV containing nucleotide frequencies')
parser.add_argument('output_amino', help='<output> CSV containing amino frequencies')
parser.add_argument('output_indels', help='<output> CSV containing insertions')
parser.add_argument('output_conseq', help='<output> CSV containing consensus sequences')

args = parser.parse_args()

# check that the amino acid reference input exists
if not os.path.exists(args.input_amino_ref):
    print 'No input amino reference sequences found at', args.input_amino_ref
    sys.exit(1)

# read in amino acid reference sequences from file
refseqs = {}
with open(args.input_amino_ref, 'rU') as handle:
    header = handle.next()
    for line in handle:
        region, aaseq = line.strip('\n').split(',')
        refseqs.update({region: aaseq})

# check that the inputs exist
if not os.path.exists(args.input_csf):
    print 'No input CSF found at', args.input_csf
    sys.exit(1)

# check that the output paths are valid
for path in [args.output_nuc, args.output_amino, args.output_indels, args.output_conseq]:
    output_path = os.path.split(path)[0]
    if not os.path.exists(output_path) and output_path != '':
        print 'Output path does not exist:', output_path
        sys.exit(1)

# determine reading frames based on region-specific consensus sequences
# FIXME: this might not be necessary - all offset 0?
conseqs = {}
best_frames = {}
with open(args.input_conseq, 'rU') as f:
    for line in f:
        region, conseq = line.strip('\n').split(',')
        if region not in refseqs:
            print 'No reference in', args.input_amino_ref, 'for', region
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
infile = open(args.input_csf, 'rU')
aafile = open(args.output_amino, 'w')
nucfile = open(args.output_nuc, 'w')
confile = open(args.output_conseq, 'w')
indelfile = open(args.output_indels, 'w')

for region, group in groupby(infile, lambda x: x.split(',')[0]):
    if region not in refseqs:
        continue
    for qcut, group2 in groupby(group, lambda x: x.split(',')[1]):
        # gather nucleotide, amino frequencies
        nuc_counts = {}
        amino_counts = {}
        pcache = {}
        min_offset = 1e6

        for line in group2:
            region2, qcut2, rank, count, offset, seq = line.strip('\n').split(',')
            offset = int(offset)
            count = int(count)
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
        aquery, aref, ascore = pair_align(hyphy, refseqs[region], aa_max)
        left, right = get_boundaries(aref)	# Coords of first/last non-gap character

        qindex_to_refcoord = {}
        inserts = []
        qindex = 0
        rindex = 0

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

        # output amino acid frequencies
        for qindex, refcoord in qindex_to_refcoord.iteritems():
            aa_pos = qindex + min(aa_coords)
            if aa_pos in inserts:
                continue
            outstr = ','.join(map(str, [amino_counts[aa_pos].get(aa, 0) for aa in amino_alphabet]))
            aafile.write('%s,%s,%d,%d,%s\n' % (region, qcut, aa_pos, refcoord+1, outstr))

        # output nucleotide frequencies
        nuc_coords = nuc_counts.keys()
        nuc_coords.sort()
        maxcon = ''
        conseqs = dict([(cut, '') for cut in conseq_mixture_cutoffs])

        for query_nuc_pos in nuc_coords:
            # FIXME: there MUST be a better way to do this :-P
            try:
                adjustment = best_frames[region] - (3 - min_offset%3)%3
                query_aa_pos = (query_nuc_pos - min_offset + adjustment) / 3
                query_codon_pos = (query_nuc_pos - min_offset + adjustment) % 3
                ref_aa_pos = qindex_to_refcoord[query_aa_pos]
                ref_nuc_pos = 3*ref_aa_pos + query_codon_pos
            except:
                print 'No coordinate mapping for query nuc %d (amino %d)' % (query_nuc_pos, query_aa_pos)
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
            for cut in conseq_mixture_cutoffs:
                mixture = []
                for count, nuc in intermed:
                    if float(count) / total_count > cut:
                        mixture.append(nuc)

                if 'N' in mixture:
                    if len(mixture) > 1:
                        mixture.remove('N')
                    else:
                        conseqs[cut] += 'N'
                        continue

                if '-' in mixture:
                    if len(mixture) > 1:
                        mixture.remove('-')
                    else:
                        conseqs[cut] += '-'
                        continue

                if len(mixture) > 1:
                    mixture.sort()
                    conseqs[cut] += ambig_dict[''.join(mixture)]
                elif len(mixture) == 1:
                    conseqs[cut] += mixture[0]
                else:
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

        # enumerate insertions
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



