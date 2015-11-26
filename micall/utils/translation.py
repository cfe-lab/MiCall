"""
Utility function for translating nucleotides (codons) into amino acids.
"""
import re

codon_dict = {'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
              'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
              'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
              'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
              'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
              'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
              'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
              'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
              'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
              'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
              'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
              'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
              'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
              'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
              'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
              'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
              '---': '-', 'XXX': '?'}
mixture_regex = re.compile('[WRKYSMBDHVN-]')
mixture_dict = {'W': 'AT', 'R': 'AG', 'K': 'GT', 'Y': 'CT', 'S': 'CG',
                'M': 'AC', 'V': 'AGC', 'H': 'ATC', 'D': 'ATG',
                'B': 'TGC', 'N': 'ATGC', '-': 'ATGC'}
ambig_dict = dict(("".join(sorted(v)), k)
                  for k, v in mixture_dict.iteritems()
                  if k != '-')
complement_dict = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A',
                   'W': 'S', 'R': 'Y', 'K': 'M', 'Y': 'R', 'S': 'W', 'M': 'K',
                   'B': 'V', 'D': 'H', 'H': 'D', 'V': 'B',
                   '*': '*', 'N': 'N', '-': '-'}


def reverse_and_complement(seq):
    return ''.join(complement_dict[nuc] for nuc in reversed(seq))


def translate(seq,
              offset=0,
              resolve=False,
              return_list=False,
              ambig_char='?',
              translate_mixtures=True):
    """
    Translate codon (nucleotide) sequence into amino acids.
    @param seq: the nucleotide sequence
    @param offset: prefix nucleotide sequence by X bases to shift reading frame
    @param resolve: whether to assign an arbitrary residue at ambiguous codons
    @param return_list: if True, then returns a list of lists containing all amino
        acid resolutions at each codon
    @param translate_mixtures: True if a codon should be translated when it has
        an unambiguous mixture, such as CTN translated to L.
    @return: string (AA sequence) or list of lists if return_list=True
    """

    seq = '-'*offset + seq.upper()
    aa_list = []
    aa_seq = ''  # use to align against reference, for resolving indels

    # loop over codon sites in nucleotide sequence
    for codon_site in xrange(0, len(seq), 3):
        codon = seq[codon_site:codon_site+3]

        if len(codon) < 3:
            break

        # note that we're willing to handle a single missing nucleotide as an ambiguity
        if codon.count('-') > 1 or '?' in codon:
            if codon == '---':  # don't bother to translate incomplete codons
                aa_seq += '-'
                aa_list.append(['-'])
            else:
                aa_seq += ambig_char
                aa_list.append([ambig_char])
            continue

        # look for nucleotide mixtures in codon, resolve to alternative codons if found
        num_mixtures = len(mixture_regex.findall(codon))

        if num_mixtures == 0:
            aa = codon_dict[codon]
            aa_seq += aa
            aa_list.append([aa])
        elif not translate_mixtures:
            aa_seq += ambig_char
        else:
            # expand codon into all possible resolutions of mixtures
            codons = [codon]
            while True:
                next_codons = []
                for codon in codons:
                    mixtures = mixture_regex.findall(codon)
                    if len(mixtures) == 0:
                        next_codons.append(codon)
                        continue
                    pos = codon.index(mixtures[0])
                    for nuc in mixture_dict[mixtures[0]]:
                        next_codons.append(codon[0:pos] + nuc + codon[(pos+1):])

                if len(codons) == len(next_codons):
                    # no change in number of codons, exit
                    break
                codons = next_codons

            aminos = list(set([codon_dict[codon] for codon in codons]))
            aa_list.append(aminos)

            if len(aminos) > 1:
                if resolve:
                    aa_seq += aminos[0]  # arbitrary resolution
                else:
                    aa_seq += ambig_char
            else:
                aa_seq += aminos[0]

    return aa_list if return_list else aa_seq
