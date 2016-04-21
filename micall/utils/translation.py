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
              translate_mixtures=True,
              list_ambiguous=False,
              stats=None):
    """
    Translate codon (nucleotide) sequence into amino acids.
    :param str seq: the nucleotide sequence
    :param int offset: prefix nucleotide sequence by X bases to shift reading frame
    :param bool resolve: whether to assign an arbitrary residue at ambiguous codons
    :param bool return_list: if True, then returns a list of lists containing all amino
        acid resolutions at each codon
    :param str ambig_char: character to use for ambiguous codons, unless
    translate_mixtures, return_list, or list_ambiguous are true.
    :param bool translate_mixtures: True if a codon should be translated when it has
        an unambiguous mixture, such as CTN translated to L.
    :param bool list_ambiguous: if True, then include all possible amino acids
    in brackets for an ambiguous codon.
    :param dict stats: will have 'length' set to the number of full codons,
    'ambiguous' set to the number of ambiguous positions and max_aminos set to
    the most ambiguous amino acids at any position
    :return: string (AA sequence) or list of lists if return_list=True
    """

    seq = '-'*offset + seq.upper()
    aa_list = []
    aa_seq = ''  # use to align against reference, for resolving indels
    if stats is not None:
        stats['ambiguous'] = 0
        stats['length'] = 0
        stats['max_aminos'] = 1 if seq else 0

    # loop over codon sites in nucleotide sequence
    for codon_site in xrange(0, len(seq), 3):
        codon = seq[codon_site:codon_site+3]

        if len(codon) < 3:
            break

        if stats is not None:
            stats['length'] += 1
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
        elif not translate_mixtures and not list_ambiguous:
            if stats is not None:
                stats['ambiguous'] += 1
            aa_seq += ambig_char
            aa_list.append([ambig_char])
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
                if stats is not None:
                    stats['ambiguous'] += 1
                    stats['max_aminos'] = max(len(aminos), stats['max_aminos'])
                if list_ambiguous or return_list:
                    aminos.sort()
                    aa_seq += '[{}]'.format(''.join(aminos))
                elif resolve:
                    aa_seq += aminos[0]  # arbitrary resolution
                else:
                    aa_seq += ambig_char
            else:
                aa_seq += aminos[0]
    return aa_list if return_list else aa_seq

if __name__ == '__live_coding__':
    import unittest
    from micall.tests.translation_test import TranslateTest

    suite = unittest.TestSuite()
    suite.addTest(TranslateTest("testMixturesNotTranslated2"))
    test_results = unittest.TextTestRunner().run(suite)

    print(test_results.errors)
    print(test_results.failures)
