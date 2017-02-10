"""
Reimplementation of Conan's pssm_lib.rb Ruby script in Python.
That script originally supported position-specific scoring matrix (PSSM)
and geno2pheno (G2P), but this version only supports G2P.

Based on work published at http://coreceptor.geno2pheno.org
"""

from math import exp
import os

import gotoh

from micall.utils.translation import translate


class Pssm(object):
    def __init__(self, std='g2p', path_to_lookup=None, path_to_matrix=None):
        if std == 'pssm':
            self.std_v3 = 'CTRPNNNTRKGIHIGPGRAFYATGEIIGDIRQAHC'
        elif std == 'nuc':
            self.std_v3 = ('TGTACAAGACCCAACAACAATACAAGAAAAAGTATACATATAGGACCAGGGAG'
                           'AGCATTTTATGCAACAGGAGAAATAATAGGAGATATAAGACAAGCACATTGT')
        elif std == 'g2p':
            self.std_v3 = 'CTRPNXNNTXXRKSIRIXXXGPGQXXXAFYATXXXXGDIIGDIXXRQAHC'
        else:
            raise RuntimeError('Unrecognized argument to std')

        if path_to_lookup is None:
            lookup_paths = [os.path.join(os.path.dirname(__file__), 'g2p_fpr.txt'),
                            'g2p_fpr.txt',
                            '../../g2p/g2p_fpr.txt',
                            'micall/g2p/g2p_fpr.txt']
        else:
            lookup_paths = [path_to_lookup]
        for path in lookup_paths:
            self.g2p_fpr_data = []
            try:
                with open(path, 'rU') as handle:
                    # load empirical curve of FPR to g2p scores from file
                    for line in handle:
                        g2p, fpr = map(float, line.strip('\n').split(','))
                        self.g2p_fpr_data.append((g2p, fpr))
                break
            except:
                self.g2p_fpr_data = None

        if not self.g2p_fpr_data:
            raise RuntimeError('No g2p_fpr data found in {!r}'.format(
                lookup_paths))

        self.g2p_fpr_data.sort()  # make sure the list is sorted

        if path_to_matrix is None:
            matrix_paths = [os.path.join(os.path.dirname(__file__), 'g2p.matrix'),
                            'g2p.matrix',
                            '../../g2p/g2p.matrix',
                            'micall/g2p/g2p.matrix']
        else:
            matrix_paths = [path_to_matrix]
        for path in matrix_paths:
            try:
                with open(path, 'rU') as handle:
                    # load g2p score matrix
                    self.g2p_matrix = dict([(i, {}) for i in range(len(self.std_v3))])
                    for line in handle:
                        items = line.split('\t')
                        residue = items[0]
                        scores = map(float, items[1:])  # V3 reference length
                        for i, score in enumerate(scores):
                            self.g2p_matrix[i].update({residue: score})
                break
            except:
                self.g2p_matrix = None

        if not self.g2p_matrix:
            raise RuntimeError('No g2p matrix data found in {!r}'.format(
                matrix_paths))

    def g2p(self, aa_lists):
        """
        Calculate geno2pheno coreceptor score prediction.
        :param aa_lists: a List of Lists for amino acids per position
        :return:
        """
        rho = 1.33153
        probA = -2.31191
        probB = 0.244784
        ssums = [0.]

        # Calculate sums for all possible sequences, based on ambiguous aminos.
        for i, aa_list in enumerate(aa_lists):
            num_sums = len(ssums)
            if len(aa_list) > 1:
                ssums *= len(aa_list)
            for j, aa in enumerate(aa_list):
                aa_score = self.g2p_matrix[i][aa]
                for k in range(num_sums):
                    ssums[j*num_sums + k] += aa_score

        # Calculate all possible scores, then average.
        score = 0
        for ssum in ssums:
            dv = rho - ssum
            fapb = (dv * probA) + probB
            score += 1. / (1 + exp(fapb-0.5))
        return score/len(ssums)

    def g2p_to_fpr(self, g2p):
        """
        Retrieve FPR value from empirically-derived curve recorded as finite set of values
        in file.  Use bisection search.  In case of inexact match, use midpoint FPR.
        :param g2p:
        :return:
        """
        if g2p is None or g2p < 0.0 or g2p > 1.0:
            return None

        # binary search of sorted list
        left = 0
        right = len(self.g2p_fpr_data)
        while True:
            pivot = (right+left) // 2
            pivot_g2p, pivot_fpr = self.g2p_fpr_data[pivot]
            if g2p == pivot_g2p:
                # found an exact match
                return pivot_fpr
            elif g2p < pivot_g2p:
                right = pivot
            else:
                left = pivot

            if (right-left) == 1:
                # adjacent indices, use one with closest G2P value
                left_g2p, left_fpr = self.g2p_fpr_data[left]
                right_g2p, right_fpr = self.g2p_fpr_data[right]
                if abs(right_g2p - g2p) < abs(left_g2p - g2p):
                    return right_fpr
                return left_fpr

    def align_aminos(self, seq, gapIns=3, removeinserts=False, qachecks=False):
        """
        Align amino acids to a standard reference using gotoh.cpp
        :param seq:  AA sequence in list form, to align against reference standard
        :param removeinserts:  Whether to remove insertions relative to standard
        :param qachecks:  These are not used when [seq] is a list
        :return:
        """
        std = self.std_v3

        if qachecks:
            if seq is None:
                return -1, None
            if (len(seq) % 3 != 0) or len(seq) < 96:
                return -1, None
            if seq.startswith('----') or seq.endswith('----'):
                return -1, None

        if type(seq) is list:
            aa_lists = seq  # aa_seq in pssm_lib.rb
        else:
            # assume this is a codon sequence
            aa_lists = translate(seq=seq, offset=0, resolve=False, return_list=True, ambig_char='X')

        for i, aa_list in enumerate(aa_lists):
            for j, aa in enumerate(aa_list):
                if aa == 'X':
                    aa_list[j] = '-'
            if len(aa_list) > 1 and '*' in aa_list:
                aa_lists[i] = [aa for aa in aa_list if aa != '*']

        while ['-'] in aa_lists:
            aa_lists.remove(['-'])

        # resolve into string
        aa_seq = ''.join(aa_list[0] for aa_list in aa_lists)  # aa_seq_s in pssm_lib.rb

        if qachecks:
            if any(['*' in aa_list for aa_list in aa_lists]):
                return -1, None

        std = std.replace('-', 'X')  # fix gaps in reference
        aligned_std, aligned_seq = gotoh.align_it_aa_rb(std, aa_seq, gapIns, 1)  # method_recall
        aligned_std = aligned_std.replace('X', '-')
        std = std.replace('X', '-')  # restore original state

        # apply alignment to lists
        aa_seq = aligned_seq
        for i in range(len(aa_seq)):
            if aa_seq[i] == '-':
                aa_lists.insert(i, ['-'])  # insert before index, like Ruby

        indels = False
        if removeinserts and '-' in aligned_std:
            new_aa_lists = []
            indices = range(len(aligned_std))
            indices.reverse()
            for i in indices:
                if aligned_std[i] == '-':
                    # skip positions that are insertions relative to standard
                    indels = True
                    continue
                new_aa_lists.append(aa_lists[i])
            aa_lists = new_aa_lists
        else:
            if aligned_std != std:
                # reject sequences with insertions relative to standard
                return -2, None

        return aa_lists, indels

    def run_g2p(self, seqs):
        """
        Wrapper function to calculate g2p score over a set of sequences
        @param seqs: a single sequence (str) or list of sequences
        @return: a tuple of g2p scores (float or list) and the aligned
            sequence (if a single sequence was requested). Score may also be
            None if the alignment failed.
        """
        aa_aligned = ''
        is_array = (type(seqs) is list)
        scores = []
        if not is_array:
            seqs = [seqs]

        for seq in seqs:
            aa, _indels = self.align_aminos(seq, removeinserts=False, qachecks=(type(seq) is not list))
            if not isinstance(aa, list) and aa < 0:
                # failed alignment, try higher gap insert penalty
                aa, _indels = self.align_aminos(seq,
                                                gapIns=6,
                                                removeinserts=False,
                                                qachecks=(type(seq) is not list))  # :recall6
            score = None if not isinstance(aa, list) and aa < 0 else self.g2p(aa)
            scores.append(score)
            if aa == -1:
                aa, _indels = self.align_aminos(seq,
                                                gapIns=6,
                                                removeinserts=False,
                                                qachecks=(type(seq) is not list))  # :recall6
            aa_aligned = aa

        if not is_array:
            return scores[0], aa_aligned
        else:
            return scores, None
