import sys
from tcr.igblast import igblast_seq, parse_fmt3
from tcr.idprimers import idp

def load_bins_strict(path):
    fastq = []
    for line in open(path):
        fastq.append(line.strip()) 

    tra, trb, runts = {}, {}, {}
    for _, seq, _, _ in zip(*(iter(fastq),) * 4):
        if len(seq) < 60:
            if seq not in runts:
                runts[seq] = 0
            runts[seq] += 1
            continue

        if len(seq) < 200:
            continue

        if "CCCACCCGAGGTCGCTG" in seq:
            if seq not in trb:
                trb[seq] = 0
            trb[seq] += 1

        if "ACCAGCTGAGAGACTCT" in seq:
            if seq not in tra:
                tra[seq] = 0
            tra[seq] += 1

        if "CCCACCCGAGGTCGCTG" in seq and "ACCAGCTGAGAGACTCT" in seq:
            print("Error: both constant sequences occur in", seq)
            exit(1)
    
    return tra, trb, runts

def load_single_bin(path):
    fastq = []
    for line in open(path):
        fastq.append(line.strip()) 

    bins, runts = {}, {}
    for _, seq, _, _ in zip(*(iter(fastq),) * 4):
        if len(seq) < 60:
            if seq not in runts:
                runts[seq] = 0
            runts[seq] += 1
            continue

        if len(seq) < 200:
            continue

        if seq not in bins:
            bins[seq] = 0
        bins[seq] += 1

    return bins, runts

def id_bins(bins):
    groups = {}
    group_count = {}
    group_max = {}
    cdr3s = {}
    scores = {}
    for k in bins:
        if bins[k] < 500:
            continue
        primer = idp(k)
        seq = k
        if primer is None:
            print("\tUNRECOGNIZED PRIMER", seq[:20])
            continue
        seq = k[len(primer):]
        c = parse_fmt3(igblast_seq(seq))[0]
        print(c)
        calls = c[1]
        cdr3 = c[2]
        genes = []
        for gene in calls: 
            gt = '|'.join(map(lambda x: x[0], calls[gene]))
            genes.append(gt)
        if tuple(genes) not in groups:
            groups[tuple(genes)] = 0
        groups[tuple(genes)] += bins[k]
        cdr3s[tuple(genes)] = cdr3
        scores[tuple(genes)] = calls[gene][0][1]

        if tuple(genes) not in group_count:
            group_count[tuple(genes)] = bins[k]
        if bins[k] >= group_count[tuple(genes)]:
            group_max[tuple(genes)] = k

    return groups, group_max, cdr3s, scores
