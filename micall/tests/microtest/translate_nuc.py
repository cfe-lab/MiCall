from micall.core import aln2counts

""" Translate nucleotide sequence to amino acid sequence and compare the result
with an expected sequence.
You might also want to identify the amino acid sequence with BLASTp:
http://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastp&PAGE_TYPE=BlastSearch
"""
from csv import DictReader
import re
from collections import Counter

nuc_seq = ''.join([
        "CCTCAGGTCACTCTTTGGCAACGACCCCTCGTCACAATAAAGATAGGGGGGCAACTAAAGGAAGC",
        "TCTATTAGATACAGGAGCAGATGATACAGTATTAGAAGAAATGAGTTTGCCAGGAAGATGGAAAC",
        "CAAAAATGATAGGGGGAATTGGAGGTTTTATCAAAGTAAGACAGTATGATCAGATACTCATAGAA",
        "ATCTGTGGACATAAAGCTATAGGTACAGTATTAGTAGGACCTACACCTGTCAACATAATTGGAAG",
        "AAATCTGTTGACTCAGATTGGTTGCACTTTAAATTTT"
      ])
aa_compare = ''.join([
        "PQVTLWQRPLVTIKIGGQLKEALLDTGADDTVLEEMSLPGRWKPKMIGGIGGFIKVRQYDQILIE",
        "ICGHKAIGTVLVGPTPVNIIGRNLLTQIGCTLNF"
      ])
if False:
    aa_seq = aln2counts.translate(nuc_seq)
    pairs = zip(aa_seq, aa_compare)
    diffs = [' ' if a == b else '*' for a, b in pairs]
    print 'result ', aa_seq
    print 'diffs  ', ''.join(diffs) if aa_seq != aa_compare else 'no diffs'
    print 'compare', aa_compare

if False:
    fname = '/home/don/data/miseq/150529_M01841/63201A-V3-3-V3LOOP_S78_rerun.failed_read.csv'
    cmrpn_count = 0
    ctrhn_count = 0
    with open(fname, 'rU') as f:
        reader = DictReader(f)
        for i, row in enumerate(reader):
            if False and i >= 10:
                break
            for seq in (row['seq1'], row['seq2']):
                for reading_frame in range(3):
                    prot = aln2counts.translate('-'*reading_frame + seq.strip())
                    m = re.search('(CMRP.*RRAHC|CTRHN.*RKADC)', prot)
                    if m:
                        if m.group().startswith('CMRPN'):
                            cmrpn_count += 1
                        else:
                            ctrhn_count += 1
    print cmrpn_count, ctrhn_count
    
if False:
    fname = '/home/don/data/miseq/150529_M01841/63201A-V3-3-V3LOOP_S78_L001_R1_001.fastq'
    counts = Counter()
    
    print 'Starting.'
    with open(fname, 'rU') as f:
        for i, (header, seq, _plus, _qual) in enumerate(zip(f, f, f, f)):
            if False and i >= 10:
                break
            match_text = ''
            for reading_frame in range(3):
                prot = aln2counts.translate('-'*reading_frame + seq.strip())
                m = re.search('(CMRP.*AHC|CTRHN.*RKADC)', prot)
                if m:
                    counts[m.group()] += 1
    for prot, count in counts.most_common(100):
        print count, prot

if False:
    fname = '/home/don/data/miseq/150529_M01841/63221A-V3-1-V3LOOP_S19.unmapped1.fastq'
    count = 0
    
    print 'Starting.'
    with open(fname, 'rU') as f:
        for i, (header, seq, _plus, _qual) in enumerate(zip(f, f, f, f)):
            if False and i >= 10:
                break
            match_text = ''
            for reading_frame in range(3):
                prot = aln2counts.translate('-'*reading_frame + seq.strip())
                if 'CARPNNNTRKSIRIHRGPGRAFFTSGIKGDIRQAHC' in prot:
                    count += 1
    print count

if False:
    fname = '/home/don/data/miseq/150529_M01841/63221A-V3-1-V3LOOP_S19.remap.csv'
    
    print 'Starting remap.'
    count = 0
    with open(fname, 'rU') as f:
        reader = DictReader(f)
        for i, row in enumerate(reader):
            if False and i >= 10:
                break
            for reading_frame in range(3):
                prot = aln2counts.translate('-'*reading_frame + row['seq'])
                if 'CARPNNNTRKSIRIHRGPGRAFFTSGIKGDIRQAHC' in prot:
                    print row
                    exit()
                    count += 1
    print count

if True:
    fname = '/home/don/data/miseq/150529_M01841/63221A-V3-1-V3LOOP_S19_rerun.aligned.csv'
    
    print 'Starting aligned and stripped.'
    count = 0
    with open(fname, 'rU') as f:
        reader = DictReader(f)
        for i, row in enumerate(reader):
            if False and i >= 10:
                break
            seq = row['seq'].replace('-', '')
            for reading_frame in range(3):
                prot = aln2counts.translate('-'*reading_frame + seq)
                if 'CARPNNNTRKSIHRGPGRAFFTSGIKGDIRQAHC' in prot:
                    count += int(row['count'])
    print count

if False:
    seq0 = 'AATGCCAAAACCATAATAGTACAGCTGGACAAACCTGTAGAAATTAATTGTGCAAGACCCAACAACAATACAAGAAAAAGTATACATCGAGGACCAGGAAGAGCATTTTTTACATC---AGGCATAAAAGGAGATATAAGACAAGCACACTGTAACATTAGTGGAACAGAATGGAATAACACTTTACAAAAGGTAGCAGAAAAATTAAGAGAACAATTTAAAAGTGAAAATATAACCTTTAGGCCATCCTCAGGAGGAGACCCAGAAATTGTCATGCACAATTTTAATTGTGGAAGGGAATTTTTC'
    seq1 = 'AATGCCAAAACCATAATAGTACAGCTGGACAAACCTGTAGAAATTAATTGTGCAAGACCCAACAACAATACAAGAAAAAGTATACATCGAGGACCAGGAAGAGCATTTTTTACATC---AGGCATAAAAGGAGATATAAGACAAGCACACTGTAACATTAGTGGAACAGAATGGAATAACACTTTACAAAAGGTAGCAGAAAAATTAAGAGAACAATTTAAAAGTGAAAATATAACCTTCAGGCCATCCTCAGGAGGAGACCCAGAAATTGTCATGCACAATTTTAATTGTGGAAGGGAATTTTTC'
    seq2 = 'AATGCCAAAACCATAATAGTACAGCTGGACAAACCTGTAGAAATTAATTGTGCAAGACCCAACAACAATACAAGAAAAAGTATACATCGAGGACCAGGAAGAGCATTTTTTACATC---AGGCATAAAAGGAGATATAAGACAAGCACACTGTAACATTAGTGGAACAGAATGGAATAACACTTTACAAAAGGTAGCAGAAAAATTAAGAGAACAATTTAAAAGTAAAAATATAACCTTTAGGCCATCCTCAGGAGGAGACCCAGAAATTGTCATGCACAATTTTAATTGTGGAAGGGAATTTTTC'
    count = 0
    for reading_frame in range(3):
        pass
        prot0 = aln2counts.translate('-'*reading_frame + seq0)
        prot1 = aln2counts.translate('-'*reading_frame + seq1)
        prot2 = aln2counts.translate('-'*reading_frame + seq2)
        print prot0
        print prot1
        print prot2
        if 'CARPNNNTRKSIRIHRGPGRAFFTSGIKGDIRQAHC' in prot0:
            count += 1
"""
CARPN-NNT--RKSIRI---HRGPGR-AFFTS-----GIKGDI--RQAHC
CARPNNNTRKSIRIHRGPGRAFFTSGIKGDIRQAHC
CARPNNNTRKSIRIHRGPGRAFFTSGIKGDIRQAHC (14423)

prot0 = 'NAKTIIVQLDKPVEINCARPNNNTRKSIRIHRGPGRAFFTSGIKGDIRQAHCNISGTEWNNTLQKVAEKLREQFKSENITFRP' 
prot1 = 'NAKTIIVQLDKPVEINCARPNNNTRKSIHRGPGRAFFTS?GIKGDIRQAHCNISGTEWNNTLQKVAEKLREQFKSENITFRP' 

after alignment:
CARPNNNTRKSIRIHRGPGRAFFTS-GIKGDIRQAHC target
CARPNNNTRKSI--HRGPGRAFFTS?GIKGDIRQAHC 2608
CARPNNNTRKSI--HRGPGRAFFTS?GIKGDIRQAHC 2175
CARPNNNTRKSI--HRGPGRAFFTS?GIKGDIRQAHC 626
CARPNNNTRKSIHRGPGRAFFTS?GIKGDIRQAHC
mseq counts from sam2aln:
4102 failed
2608 ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------AATGCCAAAACCATAATAGTACAGCTGGACAAACCTGTAGAAATTAATTGTGCAAGACCCAACAACAATACAAGAAAAAGTATACATCGAGGACCAGGAAGAGCATTTTTTACATC---AGGCATAAAAGGAGATATAAGACAAGCACACTGTAACATTAGTGGAACAGAATGGAATAACACTTTACAAAAGGTAGCAGAAAAATTAAGAGAACAATTTAAAAGTGAAAATATAACCTTTAGGCCATCCTCAGGAGGAGACCCAGAAATTGTCATGCACAATTTTAATTGTGGAAGGGAATTTTTC
2175 ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------AATGCCAAAACCATAATAGTACAGCTGGACAAACCTGTAGAAATTAATTGTGCAAGACCCAACAACAATACAAGAAAAAGTATACATCGAGGACCAGGAAGAGCATTTTTTACATC---AGGCATAAAAGGAGATATAAGACAAGCACACTGTAACATTAGTGGAACAGAATGGAATAACACTTTACAAAAGGTAGCAGAAAAATTAAGAGAACAATTTAAAAGTGAAAATATAACCTTCAGGCCATCCTCAGGAGGAGACCCAGAAATTGTCATGCACAATTTTAATTGTGGAAGGGAATTTTTC
626  ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------AATGCCAAAACCATAATAGTACAGCTGGACAAACCTGTAGAAATTAATTGTGCAAGACCCAACAACAATACAAGAAAAAGTATACATCGAGGACCAGGAAGAGCATTTTTTACATC---AGGCATAAAAGGAGATATAAGACAAGCACACTGTAACATTAGTGGAACAGAATGGAATAACACTTTACAAAAGGTAGCAGAAAAATTAAGAGAACAATTTAAAAGTAAAAATATAACCTTTAGGCCATCCTCAGGAGGAGACCCAGAAATTGTCATGCACAATTTTAATTGTGGAAGGGAATTTTTC
201  ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------AATGCCAAAACCATAATAGTACAGCTGGATAAACCTGTAGAAATTAATTGTGCAAGACCCAACAACAATACAAGAAAAAGTATACATCGAGGACCAGGAAGAGCATTTTTTACATC---AGGCATAAAAGGAGATATAAGACAAGCACACTGTAACATTAGTGGAACAGAATGGAATAACACTTTACAAAAGGTAGCAGAAAAATTAAGAGAACAATTTAAAAGTGAAAATATAACCTTTAGGCCATCCTCAGGAGGAGACCCAGAAATTGTCATGCACAATTTTAATTGTGGAAGGGAATTTTTC
181  ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------AATGCCAAAACCATAATAGTACAGCTGGACAAACCTGTAGAAATTAATTGTGCAAGACCCAACAACAATACAAGAAAAAGTATACATCGAGGACCAGGAAGAGCATTTTTTACATC---AGGCATAAAAGGAGATATAAGACAAGCACACTGTAACATTAGTGGAACAGAATGGAATAACACTTTACAAAAGGTAGCAGAAAAATTAAGAGAACAATTTAAAAGTAGAAATATAACCTTTCAGCCATCCTCAGGAGGAGACCCAGAAATTGTCATGCACAGTTTTAATTGTGGAAGGGAATTTTTC
126  ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------AATGCCAAAACCATAATAGTACAGCTGGATAAACCTGTAGAAATTAATTGTGCAAGACCCAACAACAATACAAGAAAAAGTATACATCGAGGACCAGGAAGAGCATTTTTTACATC---AGGCATAAAAGGAGATATAAGACAAGCACACTGTAACATTAGTGGAACAGAATGGAATAACACTTTACAAAAGGTAGCAGAAAAATTAAGAGAACAATTTAAAAGTGAAAATATAACCTTCAGGCCATCCTCAGGAGGAGACCCAGAAATTGTCATGCACAATTTTAATTGTGGAAGGGAATTTTTC
102  ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------AATGCCAAAACCATAATAGTACAGCCGGACAAACCTGTAGAAATTAATTGTGCAAGACCCAACAACAATACAAGAAAAAGTATACATCGAGGACCAGGAAGAGCATTTTTTACATC---AGGCATAAAAGGAGATATAAGACAAGCACACTGTAACATTAGTGGAACAGAATGGAATAACACTTTACAAAAGGTAGCAGAAAAATTAAGAGAACAATTTAAAAGTGAAAATATAACCTTCAGGCCATCCTCAGGAGGAGACCCAGAAATTGTCATGCACAATTTTAATTGTGGAAGGGAATTTTTC
74   ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------AATGCCAAAACCATAATAGTACAGCTGGATAAACCTGTAGAAATTAATTGTGCAAGACCCAACAACAATACAAGAAAAAGTATACATCGAGGACCAGGAAGAGCATTTTTTACATC---AGGCATAAAAGGAGATATAAGACAAGCACACTGTAACATTAGTGGAACAGAATGGAATAACACTTTACAAAAGGTAGCAGAAAAATTAAGAGAACAATTTAAAAGTAAAAATATAACCTTTAGGCCATCCTCAGGAGGAGACCCAGAAATTGTCATGCACAATTTTAATTGTGGAAGGGAATTTTTC
62   ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------AATGCCAAAACCATAATAGTACAGCTGGACAAACCTGTAGAAATTAATTGTGCAAGACCCAACAACAATACAAGAAAAAGTATACATCGAGGACCAGGAAGAGCATTTTTTACATC---AGGCATAAAAGGAGATATAAGACAAGCACACTGTAACATTAGTGGAACAGAATGGAATAACACTTTACAAAAGGTAGCAGAAAAACTAAGAGAACAGTTTAAAAGTAAAAATATAACCTTTAGGCCATCCTCAGGAGGAGACCCAGAAATTGTCATGCACAATTTTAATTGTGGAAGGGAATTTTTC
50 ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------AATGCCAAAACCATAATAGTACAGCTGGACAAACCTGTAGAAATTAATTGTGCAAGACCCAACAACAATACAAGAAAAAGTATACATCGAGGACCAGGAAGAGCATTTTTTACATC---AGGCATAAAAGGAGATATAAGACAAGCACACTGTAACATTAGTGGAACAGAATGGAATAACACTTTACAAAAGGTAGCAGAAAAATTAAGAGAACAATTTAAAAGTGAAAATATAACCTTTAGGCCATCCTCAGGAGGAGACCCAGAAATTGTCATGCACAGTTTTAATTGTGGAAGGGAATTTTTC
36 ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------AATGCCAAAACCATAATAGTACAGCTGGACAAACCTGTAGAAATTAATTGTGCAAGACCCAACAACAATACAAGAAAAAGTATACATCGAGGACCAGGAAGAGCATTTTTTACATC---AGGCATAAAAGGAGATATAAGACAAGCACACTGTAACATTAGTGGAGCAGAATGGAATAACACTTTACAAAAGGTAGCAGAAAAATTAAGAGAACAATTTAAAAGTAAAAATATAACCTTTAGGCCATCCTCAGGAGGAGACCCAGAAATTGTCATGCACAATTTTAATTGTGGAAGGGAATTTTTC
34 ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------AATGCCAAAACCATAATAGTACAGCTGGACAAACCTGTAGAAATTAATTGTGCAAGACCCAACAACAATACAAGAAAAAGTATACATCGAGGACCAGGAAGAGCATTTTTTACATC---AGGCATAAAAGGAGATATAAGACAAGCACACTGTAACATTAATGGAACAGAATGGAATGACACTTTACAAAAGGTAGCTGAAAAATTAAGAGAACAATTTAAAGTTAAAAATATAACCTTTAGGCCATCCTCAGGAGGAGACCCAGAAATTGTCATGCACAATTTTAATTGTGGAAGGGAATTTTTC
27 ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------AATGCCAAAACCATAATAGTACAGCTGGACAAACCTGTAGAAATTAATTGTGCAAGACCCAACAACAATACAAGAAAAAGTATACATCGAGGACCAGGAAGAGCATTTTTTACATC---AGGCATAAAAGGAGATATAAGACAAGCACACTGTAACATTAGTGGAACAGAATGGAATAACACTTTACAAAAGGTAGCAGAAAAATTAAGAGAACAATTTAAAAGTGAAAATATAACCTTTAGGCCATCCTCAGGAGGAGACCCAGAAATTATCATGCACAATTTTAATTGTGGAAGGGAATTTTTC
25 ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------AATGCCAAAACCATAATAGTACAGCTGGACAAACCTGTAGAAATTAATTGTGCAAGACCCAACAACAATACAAGAAAAAGTATACATCGAGGACCAGGAAGAGCATTTTTTACATC---AGGCATAAAAGGAGATATAAGACAAGCACACTGTAACATTAGTGGAACAGAATGGAATAACACTTTACAAAAGGTAGCAGAAAAATTAAGAGAACAATTTAAAAGTAGAAATATAACCTTTCAGCCATCCTCAGGAGGAGACCCAGAAATTGTCATGCACAATTTTAATTGTGGAAGGGAATTTTTC
23 ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------AATGCCAAAACCATAATAGTACAGCCGGACAAACCTGTAGAAATTAATTGTGCAAGACCCAACAACAATACAAGAAAAAGTATACATCGAGGACCAGGAAGAGCATTTTTTACATC---AGGCATAAAAGGAGATATAAGACAAGCACACTGTAACATTAGTGGAACAGAATGGAATAACACTTTACAAAAGGTAGCAGAAAAATTAAGAGAACAATTTAAAAGTGAAAATATAACCTTTAGGCCATCCTCAGGAGGAGACCCAGAAATTGTCATGCACAATTTTAATTGTGGAAGGGAATTTTTC
23 ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------AATGCCAAAACCATAATAGTACAGCTGGACAAACCTGTAGAAATTAATTGTGCAAGACCCAACAACAATACAAGAAAAAGTATACATCGAGGACCAGGAAGAGCATTTTTTACATC---AGGCATAAAAGGAGATATAAGACAAGCACACTGTAACATTAGTGGAACAGAATGGAATAACACTTTACAAAAGGTAGCAGAAAAATTAAGAGAACAATTTAAAAGTGAAAATATAACCTTTAGGCCATCCTCAGGAGGAGACCCAGAAATTGTCATGCATAATTTTAATTGTGGAAGGGAATTTTTC
21 ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------AATGCCAAAACCATAATAGTACAGCTGGACAAACCTGTAGAAATTAATTGTGCAAGACCCAACAACAATACAAGAAAAAGTATACGTCGAGGACCAGGAAGAGCATTTTTTACATC---AGGCATAAAAGGAGATATAAGACAAGCACACTGTAACATTAGTGGAACAGAATGGAATAACACTTTACAAAAGGTAGCAGAAAAATTAAGAGAACAATTTAAAAGTGAAAATATAACCTTTAGGCCATCCTCAGGAGGAGACCCAGAAATTGTCATGCACAATTTTAATTGTGGAAGGGAATTTTTC
21 ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------AATGCCAAAACCATAATAGTACAGCTGGACAAACCTGTAGAAATTAATTGTGCAAGACCCAACAACAATACAAGAAAAAGTATACATCGAGGACCAGGAAGAGCATTTTTTACATC---AGGCATAAAAGGAGATATAAGACAAGCACACTGTAACATTAGTGGAACAAAATGGAATAACACTTTACAAAAGGTAGCAGAAAAATTAAGAGAACAATTTAAAAGTAAAAATATAACCTTTAGGCCATCCTCAGGAGGAGACCCAGAAATTGTCATGCACAATTTTAATTGTGGAAGGGAATTTTTC
21 ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------AATGCCAAAACCATAATAGTACAGCTGGACAAACCTGTAGAAATTAATTGTGCAAGACCCAACAACAATACAAGAAAAAGTATACATCGAGGACCAGGAAGAGCATTTTTTACATC---AGGCATAAAAGGAGATATAAGACAAGCACACTGTAACATTAGTGGAACAGAATGGAATAACACTTTACAAAAGGTAGCAGAAAAATTAAGAGAACAATTTAAAGTTAAAAATATAACCTTTAGGCCATCCTCAGGAGGAGACCCAGAAATTGTCATGCACAATTTTAATTGTGGAAGGGAATTTTTC
21 ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------AATGCCAAAACCATAATAGTACAGCTGGATAAACCTGTAGAAATTAATTGTGCAAGACCCAACAACAATACAAGAAAAAGTATACATCGAGGACCAGGAAGAGCATTTTTTACATC---AGGCATAAAAGGAGATATAAGACAAGCACACTGTAACATTAGTGGAACAGAATGGAATAACACTTTACAAAAGGTAGCAGAAAAATTAAGAGAACAATTTAAAAGTAGAAATATAACCTTTCAGCCATCCTCAGGAGGAGACCCAGAAATTGTCATGCACAGTTTTAATTGTGGAAGGGAATTTTTC
19 ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------AATGCCAAAACCATAATAGTACAGCTGGACAAACCTGTAGAAATTAATTGTGCAAGACCCAACAACAATACAAGAAAAAGTATACATCGAGGACCAGGAAGAGCATTTTTTACATC---AGGCATAAAAGGAGATATAAGACAAGCACACTGTAACATTAGTGGAACAGAATGGAATAACACTTTACAAAAGGTAGCAGAAAAATTAAGAGAACAATTTAAAAGTGAAAATATAACCTTTAGGCCATCCTCAGGAGGAGACCCAGAAATTGCCATGCACAATTTTAATTGTGGAAGGGAATTTTTC
19 ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------AATGCCAAAACCATAATAGTACAGCTGGACAAACCTGTAGAAATTAATTGTGCAAGACCCAACAACAATACAAGAAAAAGTATACATCGAGGACCAGGAAGAGCATTTTTTACATC---AGGCATAAAAGGAGATATAAGACAAGCACACTGTAACATTAGTGGAACAGAATGGAATAACACTTTACAAAAGGTAGCAGAAAAATTAAGAGAACAATTTAAAGTTAAAAATATAACCTTTAAGCCATCCTCAGGAGGAGACCCAGAAATTGTCATGCATAATTTTAATTGTGGAAGGGAATTTTTC
17 ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------AATGCCAAAACCATAATAGTACAGCTGGACAAACCTGTAGAAATTAATTGTGCAAGACCCAACAACAATACAAGAAAAAGTATACATCGAGGACCAGGAAGAGCATTTTTTACATC---AGGCATAAAAGGAGATATAAGACAAGCACACTGTAACATTAGTGGAACAGAATGGAATAACACTTTACAAAAGGTAGCTGAAGAATTAAGAGAACAATTTAAAGTTAAAAATATAACCTTTAAGCCATCCTCAGGAGGAGACCCAGAAATTGTCATGCACAATTTTAATTGTGGAAGGGAATTTTTC
16 ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------AATGCCAAAACCATAATAGTACAGCTGGACAAACCTGTAGAAATTAATTGTGCAAGACCCAACAACAATACAAGAAAAAGTATACATCGAGGACCAGGAAGAGCATTTTTTACATC---AGGCATAAAAGGAGATATAAGACAAGCACACTGTAACATTAGTGGAACAGAATGGAATAACACTTTACAAAAGGTAGCAGAAAAATTAAGAGAACAATTTAAAAGTAAAAATATAACCTTTAGGCCATCCTCAGGAGGAGACCCAGAAATTGTCATGCACAGTTTTAATTGTGGAAGGGAATTTTTC
16 ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------AATGCCAAAACCATAATAGTACAGCTGGACAAACCTGTAGAAATTAATTGTGCAAGACCCAACAACAATACAAGAAAAAGTATACATCGAGGACCAGGAAGAGCATTTTTTACATC---AGGCATAAAAGGAGATATAAGACAAGCACACTGTAACATTAGTGGAACAGAATGGAATAACACTTTACAAAAGGTAGCTGAAAAATTAAGAGAACAATTTAAAGTTAAAAATATAACCTTTAGGCCATCCTCAGGAGGAGACCCAGAAATTGTCATGCACAATTTTAATTGTGGAAGGGAATTTTTC
16 ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------AATGCCAAAACCATAATAGTACAGCTGGACAAACCTGTAGAAATTAATTGTGCAAGACCCAACAACAATACAAGAAAAAGTATACATCGAGGACCAGGAAGAGCATTTTTTACATC---AGGCATAAAAGGAGATATAAGACAAGCACACTGTAACATTAGTGGAACAGAATGGAATAACACTTTACAAAAGGTAGCAGAAAAATTAAGAGAACAATTTAAAAGTAAAAATATAACCTTCAGGCCATCCTCAGGAGGAGACCCAGAAATTGTCATGCACAATTTTAATTGTGGAAGGGAATTTTTC
16 ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------AATGCCAAAACCATAATAGTACAGCTGGACAAACCTGTAGAAATTAATTGTGCAAGACCCAACAACAATACAAGAAAAAGTATACATCGAGGACCAGGAAGAGCATTTTTTACATC---AGGCATAAAAGGAGATATAAGACAAGCACACTGTAACATTAGTGGAACAGAATGGAATAACACTTTACAAAAGGTAGCAGAAAAATTAAGAGAACAATTTAAAAGTGAAAATATAACCTTCAGGCCATCCTCAGGAGGAGACCCAGAAATTGTCATGCACAGTTTTAATTGTGGAAGGGAATTTTTC
15 ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------AATGCCAAAACCATAATAGTACAGCTGGACAAACCTGTAGAAATTAATTGTGCAAGACCCAACAACAATACAAGAAAAAGTATACATCGAGGACCAGGAAGAGCATTTTTTACATC---AGGCATAAAAGGAGATATAAGACAAGCACACTGTAACATTAGTGGAACAGAATGGAATAACACTTTACAAAAGGTAGCTGAAAAATTAAGAGAACAATTTAAAGTTAAAAATATAACCTTTAAGCCATCCTCAGGAGGAGACCCAGAAATTGTCATGCATAATTTTAATTGTGGAAGGGAATTTTTC
14 ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------AATGCCAAAACCATAATAGTACAGCCGGACAAACCTGTAGAAATTAATTGTGCAAGACCCAACAACAATACAAGAAAAAGTATACATCGAGGACCAGGAAGAGCATTTTTTACATC---AGGCATAAAAGGAGATATAAGACAAGCACACTGTAACATTAGTGGAACAGAATGGAATAACACTTTACAAAAGGTAGCAGAAAAATTAAGAGAACAATTTAAAAGTAAAAATATAACCTTTAGGCCATCCTCAGGAGGAGACCCAGAAATTGTCATGCACAATTTTAATTGTGGAAGGGAATTTTTC
total mseq counts: 14238

prot counts from sam_g2p
prot counts:
2607 CARPNNNTRKSIXXXXGXXXXXXTXIXXXXXXXXXHC
787 CARPNNNTRKSIRXXXGXXXXXXTXIXXXXXXXXXHC
395 CARPNNNTRKSIRIXXGXXXXXXTXIXXXXXXXXXHC
259 CARPNNNTRKSIHRGPGRAFFTSGIKGDIRQAHC
244 CARPNNNTRKSIXXXXGXGXXXXTXIXXXXXXXXXHC
174 CARPNNNTRKSIXXXPGXXXXXXTXIXXXXXXXXXHC
156 CARPNNNTRKSIXXXXGXXXXXXTXIXXXXXXXQXHC
126 CARPNNNTRKSIXIXXGXXXXXXTXIXXXXXXXXXHC
94 CARPNNNTRKSIRXXXGXGXXXXTXIXXXXXXXXXHC
84 CARPNNNTRKSIXXXRGXXXXXXTXIXXXXXXXXXHC
76 CARPNNNTRKSIRXXPGXXXXXXTXIXXXXXXXXXHC
65 CARPNNNTRKSIXXXXGXXXXXXTSIXXXXXXXXXHC
63 CARPNNNTRKSIRXXXGXXXXXXTXIXXXXXXXQXHC
53 CARPNNNTRKSIXXXXGXXXXXFTXIXXXXXXXXXHC
51 CARPNNNTRKSIRIXXGXXXXXXTXIXXXXXXXQXHC
51 CARPNNNTRKSIRIXXGXGXXXXTXIXXXXXXXXXHC
46 CARPNNNTRKSIRXXRGXXXXXXTXIXXXXXXXXXHC
40 CARPNNNTRKSIXXXXGXAXXXXTXIXXXXXXXXXHC
36 CARPNNNTRKSIXXXXGXXXXXXTXIXXXXXQXXXHC
35 CARPNNNTRKSIRIXPGXXXXXXTXIXXXXXXXXXHC
34 CARPNNNTRKSIXXXXGXXXAXXTXIXXXXXXXXXHC
32 CARPNNNTRKSIXXXXGPXXXXXTXIXXXXXXXXXHC
32 CARPNNNTRKSIXXXXGXXRXXXTXIXXXXXXXXXHC
30 CARPNNNTRKSIRXXXGXXXXXXTSIXXXXXXXXXHC
27 CARPNNNTRKSIRRGPGRAFFTSGIKGDIRQAHC
27 CARPNNNTRKSIRIXRGXXXXXXTXIXXXXXXXXXHC
26 CARPNNNTRKSIXRGPGRAFFTSGIKGDIRQAHC
24 CARPNNNTRKSIRXXXGXXXAXXTXIXXXXXXXXXHC
24 CARPNNNTRKSIRIXXGXXXXXXTSIXXXXXXXXXHC
22 CARPNNNTRKSIXIXXGXGXXXXTXIXXXXXXXXXHC
total prot count: 10289

"""