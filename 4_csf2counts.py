"""
Calculate nucleotide and amino acid counts from a FASTA or CSF file
"""
import csv,os,sys
from miseqUtils import *
from hyphyAlign import *

hyphy = HyPhy._THyPhy (os.getcwd(), 1)
change_settings(hyphy)

amino_alphabet = 'ACDEFGHIKLMNPQRSTVWY*'
cutoffs = [0.01, 0.02, 0.05, 0.1, 0.2, 0.25]        # What cutoffs are these???

# Load in the HXB2 amino reference sequences
f = open("HXB2_amino_sequences.csv", "rb")
input_file = csv.reader(f)
hxb2 = {}
for row in input_file:
    region, amino = row
    hxb2[region] = amino
f.close()

if len(sys.argv) < 3:
	print 'Usage: python 4_csf2counts.py /path/to/.csf|.fasta Nextera|Amplicon'
	sys.exit()

path = sys.argv[1]
filename = path.split('/')[-1]
sample, ref = filename.split('.')[:2]

mode = sys.argv[2]
assert mode in ['Nextera', 'Amplicon'], 'ERROR: Unrecognized mode (%s) in 4_csf2counts.py' % mode


# make the output stem by removing the extension of the filename
root = '/'.join(path.split('/')[:-1])
if root == '':
    # in case script is executed on file in cwd
    root = '.'

outpath = root + '/' + (filename.replace('.fasta', '') if mode == 'Amplicon' else filename.replace('.csf', ''))


# output to files and compute consensus
nucfile = open(outpath+'.nuc.csv', 'w')
nucfile.write('hxb2.nuc.pos,A,C,G,T\n')

aafile = open(outpath+'.amino.csv', 'w')
aafile.write('hxb2.aa.pos,%s\n' % ','.join(list(amino_alphabet)))

confile = open(outpath+'.conseq', 'w')
indelfile = open(outpath+'.indels.csv', 'w')


#assert hxb2.has_key(ref), 'Unknown reference sequence, expecting HXB2 gene'
if not hxb2.has_key(ref):
    sys.exit()

refseq = hxb2[ref]


infile = open(path, 'rU')
if mode == 'Nextera':
    fasta, lefts, rights = convert_csf(infile.readlines())
elif mode == 'Amplicon':
    fasta = convert_fasta(infile.readlines())
else:
    print 'This should never happen'
    sys.exit()

infile.close()


# use the first read to determine reading frame
max_score = 0
best_frame = 0
for frame in range(3):
    prefix = ('-'*lefts[fasta[0][0]] if mode == 'Nextera' else '')
    p = translate_nuc(prefix + fasta[0][1], frame)
    aquery, aref, ascore = pair_align(hyphy, refseq, p)
    if ascore > max_score:
        best_frame = frame # the reading frame of left = 0
        max_score = ascore

print "Best ORF = %d" % best_frame

# Iterate through reads and count WHAT?
nucs = {}
aminos = {}

# Cache protein sequences
pcache = []

# For each sequence
for i, (h, s) in enumerate(fasta):
    left = lefts[h] if mode == 'Nextera' else 0

    # FASTA headers contain read count in last entry
    count = 1 if mode == 'Nextera' else int(h.split('_')[-1])

    for j, nuc in enumerate(s):
        pos = left + j
        if not nucs.has_key(pos):
            nucs.update({pos: {}})
        if not nucs[pos].has_key(nuc):
            nucs[pos].update({nuc: 0})
        nucs[pos][nuc] += count
    
    p = translate_nuc('-'*left + s, best_frame)
    pcache.append(p)
    
    for pos, aa in enumerate(p):
        if aa == '-':
            continue
        #pos = left/3 + j
        if not aminos.has_key(pos):
            aminos.update({pos: {}})
        if not aminos[pos].has_key(aa):
            aminos[pos].update({aa: 0})
        aminos[pos][aa] += count


#print "Finished translating + determining nucleotide/amino counts"

# Generate AA plurality (max) consensus
keys = aminos.keys()
keys.sort()
aa_max = ''

for pos in keys:
    intermed = [(v, k) for k, v in aminos[pos].iteritems()]
    intermed.sort(reverse=True)
    aa_max += intermed[0][1]



#print "Generated protein plurality for {}: {}".format(sample,aa_max)

# Align consensus against HXB2
if ref == 'HIV1B-pol':
    # use PR-RT as reference instead of full length pol
    refseq = hxb2['HIV1B-prrt']

aquery, aref, ascore = pair_align(hyphy, refseq, aa_max)


#print "HXB2 aligned the amino plurality: {}".format(aquery)

# Ignore parts of query that fall outside our reference
# this will be important for pol since we are using shorter PR-RT as ref
left, right = get_boundaries(aref)
qindex_to_hxb2 = {} # Map from query to HXB2 coordinates
inserts = []        # Leep track of which aa positions are insertions
qindex = 0
rindex = 0
for i in range(len(aref)):
    # ignore parts of query that do not overlap reference
    if i < left:
        qindex += 1
        continue
    if i >= right:
        break
    
    if aref[i] == '-':
        # insertion in query
        inserts.append(qindex)
        qindex += 1
    elif aquery[i] == '-':
        # deletion in query
        # do not increment qindex
        rindex += 1
        continue
    else:
        qindex_to_hxb2.update({qindex: rindex})
        qindex += 1
        rindex += 1

#print "Determined query-to-hxb2 index: {}".format(qindex_to_hxb2)

# then reiterate through sequences to capture indels
if len(inserts) > 0:
    indelfile.write('insertion,count\n')
    indel_counts = {}
    
    for p in pcache:
        ins_str = str(inserts[0])
        last_i = -1
        for i in inserts:
            if last_i > -1 and i - last_i > 1:
                # end of a contiguous indel
                ins_str += ',%d' % i
            try:
                ins_str += p[i]
            except IndexError:
                break
            last_i = i
        
        if not indel_counts.has_key(ins_str):
            indel_counts.update({ins_str: 0})
        
        indel_counts[ins_str] += 1
    
    for ins_str, count in indel_counts.iteritems():
        indelfile.write('%s,%d\n' % (ins_str, count))

indelfile.close()


# output nucleotide and amino acid counts in HXB2 coordinates
# also output consensus sequences at varying thresholds
keys = nucs.keys()
keys.sort()
codon_pos = 0
maxcon = ''
conseqs = ['' for cut in cutoffs]

for pos in keys:
    aapos = pos/3
    """
    if aapos in inserts:
        continue
    """
    try:
        hxb2_pos = qindex_to_hxb2[aapos] + 1
    except KeyError:
        continue
    
    codon_pos += 1
    if codon_pos > 2:
        codon_pos = 0
    
    nucfile.write('%d,%s\n' % (3*hxb2_pos + codon_pos, 
        ','.join(map(str, [nucs[pos].get(nuc, 0) for nuc in 'ACGT']))))
    
    # plurality consensus
    intermed = [(count, nuc) for nuc, count in nucs[pos].iteritems()]
    intermed.sort(reverse=True)
    maxcon += intermed[0][1]
    
    # consensuses with mixtures determined by frequency cutoff
    total_count = sum([count for count, nuc in intermed])
    
    for ci, cutoff in enumerate(cutoffs):
        mixture = []
        for count, nuc in intermed:
            if float(count) / total_count > cutoff:
                mixture.append(nuc)
        
        if 'N' in mixture:
            if len(mixture) > 1:
                mixture.remove('N')
            else:
                # completely ambiguous
                conseqs[ci] += 'N'
                continue
        
        if '-' in mixture:
            if len(mixture) > 1:
                mixture.remove('-')
            else:
                conseqs[ci] += '-'
                continue
        
        # encode nucleotide mixture
        if len(mixture) > 1:
            mixture.sort()
            conseqs[ci] += ambig_dict[''.join(mixture)]
        elif len(mixture) == 1:
            conseqs[ci] += mixture[0]
        else:
            # mixture of length zero, no bases exceed cutoff
            conseqs[ci] += 'N'

# output consensus sequences
confile.write('>%s_MAX\n%s\n' % (sample, maxcon))

for ci, cutoff in enumerate(cutoffs):
    confile.write('>%s_%1.3f\n%s\n' % (sample, cutoff, conseqs[ci]))

nucfile.close()
confile.close()


#print "Wrote consensus"

# output amino acid counts
keys = aminos.keys()
keys.sort()
for aapos in keys:
    if aapos in inserts:
        continue
    try:
        hxb2_pos = qindex_to_hxb2[aapos] + 1
    except KeyError:
        continue
        
    aafile.write('%d,%s\n' % (hxb2_pos, 
        ','.join(map(str, [aminos[aapos].get(aa, 0) for aa in amino_alphabet]))))
    

aafile.close()

#print "Wrote amino counts"

