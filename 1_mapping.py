"""
0) Report raw counts for FASTQs
1) Map fastq reads to consensus B reference sequences
    --> Results stored in prelim.sam
2) Split the prelim.sam into region-specific SAMs
3) Make a sample/region specific consensus
4) Remap fastq reads against sample/region specific references
"""

import os
import sys
from miseqUtils import samBitFlag
import subprocess


## settings
REMAP_THRESHOLD = 0.95 # what fraction of the data we should be able to remap
MAX_REMAPS = 3


## arguments
refpath = sys.argv[1]
f = sys.argv[2] # should be an absolute path
qCutoff = sys.argv[3]
mode = sys.argv[4]
is_t_primer = (sys.argv[5]=='1')


# report counts to file
outfile = open(f.replace('.fastq', '.counts'), 'w')
contam_file = open(f.replace('.fastq', '.Tcontaminants.fastq'), 'w')


def remap (f1, f2, samfile, ref, qCutoff=30):
    """
    Generate a sample-specific consensus sequence from a samtools
    PILEUP file, and then remap everything to this consensus as a
    reference sequence.  Returns paths to the SAM output and
    consensus sequence file.
    [samfile] = SAM file containing last mapping
    [ref] = last reference sequence used to generate SAM file
    """
    
    bamfile = samfile.replace('.sam', '.bam')
    confile = bamfile+'.pileup.conseq'
    
    # will overwrite last *.remap.sam file
    remapped_sam = samfile if samfile.endswith('.remap.sam') else samfile.replace('.sam', '.remap.sam')
    
    if ref == refpath:
        # using the default reference file (first run)
        cmd = 'samtools view -bt %s.fasta.fai %s > %s 2>/dev/null' % (ref, samfile, bamfile)
        os.system(cmd) # isolating 'cmd' string for debugging
    else:
        # make new samtools index file
        cmd = 'samtools faidx %s' % ref
        os.system(cmd) # creates [ref].fai
        cmd = 'samtools view -bt %s.fai %s > %s 2>/dev/null' % (ref, samfile, bamfile)
        os.system(cmd)
    
    cmd = 'samtools sort %s %s.sorted' % (bamfile, bamfile)
    os.system(cmd)
    cmd = 'samtools mpileup -A %s.sorted.bam > %s.pileup 2>/dev/null' % (bamfile, bamfile)
    os.system(cmd)
    
    # create new consensus sequence
    cmd = 'python2.7 pileup2conseq_v2.py %s.pileup %s' % (bamfile, qCutoff)
    os.system(cmd)
    
    # convert into bowtie2 reference files (*.bt2)
    cmd = 'bowtie2-build -q -f %s %s' % (confile, confile)
    os.system(cmd)
    
    # map original data to new reference
    cmd = 'bowtie2 --quiet -p 1 --local -x %s -1 %s -2 %s -S %s --no-unal --met-file %s --un %s --un-conc %s' % (confile, 
            f1, f2, remapped_sam, remapped_sam.replace('.sam', '.bt2_metrics'),
            remapped_sam.replace('.sam', '.bt2_unpaired_noalign.fastq'), 
            remapped_sam.replace('.sam', '.bt2_paired_noalign.fastq'))
    os.system(cmd)
            
    return remapped_sam, confile



# Get region codes from the con B reference fasta headers
infile = open(refpath+'.fasta', 'rU')
refnames = []
for line in infile:
    if line.startswith('>'): refnames.append(line.strip('>\n'))
infile.close()



# Path to R1 FASTQ
root = '/'.join(f.split('/')[:-1])
filename = f.split('/')[-1]
prefix = filename.split('.')[0].split('_')[0]
f1 = f
f2 = f.replace('R1', 'R2')


# raw count from FASTQs
p = subprocess.Popen(['wc', '-l', f], stdout = subprocess.PIPE)
stdout, stderr = p.communicate()
count1 = int(stdout.split()[0])/4

p = subprocess.Popen(['wc', '-l', f2], stdout = subprocess.PIPE)
stdout, stderr = p.communicate()
count2 = int(stdout.split()[0])/4

outfile.write('FASTQ,%d,%d\n' % (count1, count2))


# Initial consensus B mapping
samfile = '{}/{}.prelim.sam'.format(root, prefix)
command = 'bowtie2 --quiet -p 1 --local -x %s -1 %s -2 %s -S %s' % (refpath, f1, f2, samfile)
os.system(command)


# Create region-specific SAMs
refsams = {}
for i, refname in enumerate(refnames):
    refsams.update({refname: {'handle': open('%s/%s.%s.sam' % (root, prefix, refname), 'w'),'count': [0,0]}})
refsams.update({'*': {'handle': open('%s/%s.unmapped.sam' % (root, prefix), 'w'),'count': [0,0]}})


# Then, subdivide the prelim SAM into these region-specific SAMs
infile = open(samfile, 'rU')
line_counts = [0, 0]
t_counts = [0, 0]

for line in infile:
    # Copy the SAM header into each file
    if line.startswith('@'):
        for refname in refnames: refsams[refname]['handle'].write(line)
        continue
    
    qname, flag, refname, pos, mapq, cigar, rnext, pnext, tlen, seq, qual = line.strip('\n').split('\t')[:11]
    bitinfo = samBitFlag(flag)
    
    if bitinfo['is_unmapped']:
        # a read that is not mapped is always is_reverse == False
        continue
    
    if mode == 'Amplicon':
        # filter T contaminants
        if is_t_primer:
            # Trim T(A) when known to be present due to primer
            if bitinfo['is_reverse'] and seq.endswith('A'):
                seq = seq[:-1]
            elif not bitinfo['is_reverse'] and seq.startswith('T'):
                seq = seq[1:]
            else:
                t_counts[bitinfo['is_reverse']] += 1 # missing A/T
                contam_file.write('@%s\n%s\n+\n%s\n' % (qname, seq, qual))
                continue
        else:
            # look for contaminating T primed reads in non-T run
            if bitinfo['is_reverse'] and seq.endswith('A') and cigar.endswith('1S'):
                t_counts[bitinfo['is_reverse']] += 1 # unexpected A
                contam_file.write('@%s\n%s\n+\n%s\n' % (qname, seq, qual))
                continue
            elif not bitinfo['is_reverse'] and seq.startswith('T') and cigar.startswith('1S'):
                t_counts[bitinfo['is_reverse']] += 1 # unexpected T
                contam_file.write('@%s\n%s\n+\n%s\n' % (qname, seq, qual))
                continue
            else:
                pass
        
    # write line out to respective region-specific file
    try:
        refsams[refname]['handle'].write(line)
        refsams[refname]['count'][bitinfo['is_reverse']] += 1
    except:
        print "{} appeared to be an incorrect refname".format(refname)
    
    line_counts[bitinfo['is_reverse']] += 1

infile.close()
contam_file.close()


outfile.write('preliminary map,%d,%d\n' % tuple(line_counts))
if is_t_primer:
    outfile.write('missing T,%d,%d\n' % tuple(t_counts))
else:
    outfile.write('unexpected T,%d,%d\n' % tuple(t_counts))


for refname in refnames:
    refsams[refname]['handle'].close()
    outfile.write('prelim %s,%d,%d\n' % (refname, refsams[refname]['count'][0],
        refsams[refname]['count'][1]))


# FIXME: Q CUTOFFS NOT WORKING RIGHT NOW
# Remap fastqs using sample/region specific consensus reference
total_remap = 0

for refname in refnames:
    if sum(refsams[refname]['count']) == 0 or refname == 'phiX174' or refname == '*':
        continue
    
    samfile = refsams[refname]['handle'].name
    samfile, confile = remap(f1, f2, samfile, refpath, qCutoff)
    
    # keep track of file paths
    refsams[refname].update({'samfile': samfile, 'confile': confile})
    
    # get number of reads in each remap
    p = subprocess.Popen(['wc', '-l', samfile], stdout = subprocess.PIPE)
    stdout, stderr = p.communicate()
    count = (int(stdout.split()[0]) - 3) / 2. # first three lines contain comments
    total_remap += count
    refsams[refname]['count'][0] = count # reuse dictionary
    
    outfile.write('remap %s,%d\n' % (refname, int(count)))


# continue to remap if we've failed to map enough reads
if total_remap / count1 < REMAP_THRESHOLD:
    break_out = False
    for iter in range(MAX_REMAPS):
        total_remap = 0
        for refname in refnames:
            if sum(refsams[refname]['count']) == 0 or refname == 'phiX174' or refname == '*':
                continue
            
            samfile = refsams[refname]['samfile']
            confile = refsams[refname]['confile']
            samfile, confile = remap(f1, f2, samfile, confile, qCutoff)
            refsams[refname]['samfile'] = samfile
            refsams[refname]['confile'] = confile
            
            p = subprocess.Popen(['wc', '-l', samfile], stdout = subprocess.PIPE)
            stdout, stderr = p.communicate()
            count = (int(stdout.split()[0]) - 3) / 2. # first three lines contain comments
            
            # if number of reads declines with mapping, let's quit
            if count < refsams[refname]['count'][0]:
                break_out = True
                break
            
            total_remap += count
            refsams[refname]['count'][0] = count
            outfile.write('remap %d %s,%d\n' % (iter, refname, int(count)))
        
        if break_out:
            break
        
        if total_remap / count1 >= REMAP_THRESHOLD:
            break


outfile.close()
