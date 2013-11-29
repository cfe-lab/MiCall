"""
1) Report raw fastq read counts
2) Map fastq reads to consensus B reference sequences (Store in prelim.sam)
3) Split prelim.sam into region-specific SAMs
4) Make a sample/region specific consensus
5) Remap fastq reads against sample/region specific references
"""

import os,subprocess,sys
from miseqUtils import samBitFlag

## Arguments
conbrefpath = sys.argv[1]             # Consensus B reference path
R1_fastq = sys.argv[2]                # Absolute path to R1 fastq file
qCutoff = sys.argv[3]                 # Q cutoff for pileup2conseq
mode = sys.argv[4]
is_t_primer = (sys.argv[5]=='1')      # Used for contamination detection
REMAP_THRESHOLD = float(sys.argv[6])  # Minimum tolerated fraction of mapped fastq reads
MAX_REMAPS = int(sys.argv[7])         # Number of attempts at remapping before giving up

def remap (R1_fastq, R2_fastq, samfile, ref, qCutoff=30):
    """
    1) Generate sample-specific consensus from a samtools pileup.
    2) Remap everything to this consensus as a ref seq
    3) Returns paths to the SAM output and consensus sequence file
    """
    
    bamfile = samfile.replace('.sam', '.bam')
    confile = bamfile+'.pileup.conseq'

    # Overwrite previous iterations of the remapped sam (Only keep the original specific prelim sam)
    remapped_sam = samfile if samfile.endswith('.remap.sam') else samfile.replace('.sam', '.remap.sam')

    # If this is the first run, use the default con B reference
    if ref == conbrefpath:
        cmd = 'samtools view -bt %s.fasta.fai %s > %s 2>/dev/null' % (ref, samfile, bamfile)
        os.system(cmd)
    else:
        # make new samtools index file
        cmd = 'samtools faidx %s' % ref
        os.system(cmd) # creates [ref].fai
        cmd = 'samtools view -bt %s.fai %s > %s 2>/dev/null' % (ref, samfile, bamfile)
        os.system(cmd)

    # Sort the bam file by leftmost position on the reference assembly
    cmd = 'samtools sort %s %s.sorted' % (bamfile, bamfile)
    os.system(cmd)

    # Make a pileup from the sorted bam
    cmd = 'samtools mpileup -A %s.sorted.bam > %s.pileup 2>/dev/null' % (bamfile, bamfile)
    os.system(cmd)
    
    # Create new consensus sequence from the pileup
    cmd = 'python2.7 pileup2conseq_v2.py %s.pileup %s' % (bamfile, qCutoff)
    os.system(cmd)
    
    # Convert consensus into bowtie2 reference files (*.bt2)
    cmd = 'bowtie2-build -q -f %s %s' % (confile, confile)
    os.system(cmd)
    
    # Map original fastq reads to new reference
    cmd = 'bowtie2 --quiet -p 1 --local -x %s -1 %s -2 %s -S %s --no-unal --met-file %s --un %s --un-conc %s' % (confile, 
            R1_fastq, R2_fastq, remapped_sam, remapped_sam.replace('.sam', '.bt2_metrics'),
            remapped_sam.replace('.sam', '.bt2_unpaired_noalign.fastq'), 
            remapped_sam.replace('.sam', '.bt2_paired_noalign.fastq'))
    os.system(cmd)
            
    return remapped_sam, confile



# Store region codes in con B reference fasta headers in refnames
with open(conbrefpath+'.fasta', 'rU') as conb_ref_fasta:
    refnames = []
    for line in conb_ref_fasta:
        if line.startswith('>'): refnames.append(line.strip('>\n'))

# Deduce R1/R2 file pairing
root = '/'.join(R1_fastq.split('/')[:-1])
filename = os.path.basename(R1_fastq)			# Filename of R1 fastq...
file_prefix = filename.split('.')[0]			# Has a prefix containing...
sample_name, sample_well = file_prefix.split('_')[:2]	# Sample name and well
prefix = '_'.join([sample_name, sample_well])		# Change this prefix to be joined by _ (??)
count_file = open(R1_fastq.replace('.fastq', '.counts'), 'w')
R2_fastq = R1_fastq.replace('R1', 'R2')

# Determine the number of reads in both (R1 + R2) fastq files, store in the .count file
p = subprocess.Popen(['wc', '-l', R1_fastq], stdout = subprocess.PIPE)
stdout, stderr = p.communicate()
total_reads_R1 = int(stdout.split()[0])/4
p = subprocess.Popen(['wc', '-l', R2_fastq], stdout = subprocess.PIPE)
stdout, stderr = p.communicate()
total_reads_R2 = int(stdout.split()[0])/4
count_file.write('Raw FASTQ,R1,%d,R2,%d\n' % (total_reads_R1, total_reads_R2))

# Initial consensus B mapping
prelim_samfile = '{}/{}.prelim.sam'.format(root, prefix)
command = 'bowtie2 --quiet -p 1 --local -x %s -1 %s -2 %s -S %s' % (conbrefpath, R1_fastq, R2_fastq, prelim_samfile)
os.system(command)

# Define region-specific SAMs: refsams[refname] points to a nested dict which includes a file handle to each specific SAM
refsams = {}
for i, refname in enumerate(refnames):
    region_specific_sam = '%s/%s.%s.sam' % (root, prefix, refname)
    refsams.update({refname: {'sam_file_handle': open(region_specific_sam, 'w'),'count': [0,0]}})
refsams.update({'*': {'sam_file_handle': open('%s/%s.unmapped.sam' % (root, prefix), 'w'),'count': [0,0]}})

# Subdivide prelim SAMs into region-specific SAMs
prelim_sam_infile = open(prelim_samfile, 'rU')
line_counts = [0, 0]
t_counts = [0, 0]
contam_file = open(R1_fastq.replace('.fastq', '.Tcontaminants.fastq'), 'w')

for line in prelim_sam_infile:

    # Copy the original SAM header into each region-specific SAM
    if line.startswith('@'):
        for refname in refnames: refsams[refname]['sam_file_handle'].write(line)
        continue

    # SAM documentation explains what these fields mean - http://samtools.sourceforge.net/SAMv1.pdf
    qname, flag, refname, pos, mapq, cigar, rnext, pnext, tlen, seq, qual = line.strip('\n').split('\t')[:11]

    # SAM bitwise flag variable specifies whether read is paired, successfully mapped, etc
    bitinfo = samBitFlag(flag)

    # Ignore unmapped reads
    if bitinfo['is_unmapped']:
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

                # Dump in fastq format to contaminant file
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
        
    # Write each line into respective region-specific SAMs
    try:
        refsams[refname]['sam_file_handle'].write(line)
        refsams[refname]['count'][bitinfo['is_reverse']] += 1
    except:
        print "{} appeared to be an incorrect refname".format(refname)

    # I think is_reverse identifies R1 reads from R2 in the prelim SAM...? 
    line_counts[bitinfo['is_reverse']] += 1

prelim_sam_infile.close()
contam_file.close()

# FIXME: CONFIRM WITH ART IF THIS IS TRUE...?
# Show the number of reads that made it to preliminary mapping in the count file
count_file.write('Preliminary map,{},{}\n'.format(line_counts[0], line_counts[1]))
if is_t_primer: error_type = "Missing T"
else: error_type = "Unexpected T"
count_file.write('{} (Detected at prelim map),R1,{},R2,{}\n'.format(error_type, t_counts[0], t_counts[1]))

# Show the reads mapped at preliminary mapping with respect to region
for refname in refnames:
    refsams[refname]['sam_file_handle'].close()
    count_file.write('prelim %s,R1,%d,R2,%d\n' % (refname, refsams[refname]['count'][0],refsams[refname]['count'][1]))

# Track the number of reads that have been successfully remapped in total
total_remap = 0

# Remap the fastqs using sample/region specific conseqs
for refname in refnames:

    # Completely ignore phiX, unmapped reads, and regions which had no mapping at the prelim mapping stage
    if sum(refsams[refname]['count']) == 0 or refname == 'phiX174' or refname == '*':
        continue

    # Run remap on the region-specific sam, and get the remapped sam and consensus pileup used to generate it
    samfile = refsams[refname]['sam_file_handle'].name
    samfile, confile = remap(R1_fastq, R2_fastq, samfile, conbrefpath, qCutoff)
    
    # Track file paths
    refsams[refname].update({'samfile': samfile, 'confile': confile})
    
    # Track the number of mapped reads in each region
    p = subprocess.Popen(['wc', '-l', samfile], stdout = subprocess.PIPE)
    stdout, stderr = p.communicate()
    region_specific_count = (int(stdout.split()[0]) - 3) / 2. # First 3 lines contain comments
    total_remap += region_specific_count
    refsams[refname]['count'][0] = region_specific_count
    count_file.write('remap %s,%d\n' % (refname, int(region_specific_count)))

# Continue to remap if we've failed to map enough total reads
if total_remap / total_reads_R1 < REMAP_THRESHOLD:
    break_out = False

    # Repeat remapping up to the number of MAX_REMAPS permitted
    for iter in range(MAX_REMAPS):
        total_remap = 0
        for refname in refnames:
            if refsams[refname]['count'][0] == 0 or refname == 'phiX174' or refname == '*':
                continue
            
            samfile = refsams[refname]['samfile']
            confile = refsams[refname]['confile']
            samfile, confile = remap(R1_fastq, R2_fastq, samfile, confile, qCutoff)
            refsams[refname]['samfile'] = samfile
            refsams[refname]['confile'] = confile
            
            # Continue to determine the number of reads mapped in this region
            p = subprocess.Popen(['wc', '-l', samfile], stdout = subprocess.PIPE)
            stdout, stderr = p.communicate()
            region_specific_count = (int(stdout.split()[0]) - 3) / 2.
            
            # If number of reads is declining, that's really bad
            if region_specific_count < refsams[refname]['count'][0]:
                break_out = True
                break
            
            total_remap += region_specific_count
            refsams[refname]['count'][0] = region_specific_count

            # Write each remap iteration result to the file
            count_file.write('remap %d %s,%d\n' % (iter, refname, int(region_specific_count)))

        if break_out or total_remap / total_reads_R1 >= REMAP_THRESHOLD:
            break

count_file.close()
