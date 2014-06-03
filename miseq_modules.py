def collate_counts(run_path,output_path):
    """
    Collate .counts files into a single CSV summary file.
    """
    import glob,os

    with open(output_path, "w") as collated_file:
        for count_file in [f for f in glob.glob("{}/*.counts".format(run_path))]:
            prefix = (os.path.basename(count_file)).rstrip(".counts")
            with open(count_file,"r") as f_in:
                for i, line in enumerate(f_in.readlines()):
                    collated_file.write("{},{}".format(prefix,line))

def collate_conseqs(run_path,output_path):
    """
    Collate .conseq files into a single CSV summary file.
    """
    import glob,os
    from miseqUtils import convert_fasta

    files = [f for f in glob.glob("{}/*.conseq".format(run_path)) if 'pileup' not in f]

    with open(output_path,"w") as f_out:
        f_out.write("sample,region,q-cutoff,s-number,consensus-percent-cutoff,sequence\n")

        for path in files:
            prefix = (os.path.basename(path)).rstrip(".conseq")
            sample, region, q = prefix.split(".")[:3]

            with open(path,"r") as f:
                fasta = convert_fasta(f.readlines())

            for header, sequence in fasta:
                fasta_sample, s_number, consensus_percentage = header.split("_")
                f_out.write("{},{},{},{},{},{}\n".format(fasta_sample,region,q,s_number,consensus_percentage,sequence))

def collate_frequencies (run_path, output_path, type):
    """
    Collate amino/nuc .freq files into a single summary file.
    """
    import glob,os

    if type == "amino":
        file_extension = "amino.freqs"
        header = "sample,region,q-cutoff,query.aa.pos,refseq.aa.pos,A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,*"
    elif type == "nuc":
        file_extension = "nuc.freqs"
        header = "sample,region,q-cutoff,query.nuc.pos,refseq.nuc.pos,A,C,G,T"
    else:
        raise Exception("Incorrect type parameter")

    # Summarize uncleaned freqs
    with open(output_path, "w") as f_out:
        f_out.write("{}\n".format(header))
        for path in [x for x in glob.glob("{}/*.{}".format(run_path,file_extension)) if "clean.{}".format(file_extension) not in x]:
            prefix = (os.path.basename(path)).rstrip(".{}".format(file_extension))
            sample, region, q = prefix.split(".")[:3]
            with open(path,"r") as f:
                lines = f.readlines()

            for j, line in enumerate(lines):
                if j == 0:
                    continue
                f_out.write("{},{},{},{}\n".format(sample, region, q, line.rstrip("\n")))

    # Summarize clean freqs
    file_extension = "clean.{}".format(file_extension)

    output_path = output_path.replace("_frequencies.csv", "_cleaned_frequencies.csv")
    with open(output_path, "w") as f_out:
        f_out.write("{}\n".format(header))
        for path in glob.glob("{}/*.{}".format(run_path,file_extension)):
            prefix = (os.path.basename(path)).rstrip(".{}".format(file_extension))
            sample, region, q = prefix.split(".")[:3]
            with open(path,"r") as f:
                lines = f.readlines()

            for j, line in enumerate(lines):
                if j == 0:
                    continue
                f_out.write("{},{},{},{}\n".format(sample, region, q, line.rstrip("\n")))


def csf2nuc(path, mode, nuc_reference_file, min_avg_score=2.):
    """
    Generate a FASTA where reads have been aligned and trimmed to some
    nucleotide reference sequence.  FASTA headers contain index and count
    of reads.
    """
    import os, csv, logging, HyPhy, hyphyAlign
    from miseqUtils import parse_csf
    logger = logging.getLogger()

    hyphy = HyPhy._THyPhy(os.getcwd(), 1)
    hyphyAlign.change_settings(hyphy, alphabet=hyphyAlign.nucAlphabet,
                               scoreMatrix=hyphyAlign.nucScoreMatrix,
                               gapOpen=20, gapOpen2=20,
                               gapExtend=10, gapExtend2=10,
                               noTerminalPenalty=1)

    # set up file paths
    root = os.path.dirname(path)
    if root == '':
        root = '.'

    filename = os.path.basename(path)
    if not filename.endswith('.csf'):
        return logger.error('{} is not a CSF file - exiting csf_to_fasta'.format(path))
    file_prefix = filename.replace('.csf', '')
    outpath = root+'/'+file_prefix

    # 'ref' corresponds to one of the target/region-specific references, e.g., HIV1-nef
    sample, ref = filename.split('.')[:2]

    # parse nucleotide reference sequences
    refseqs = {}
    with open(nuc_reference_file, 'rb') as f:
        rows = csv.reader(f)
        for region, variant, subregion, sequence in rows:
            if region not in refseqs:
                refseqs.update({region: {}})
            # region should correspond to CfE reference notation
            refseqs[region][subregion] = sequence

    if ref not in refseqs:
        # no nuc ref specified for this region
        return

    # prepare output files
    outfiles = {}
    for subregion in refseqs[ref].iterkeys():
        handle = open(outpath+'.'+subregion+'.nuc', 'w')
        handle.write('>Reference_%s_%s\n%s\n' % (ref, subregion, refseqs[region][subregion]))
        outfiles.update({subregion: {'handle': handle, 'count': 0}})

    # do quick and dirty alignment of each read
    with open(path, 'rU') as f:
        csf = parse_csf(f, mode)
        for items in csf:
            header, seq = items[:2]
            for subregion, refseq in refseqs[ref].iteritems():
                aquery, aref, ascore = hyphyAlign.pair_align(hyphy, refseq, seq)
                if float(ascore)/len(aref) < min_avg_score:
                    continue

                # output portion of read that overlaps this subregion
                left, right = hyphyAlign.get_boundaries(aref)
                outfiles[subregion]['handle'].write('>%s\n%s\n' % (header, aquery[left:right]))
                outfiles[subregion]['count'] += 1

    for val in outfiles.itervalues():
        val['handle'].close()
        if val['count'] == 0:
            # empty file except for reference sequence, delete
            os.remove(val['handle'].name)



def csf2counts(path, mode, mixture_cutoffs, amino_reference_sequence):
    """
    Calculate HXB2-aligned nucleotide and amino acid counts from a CSF.
    Assumes that there is a consistent reading frame over the entire alignment.
    """

    import csv, logging, HyPhy, os
    from hyphyAlign import change_settings, get_boundaries, pair_align
    from miseqUtils import ambig_dict, parse_csf, convert_fasta, mixture_dict, translate_nuc

    logger = logging.getLogger()
    hyphy = HyPhy._THyPhy (os.getcwd(), 1)
    change_settings(hyphy) # default gap open penalty 40(20), extension penalty 10(5) - we may need to change these

    amino_alphabet = 'ACDEFGHIKLMNPQRSTVWY*'

    if mode not in ['Amplicon', 'Nextera']:
        return logger.error("{} is an unsupported mode - halting csf2counts".format(mode))

    # set up file paths
    filename = os.path.basename(path)
    root = os.path.dirname(path) if os.path.dirname(path) != '' else '.'
    file_prefix = filename.replace('.csf', '')
    outpath = "{}/{}".format(root, file_prefix)

    # CSF contains sample + region in filename (Ex: F00844_S68.HIV1B-pol.0.csf)
    sample, ref = filename.split('.')[:2]

    # Amino reference sequences in refseqs is used to coordinate normalize our samples
    with open(amino_reference_sequence, "rb") as f:
        input_file = csv.reader(f)
        refseqs = {}
        for row in input_file:
            region, amino = row
            refseqs[region] = amino

    # If we have no reference sequence, we can't align the input sequences
    if ref not in refseqs:
        logger.error("No reference for {} - halting csf2counts".format(ref))
        return

    refseq = refseqs[ref]

    # load region- and sample-specific reference
    try:
        with open(root+'/'+sample+'.'+ref+'.bam.pileup.conseq', 'rU') as handle:
            pileup_conseq = convert_fasta(handle.readlines())[0][1].upper()
    except:
        logger.error('No sample/region-specific consensus sequence for %s, %s' % (sample, ref))
        return

    # use reference to determine reading frame
    max_score = -999
    best_frame = 0
    for frame in range(3):
        p = translate_nuc(pileup_conseq, frame)
        aquery, aref, ascore = pair_align(hyphy, refseq, p)
        if ascore > max_score:
            best_frame = frame
            max_score = ascore

    logging.debug('Best ORF = %d' % best_frame)

    nuc_counts = {} # Base counts by self-consensus coordinate
    aa_counts = {}	# Amino counts by self-consensus coordinate
    pcache = []		# Cache protein sequences


    # CSF reads aligned against self-consensus: offset is with respect to self
    # Load CSF (CSF header, offset, sequence) into fasta data structure
    min_left = 1e6
    with open(path, 'rU') as f:
        csf = parse_csf(f, mode)
        for header, seq, count, left, right in csf:
            if left < min_left:
                min_left = left

            # Determine nuc counts with respect to self-consensus coordinates
            for j, nuc in enumerate(seq):
                pos = left + j
                if pos not in nuc_counts:
                    nuc_counts.update({pos: {}})
                if nuc not in nuc_counts[pos]:
                    nuc_counts[pos].update({nuc: 0})
                nuc_counts[pos][nuc] += count

            # Determine amino counts with respect to self-consensus coordinates
            p = translate_nuc('-'*left + seq, best_frame)
            pcache.append(p)

            # boundaries of read in amino acid space
            aa_left = (left + best_frame) / 3
            aa_right = (right + left + best_frame) / 3

            for pos, aa in enumerate(p):
                # Do not store gap information
                #if aa == '-':
                #	continue
                if pos < aa_left or pos >= aa_right:
                    continue
                if pos not in aa_counts:
                    aa_counts.update({pos: {}})
                if aa not in aa_counts[pos]:
                    aa_counts[pos].update({aa: 0})
                aa_counts[pos][aa] += count

    if len(nuc_counts) == 0:
        # skip empty file
        logger.error('{} is an empty file'.format(filename))
        return

    # Generate amino plurality consensus for query to reference coordinate mapping
    aa_coords = aa_counts.keys()
    aa_coords.sort()
    aa_max = ''
    for pos in range(min(aa_coords), max(aa_coords)+1):
        if pos in aa_coords:
            intermed = [(aa_count, amino) for amino, aa_count in aa_counts[pos].iteritems()]
            intermed.sort(reverse=True)
            aa_max += intermed[0][1]
        else:
            aa_max += '?'  # no coverage but not a gap

    logger.debug("Amino plurality consensus = {}".format(aa_max))

    aquery, aref, ascore = pair_align(hyphy, refseq, aa_max)
    left, right = get_boundaries(aref)	# Coords of first/last non-gap character

    logger.debug("Aligned amino plurality conseq = {}".format(aquery))
    logger.debug("Aligned reference sequence = {}".format(aref))

    ## Generate coordinate map to HXB2 reference
    qindex_to_refcoord = {}			# Query <-> reference coordinate mapping
    inserts = []					# Keep track of which aa positions are insertions
    qindex = 0						# Where we are in the query?
    rindex = 0						# Where we are in the reference?
    ref_coords = range(len(aref))

    logger.debug("Making qindex_to_refcoord")

    # For each coordinate on the reference, create a mapping to the query
    for i in ref_coords:
        # Do not consider parts of the query outside of the reference
        if i < left:
            qindex += 1
        elif i >= right:
            break
        # A gap in the reference is an insertion in the query which we want to skip in the mapping
        elif aref[i] == '-':
            inserts.append(qindex)  # Store insert location in query coordinate space
            qindex += 1				# Track along the query
        # If theres a gap in the query we are only effectively tracking along the pre-alignment reference
        elif aquery[i] == '-':
            rindex += 1
        # Normal case: tracking forward on both sequences
        else:
            qindex_to_refcoord[qindex] = rindex
            qindex += 1
            rindex += 1

    logger.debug('qindex_to_refcoord: ' + str(qindex_to_refcoord))#"qindex_to_refcoord {}".format(qindex_to_refcoord))

    # Write inserts to an indels.csv file
    if len(inserts) > 0:
        with open(outpath+".indels.csv", 'w') as indelfile:
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

    # Initialize initial (blank) consensus sequence for each mixture rule
    maxcon = ''
    conseqs = ['' for cut in mixture_cutoffs]
    query_codon_pos = 0
    nuc_coords = nuc_counts.keys()		  # nucs[self-coord][nuc] = count
    nuc_coords.sort()

    # account for assembly offset due to extra bases in sample-specific consensus
    nuc_assembly_offset = min_left

    # Output nucleotide counts in reference coordinate space to nuc.csv files
    logging.debug("Making {}".format(outpath+'.nuc.freqs'))
    nucfile = open(outpath+'.nuc.freqs', 'w')#open("{}.nuc.csv".format(outpath), 'w')
    nucfile.write("query.nuc.pos,refSeq.nuc.pos,A,C,G,T\n")

    for query_nuc_pos in nuc_coords:
        nucleotide_counts = [nuc_counts[query_nuc_pos].get(nuc, 0) for nuc in 'ACGT']
        nucleotide_counts_string = ','.join(map(str, nucleotide_counts))

        # Convert nucleotide query index into reference index
        try:
            # best frame is adjusted by shift from query to assembly coordinates
            adjustment = best_frame - (3 - nuc_assembly_offset%3)%3
            query_aa_pos = (query_nuc_pos - nuc_assembly_offset + adjustment) / 3
            query_codon_pos = (query_nuc_pos - nuc_assembly_offset + adjustment) % 3
            ref_aa_pos = qindex_to_refcoord[query_aa_pos]
            ref_nuc_pos = 3*ref_aa_pos + query_codon_pos

            nucfile.write(','.join(map(str, [query_nuc_pos+1, ref_nuc_pos+1, nucleotide_counts_string])))#"{},{},{}\n".format(query_nuc_pos, ref_nuc_pos, nucleotide_counts_string))
            nucfile.write('\n')

        except KeyError:
            logger.debug('No coordinate mapping for query nuc %d / amino %d (%s)' % (query_nuc_pos, query_aa_pos, filename))
            continue


        # Store self-aligned nucleotide plurality conseqs
        intermed = [(count, nuc) for nuc, count in nuc_counts[query_nuc_pos].iteritems()]
        intermed.sort(reverse=True)
        maxcon += intermed[0][1]

        # Determine the number of bases in total at this query position
        total_count = sum([count for count, nuc in intermed])

        for ci, mixture_cutoff in enumerate(mixture_cutoffs):
            mixture = []

            # If a base is greater than the proportion cutoff, the base contributes
            for count, nuc in intermed:
                if float(count) / total_count > mixture_cutoff:
                    mixture.append(nuc)

            # If an N exists with other bases, those bases take precedence
            if 'N' in mixture:
                if len(mixture) > 1:
                    mixture.remove('N')
                else:
                    conseqs[ci] += 'N'
                    #logger.debug("N was the majority base at position {} - {} (mixture_cutoff = {})".format(query_nuc_pos, filename, mixture_cutoff))
                    logger.debug("N was the majority base at position %d - %s (mixture_cutoff = %f)" % (query_nuc_pos, filename, mixture_cutoff))
                    continue

            # If there is a gap, but also bases, those bases take precedence
            if '-' in mixture:
                if len(mixture) > 1:
                    mixture.remove('-')
                else:
                    conseqs[ci] += '-'
                    continue

            # Attach mixture (If one exists) to the conseq with appropriate mixture cutoff rule
            if len(mixture) > 1:
                mixture.sort()
                conseqs[ci] += ambig_dict[''.join(mixture)]
            elif len(mixture) == 1:
                conseqs[ci] += mixture[0]
            else:
                # Mixture of length zero, no bases exceed cutoff
                conseqs[ci] += 'N'
    nucfile.close()

    # Store self-aligned plurality amino sequences in .conseq files
    #with open("{}.conseq".format(outpath), 'w') as confile:
    with open(outpath+'.conseq', 'w') as confile:
        confile.write('>%s_MAX\n%s\n' % (sample, maxcon))
        for ci, cutoff in enumerate(mixture_cutoffs):
            confile.write('>%s_%1.3f\n%s\n' % (sample, cutoff, conseqs[ci]))


    # Write amino acid counts in reference coordinate space in amino.csv files
    #with open("{}.amino.csv".format(outpath), 'w') as aafile:
    with open(outpath+".amino.freqs", 'w') as aafile:
        aafile.write("query.aa.pos,refseq.aa.pos,%s\n" % (','.join(list(amino_alphabet))))

        for qindex, ref_aa_pos in qindex_to_refcoord.iteritems():
            # adjust for assembly offset
            aa_pos = qindex + min(aa_coords)

            # Ignore query inserts
            if aa_pos in inserts:
                logger.debug("%d is an insert - ignoring" % (aa_pos))
                continue

            try:
                #ref_aa_pos = qindex_to_refcoord[aa_pos] + 1	 # FIXME: DO WE NEED TO ADD 1?
                aa_counts_string = ','.join(map(str, [aa_counts[aa_pos].get(aa, 0) for aa in amino_alphabet]))
                # note that we are subtracting the minimum aa_counts key
                aafile.write('%d,%d,%s\n' % (aa_pos, ref_aa_pos+1, aa_counts_string))

            except KeyError:
                logger.debug("No query-ref mapping available for aapos=%d (%s)" % (aa_pos, filename))

def system_call(command):
    import logging, subprocess
    logger = logging.getLogger()
    logger.info(command)
    subprocess.call(command, shell=True)

def remap (R1_fastq, R2_fastq, samfile, ref, original_reference, conseq_qCutoff=30, num_threads=1):
    """
    1) Generate sample-specific consensus from a samtools pileup.
    2) Remap everything to this consensus as a ref seq
    3) Returns paths to the SAM output and consensus sequence file
    """
    import logging
    from miseq_modules import system_call
    from miseqUtils import pileup_to_conseq

    logger = logging.getLogger()
    bamfile = samfile.replace('.sam', '.bam')
    #confile = "{}.pileup.conseq".format(bamfile)

    # Overwrite previous iterations of the remapped sam (Only keep the original specific prelim sam)
    remapped_sam = samfile if samfile.endswith('.remap.sam') else samfile.replace('.sam', '.remap.sam')

    # If this is the first run, use the static reference
    if ref == original_reference:
        system_call('samtools view -bt {}.fasta.fai {} > {} 2>/dev/null'.format(ref, samfile, bamfile))
    else:
        # Make new samtools index file (Creates [ref].fai)
        system_call('samtools faidx {}'.format(ref))
        system_call('samtools view -bt {}.fai {} > {} 2>/dev/null'.format(ref, samfile, bamfile))

    # Sort the bam file by leftmost position on the reference assembly
    system_call('samtools sort {} {}.sorted'.format(bamfile, bamfile))

    # Make a pileup from the sorted bam
    system_call('samtools mpileup -A {}.sorted.bam > {}.pileup 2>/dev/null'.format(bamfile, bamfile))

    # Create new consensus sequence from the pileup
    #system_call('python2.7 pileup2conseq_v2.py {}.pileup {}'.format(bamfile, conseq_qCutoff))
    pileup_file = bamfile + '.pileup'
    logging.info("pileup_to_conseq({},{})".format(pileup_file, conseq_qCutoff))
    confile = pileup_to_conseq(pileup_file, conseq_qCutoff)

    # Convert consensus into bowtie2 reference files (Creates 6 files of *.bt2)
    system_call('bowtie2-build -f -q {} {}'.format(confile, confile))

    # Map original fastq reads to new reference
    cmd = 'bowtie2 --quiet -p {} --local -x {} -1 {} -2 {} -S {} --no-unal'.format(num_threads,
            confile, R1_fastq, R2_fastq, remapped_sam)
    system_call(cmd)

    return remapped_sam, confile


def mapping(refpath, R1_fastq, conseq_qCutoff, mode, is_t_primer, REMAP_THRESHOLD, MAX_REMAPS, num_threads=1, filter_t_contaminants=False):
    """
    refpath		Absolute path to static reference sequences used during original mapping
    R1_fastq	Absolute path to R1 fastq file
    conseq_qCutoff	Q-cutoff for determining if a base contributes to the consensus in pileup2conseq
    is_t_primer Used to detect contaminants from previous runs (FIXME: REMOVE?)
    REMAP_THRESHOLD Efficiency of read mapping at which we stop bothering to perform additional remapping
    MAX_REMAPS	Number of extra attempts at remapping before giving up
    """
    import logging, os, subprocess, sys
    from miseqUtils import samBitFlag, parse_fasta
    from miseq_modules import system_call
    import re

    soft_clip = re.compile('^[0-9]+S')  # CIGAR string starts with soft clip

    logger = logging.getLogger()
    original_reference = refpath		# Path to the original reference sequence

    # Store region codes of static reference fasta in refnames (ConB for HIV, H77 for HCV)
    refseq_lengths = {}
    refnames = []
    handle = open(refpath+'.fasta', 'rU')
    ref_fasta = parse_fasta(handle)
    for h, s in ref_fasta:
        refnames.append(h)
        refseq_lengths.update({h: len(s)})


    # Deduce R1/R2 file pairing
    root = '/'.join(R1_fastq.split('/')[:-1])
    filename = os.path.basename(R1_fastq)			# Filename of R1 fastq...
    file_prefix = filename.split('.')[0]			# Has a prefix containing...
    sample_name, sample_well = file_prefix.split('_')[:2]	# Sample name and well
    prefix = '_'.join([sample_name, sample_well])		# Change this prefix to be joined by _ (??)
    #count_file = open(R1_fastq.replace('.fastq', '.counts'), 'w')
    count_file = open('%s/%s.counts' % (root, prefix), 'w')
    R2_fastq = R1_fastq.replace('R1', 'R2')

    # Determine the number of reads in both (R1 + R2) fastq files, store in the .count file
    stdout, stderr = subprocess.Popen(['wc', '-l', R1_fastq], stdout = subprocess.PIPE, stderr = subprocess.PIPE).communicate()
    total_reads_R1 = int(stdout.split()[0])/4
    stdout, stderr = subprocess.Popen(['wc', '-l', R2_fastq], stdout = subprocess.PIPE, stderr = subprocess.PIPE).communicate()
    total_reads_R2 = int(stdout.split()[0])/4
    logger.debug("{} R1 and {} R2 reads in {} and {}".format(total_reads_R1, total_reads_R2, R1_fastq, R2_fastq))
    count_file.write('Raw FASTQ,R1,%d,R2,%d\n' % (total_reads_R1, total_reads_R2))

    # Initial consensus B mapping
    prelim_samfile = '{}/{}.prelim.sam'.format(root, prefix)
    system_call('bowtie2 --quiet -p {} --local -x {} -1 {} -2 {} -S {}'.format(num_threads, refpath, R1_fastq, R2_fastq, prelim_samfile))

    # Define region-specific SAMs: refsams[refname] points to a nested dict which includes a file handle to each specific SAM
    refsams = {}
    for i, refname in enumerate(refnames):
        region_specific_sam = "{}/{}.{}.sam".format(root, prefix, refname)
        refsams.update({refname: {'sam_file_handle': open(region_specific_sam, 'w'),'count': [0,0]}})
    refsams.update({'*': {'sam_file_handle': open('%s/%s.unmapped.sam' % (root, prefix), 'w'),'count': [0,0]}})

    # Subdivide prelim SAMs into region-specific SAMs
    prelim_sam_infile = open(prelim_samfile, 'rU')
    line_counts = [0, 0]
    if filter_t_contaminants:
        t_counts = [0, 0]
        contam_file = open(R1_fastq.replace('.fastq', '.Tcontaminants.fastq'), 'w')

    for line in prelim_sam_infile:

        # Copy the original SAM header into each region-specific SAM
        if line.startswith('@'):
            for refname in refnames: refsams[refname]['sam_file_handle'].write(line)
            continue

        # SAM documentation explains what these fields mean - http://samtools.sourceforge.net/SAMv1.pdf
        qname, flag, refname, pos, mapq, cigar, rnext, pnext, tlen, seq, qual = line.strip('\n').split('\t')[:11]


        # TODO: remove soft clips that extend beyond boundaries of the reference
        if int(pos) == 0:
            matches = soft_clip.findall(cigar)
            if len(matches) > 0:
                clip_length = int(matches[0].rstrip('S'))

        if refname != '*' and int(pos) == refseq_lengths[refname]: # FIXME: adjust for zero-index?
            pass

        # SAM bitwise flag variable specifies whether read is paired, successfully mapped, etc
        bitinfo = samBitFlag(flag)
        if mode == 'Amplicon' and filter_t_contaminants:
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

        # FIXME: Does this identify R1 reads from R2 in the prelim SAM...?
        line_counts[bitinfo['is_reverse']] += 1

    prelim_sam_infile.close()
    if filter_t_contaminants:
        contam_file.close()

    # Show the number of reads that made it to preliminary mapping in the count file
    count_file.write('Preliminary map,{},{}\n'.format(line_counts[0], line_counts[1]))
    if filter_t_contaminants:
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
        if refname == '*':
            # don't attempt to remap reads from the unmapped SAM file
            continue

        # Ignore phiX, unmapped reads, and regions which had no mapping at the prelim mapping stage
        if sum(refsams[refname]['count']) == 0 or refname == 'phiX174' or refname == '*':
            continue

        # Run remap on the region-specific sam, and get the remapped sam and consensus pileup used to generate it
        samfile = refsams[refname]['sam_file_handle'].name
        logging.info("remap({},{},{},{},{},{})".format(R1_fastq, R2_fastq, samfile, refpath, original_reference, conseq_qCutoff))
        samfile, confile = remap(R1_fastq, R2_fastq, samfile, refpath, original_reference, conseq_qCutoff, num_threads)

        # Track file paths
        refsams[refname].update({'samfile': samfile, 'confile': confile})

        # Track the number of mapped reads in each region
        stdout, stderr = subprocess.Popen(['wc', '-l', samfile], stdout = subprocess.PIPE, stderr = subprocess.PIPE).communicate()
        region_specific_count = (int(stdout.split()[0]) - 3) / 2. # First 3 lines contain comments
        total_remap += region_specific_count
        refsams[refname]['count'][0] = region_specific_count
        count_file.write('remap %s,%d\n' % (refname, int(region_specific_count)))

    # Continue to remap if we've failed to map enough total reads
    if total_reads_R1 > 0:
        mapping_efficiency = total_remap / total_reads_R1 
    else:
        logger.info("Empty fastq: aborting mapping")
        count_file.close()
        return
    
    if mapping_efficiency < REMAP_THRESHOLD:
        break_out = False
        logger.info("Poor remapping efficiency for {} of {:.2%} - trying more iterations".format(sample_name, mapping_efficiency))

        # Repeat remapping up to the number of MAX_REMAPS permitted
        for iter in range(MAX_REMAPS):
            logger.debug("Additional remapping for {} (iteration #{})".format(sample_name, iter))

            total_remap = 0
            for refname in refnames:
                if refsams[refname]['count'][0] == 0 or refname == 'phiX174' or refname == '*':
                    continue

                samfile = refsams[refname]['samfile']
                confile = refsams[refname]['confile']
                samfile, confile = remap(R1_fastq, R2_fastq, samfile, confile, original_reference, conseq_qCutoff, num_threads)
                refsams[refname]['samfile'] = samfile
                refsams[refname]['confile'] = confile

                # Continue to determine the number of reads mapped in this region
                stdout, stderr = subprocess.Popen(['wc', '-l', samfile], stdout = subprocess.PIPE, stderr = subprocess.PIPE).communicate()
                region_specific_count = (int(stdout.split()[0]) - 3) / 2.


                # FIXME: Something is breaking here... coverage goes from 85% to 0%!
                if region_specific_count < refsams[refname]['count'][0]:
                    logger.warn("Remapping for {} region {} resulting in LESS reads ({:.2%} vs {:.2%}) - halting iterative remapping".format(
                        sample_name, refname, float(region_specific_count / total_reads_R1), float(refsams[refname]['count'][0] / total_reads_R1)))
                    break_out = True
                    break

                total_remap += region_specific_count
                refsams[refname]['count'][0] = region_specific_count

                # Write each remap iteration result to the file
                count_file.write('remap %d %s,%d\n' % (iter, refname, int(region_specific_count)))

            mapping_efficiency = float(total_remap / total_reads_R1)
            logger.info("Mapping efficiency for {} is now {:.2%}".format(sample_name, mapping_efficiency))
            if break_out or mapping_efficiency >= REMAP_THRESHOLD:
                break

    else:
        logger.info("{} had acceptable remap efficiency ({:.2%})".format(sample_name, mapping_efficiency))
    count_file.close()

def g2p_scoring(csf_path, g2p_alignment_cutoff):
    """
    Take an env (amplicon) CSF and generate a v3prot file.

    Header: contains the G2P FPR and read count
    Sequence: protein aligned V3

    The CSF must be from an amplicon run: column 1 must contain the rank + count.
    """

    import logging,os,sys
    from hyphyAlign import apply2nuc, change_settings, get_boundaries, HyPhy, pair_align, refSeqs
    from minG2P import conan_g2p
    from miseqUtils import translate_nuc

    csf_filename = os.path.basename(csf_path)
    prefix = csf_filename.split('.')[0]
    logger = logging.getLogger()
    hyphy = HyPhy._THyPhy (os.getcwd(), 1)			# HyPhy is used for alignment
    change_settings(hyphy)					# Configure scoring matrix / gap penalties
    refseq = translate_nuc(refSeqs['V3, clinical'], 0)	# V3 ref seq is NON-STANDARD: talk to Guin

    if csf_filename.find("HIV1B-env") == -1 or not csf_path.endswith('.csf'):
        return logger.error("{} is not an HIV1B-env CSF file".format(csf_filename))

    # Store CSF in fasta-like variable called sequences
    sequences = []
    with open(csf_path, 'rU') as csf_file:
        for line in csf_file:
            header, left_offset, seq_no_gaps  = line.strip("\n").split(",")
            sequences.append((header, seq_no_gaps))

    if len(sequences) == 0:
        # skip empty file
        return logger.error('%s is an empty file' % csf_filename)

    # Determine offset from 1st sequence to correct frameshift induced by sample-specific remapping
    seq1 = sequences[0][1].strip("-")
    best_offset = 0
    best_score = -999
    possible_ORFs = [0, 1, 2]
    for offset in possible_ORFs:
        aaEnvSeq = translate_nuc(seq1, offset)
        aquery, aref, ascore = pair_align(hyphy, refseq, aaEnvSeq)
        if ascore > best_score:
            best_offset = offset
            best_score = ascore

    # For each env sequence, extract the V3 nucleotide sequence
    badfile = open(csf_path.replace('.csf', '.badV3'), 'w')
    v3nucs = {}
    for header, seq in sequences:
        count = int(header.split('_')[-1])
        seq = seq.replace("-","")					# Strip dashes at flanking regions generated by alignment
        aaEnvSeq = translate_nuc(seq, best_offset)			# Translate env on correct ORF
        aquery, aref, ascore = pair_align(hyphy, refseq, aaEnvSeq)
        left, right = get_boundaries(aref)				# Get left/right boundaries of V3 protein
        v3prot = aquery[left:right]					# Extract V3 protein
        v3nuc = apply2nuc(seq[(3*left-best_offset):], v3prot,		# Use alignment to extract V3 nuc seq
                aref[left:right], keepIns=True, keepDel=False)

        # Drop V3 data that don't satisfy quality control
        if 'N' in v3nuc or not v3prot.startswith('C') or not v3prot.endswith('C') or '*' in v3prot or ascore < g2p_alignment_cutoff or len(v3prot) < 32 or len(v3prot) > 40:
            badfile.write('>%s_reason_%s\n%s\n' % (header,
                '|'.join(['stopcodon' if '*' in v3prot else '',					# V3 can't have internal stop codon
                'lowscore' if ascore < g2p_alignment_cutoff else '',				# The G2P alignment can't be poor
                'cystines' if not v3prot.startswith('C') or not v3prot.endswith('C') else '',	# V3 must start/end with C
                'ambig' if 'N' in v3nuc else '']),seq))						# There must be no unknown bases
        else:
            # Track the count of each v3 nucleotide sequence
            if v3nucs.has_key(v3nuc):
                v3nucs[v3nuc] += count
            else:
                v3nucs.update({v3nuc: count})
    badfile.close()

    # Calculate g2p scores for each v3 nuc sequence
    v3prots = {}
    for v3nuc, count in v3nucs.iteritems():
        g2p, fpr, aligned = conan_g2p(v3nuc)

        if g2p is None:
            continue

        # Track the count of each protein sequence
        if v3prots.has_key(aligned):
            v3prots[aligned]['count'] += count
        else:
            # Dict within dict - store count and fpr for each sequence
            v3prots.update({aligned: {'count': count, 'fpr': fpr}})

    # Collect v3 prot sequences and their output (v is a dict mapping to count and fpr)
    intermed = [(v['count'], v3prot) for v3prot, v in v3prots.iteritems()]
    intermed.sort(reverse=True)

    # For this sample, write a v3prot file containing the prefix, sequence, rank, count, and fpr
    v3prot_path = csf_path.replace('.csf', '.v3prot')
    logger.info("Writing results to {}".format(v3prot_path))
    with open(v3prot_path, 'w') as v3protfile:
        for i, (count, v3prot) in enumerate(intermed):
            fpr = v3prots[v3prot]['fpr']
            v3protfile.write(">{}_variant_{}_count_{}_fpr_{}\n{}\n".format(prefix, i, count, fpr, v3prot))


def sam2csf_with_base_censoring(samfile, censoring_qCutoff, mapping_cutoff, mode, max_prop_N):
    """
    From a path to a remapped SAM containing paired reads, create a comma-delimited,
    3-column CSF file with the merged, censored sequence.

    Column 1: Either the rank + count for amplicon, or qname for Nextera
    Column 2: Left-offset of the read
    Column 3: Gap-stripped sequence (Alignment data is lost)
    """
    import logging
    from miseqUtils import len_gap_prefix, sam2fasta

    filename = samfile.split('/')[-1]
    prefix, region = filename.split('.')[:2]

    csf_filename = '.'.join(map(str,[samfile.replace('.remap.sam', ''), censoring_qCutoff, 'csf']))

    with open(samfile, 'rU') as handle:
        fasta = sam2fasta(handle, censoring_qCutoff, mapping_cutoff, max_prop_N)

        if mode == 'Amplicon':
            # compress identical reads
            d = {}
            for qname, seq in fasta:
                if d.has_key(seq):
                    d[seq] += 1
                    continue
                d.update({seq: 1})

            if len(d) == 0:
                # empty SAM file or no qualified reads
                return logging.warn("{} likely empty or invalid - halting sam2fasta".format(samfile))

            # sort by read count in decreasing order
            intermed = [(count, len_gap_prefix(seq), seq) for seq, count in d.iteritems()]
            intermed.sort(reverse=True)

            # Write CSF to disk
            with open(csf_filename, 'w') as outfile:
                for rank, (count, left_offset, seq) in enumerate(intermed):
                    outfile.write('{}_{},{},{}\n'.format(rank, count, left_offset, seq.strip('-')))

        elif mode == 'Nextera':
            # Sort csf by left-gap prefix: the offset of the read relative to the ref seq
            intermed = [(len_gap_prefix(seq), qname, seq) for qname, seq in fasta]
            intermed.sort()
            with open(csf_filename, 'w') as outfile:
                for (left_offset, qname, seq) in intermed:
                    outfile.write("{},{},{}\n".format(qname, left_offset, seq.strip('-')))

        else:
            return logging.warn("ERROR: Unrecognized mode {} in miseqUtils.sam2csf()".format(mode))


def slice_outputs(root, region_slices):
    """
    Slice matrix output of previous step (*.nuc|amino.count.csv) into sub-regions.
    """

    import os, sys, time
    from glob import glob

    # Coordinates are in nucleotide space: start/end are inclusive (Relative to HXB2 aligned sequences)
    # Ex: region_slices = [("PROTEASE", "HIV1B-pol", 1, 297), ("V3", "HIV1B-env", 887, 993)]

    # For each region slice rule
    for rule in region_slices:
        slice, region, start, end = rule
        files = glob(root + '/*.{}.*nuc.csv'.format(region))
        files += glob(root + '/*.{}.*amino.csv'.format(region))

        # Get all nuc/amino files containing the region to be sliced
        for path in files:
            fileName = os.path.basename(path)
            sample,old_region = fileName.split(".")[:2]

            with open(path, 'rU') as f:
                lines = f.readlines()

            newFileName = fileName.replace(region,slice)
            dirName = os.path.dirname(path)
            slice_filename = "{}/{}".format(dirName, newFileName)
            f = open(slice_filename, 'w')
            conseq_filename = "{}/{}".format(dirName, newFileName.replace(".csv",".conseq"))
            f_conseq = open(conseq_filename, 'w')

            conseq = ""
            is_empty = True

            # For each line in the frequency matrix file
            for i,line in enumerate(lines):
                line = line.rstrip("\n")

                # First, extract character dictionary from header
                if (i == 0):
                    f.write("{}\n".format(line))
                    dictionary = line.split(",")[2:]
                    continue

                # For nuc.csv, hxb2_pos is in nucleotide space
                query_pos, hxb2_pos = map(int, line.split(",")[:2])
                if "nuc.csv" in path:
                    if hxb2_pos < start or hxb2_pos > end+1: continue
                    region_pos = hxb2_pos - start + 1

                # For amino.csv, hxb2_pos is in amino space
                elif "amino.csv" in path:
                    if hxb2_pos < (start+2)/3 or hxb2_pos > end/3: continue
                    region_pos = hxb2_pos - (start+2)/3 + 1

                # If we reached this point, the slice contains data
                is_empty = False

                # Generate consensus sequence
                counts = line.split(",")[2:]
                max_count = max(map(int,counts))
                max_char = filter(lambda x: int(x) == max_count, counts)
                index = counts.index(max_char[0])
                majority_char = dictionary[index]
                conseq += majority_char
                f.write("{},{},{}\n".format(query_pos,region_pos,",".join(counts)))

            # Close the new slice matrix, and write the consensus sequence of it
            f.close()
            f_conseq.write(">{}_{}\n{}".format(sample, slice, conseq))
            f_conseq.close()

            # If slice contains no data, delete it
            if is_empty:
                os.remove(slice_filename)
                os.remove(conseq_filename)


def filter_cross_contaminants (root, qcutoff, dbname = 'sample_references', num_threads=1):
    """
    Use sample-specific consensus sequences to make a bowtie2 database.
    Do this in a region-specific manner, otherwise contamination
    screening is meaningless.  Apply bowtie2 (without local alignment)
    to CSF outputs by converting the latter into FASTA and use the
    newly-generated database as a reference set.
    Report the numbers of cross-contaminants in *.counts file (as in mapping()).
    """
    from glob import glob
    from miseqUtils import parse_fasta, sampleSheetParser


    # parse sample sheet
    with open(root+'/SampleSheet.csv', 'rU') as sample_sheet:
        run_info = sampleSheetParser(sample_sheet)

    files = glob(root+'/*.%d.conseq' % qcutoff)
    assert len(files) > 0, 'No *.conseq files found with qcutoff (%f)' % qcutoff

    regions = set([f.split('/')[-1].split('.')[-3] for f in files])

    for region in regions:
        outfile = open('%s/%s.%s.fas' % (root, dbname, region), 'w')
        for f in files:
            if region not in f.split('/')[-1]:
                continue
            handle = open(f, 'rU')
            fasta = parse_fasta(handle)
            handle.close()

            for h, s in fasta:
                if h.endswith('_MAX') and len(s) > 0:
                    outfile.write('>%s\n%s\n' % (h.replace('_MAX', ''), s))

        outfile.close()

        # generate .bt2 files
        system_call('bowtie2-build -f -q %s/%s.%s.fas %s/%s.%s' % (root, dbname, region,
                                                                   root, dbname, region))

        # screen all CSF files
        csf_files = glob(root+'/*.%s.%d.csf' % (region, qcutoff))
        for csf_file in csf_files:
            filename = csf_file.split('/')[-1]
            #print filename
            sample_name = filename.split('.')[0] # includes S number

            # check whether this sample should be excluded from filtering
            if run_info['Data'][sample_name]['disable_contamination_check']:
                # simply copy *.csf to *.clean.csf
                system_call('cp %s %s' % (csf_file, csf_file.replace('.csf', '.clean.csf')))
                continue

            sam_file = csf_file.replace('.csf', '.filtering.sam')

            # export contents of CSF into FASTA file, storing counts in header
            fasta_file = csf_file.replace('.csf', '.csf.fa')
            output_handle = open(fasta_file, 'w')
            handle = open(csf_file, 'rU')
            for line in handle:
                label, offset, seq = line.strip('\n').split(',') # label contains count
                output_handle.write('>%s+%s\n%s\n' % (label, offset, seq))
            handle.close()
            output_handle.close()

            # note -f flag to indicate FASTA input
            system_call('bowtie2 --quiet -p %d -x %s/%s.%s -f -U %s -S %s' %
                        (num_threads, root, dbname, region, fasta_file, sam_file))

            # partition results into new CSF file and contaminant file
            clean_file = open(csf_file.replace('.csf', '.clean.csf'), 'w')
            contam_file = open(csf_file.replace('.csf', '.contam.csf'), 'w')
            handle = open(sam_file, 'rU')
            n_clean = 0
            n_contam = 0
            for line in handle:
                if line.startswith('@'):
                    continue
                qname, flag, refname, pos, mapq, cigar, rnext, pnext, tlen, seq, qual = line.strip('\n').split('\t')[:11]
                label, offset = qname.split('+')

                if refname == sample_name:
                    clean_file.write(','.join([label, offset, seq]))
                    clean_file.write('\n')
                    n_clean += 1
                else:
                    contam_file.write(','.join([label, offset, refname, mapq, seq]))
                    contam_file.write('\n')
                    n_contam += 1

            clean_file.close()
            contam_file.close()

            # append results to file generated by mapping()
            count_file = open('%s/%s.counts' % (root, sample_name), 'a')
            count_file.write('cross-contamination filter %s,%d,%d,%d\n' % (region, qcutoff, n_clean, n_contam))
            count_file.close()






