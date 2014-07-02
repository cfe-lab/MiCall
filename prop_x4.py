def prop_x4 (v3prot_path, fpr_cutoff, min_count):
    """Determine proportion X4 from v3prot file"""

    import logging

    with open(v3prot_path, 'rU') as infile:
        try:
            logging.info("Reading {}".format(v3prot_path))
            fasta = convert_fasta(infile.readlines())
        except Exception as e:
            raise Exception("convert_fasta threw exception '{}'".format(str(e)))

    total_count = 0
    total_x4_count = 0

    # Example header: F00309-IL_variant_0_count_27_fpr_4.0
    for h, _ in fasta:
        try:
            tokens = h.split('_')
            count = int(tokens[tokens.index('count')+1])
            fpr = float(tokens[tokens.index('fpr')+1])
        except:
            continue

        if count < min_count:
            continue
        if fpr <= fpr_cutoff:
            total_x4_count += count
        total_count += count

    proportion_X4 = (float(total_x4_count) / float(total_count)) * 100
    return (proportion_X4, total_x4_count, total_count)

def convert_fasta (lines):    
    blocks = []
    sequence = ''
    h = None
    for i in lines:
        if i[0] == '$': # skip h info
            continue
        elif i[0] == '>' or i[0] == '#':
            if len(sequence) > 0:
                blocks.append([h,sequence])
                sequence = ''    # reset containers
                h = i.strip('\n')[1:]
            else:
                h = i.strip('\n')[1:]
        else:
            sequence += i.strip('\n')
    try:
        blocks.append([h,sequence])    # handle last entry
    except:
        raise Exception("convert_fasta(): Error appending to blocks [{},{}]".format(h, sequence))
    return blocks
