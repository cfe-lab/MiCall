"""
Generate gene-by-gene bowtie reference for specific HCV genotypes
"""

import os
import HyPhy
import hyphyAlign

from project_config import ProjectConfig

def convert_fasta (lines):
    blocks = []
    h = None
    sequence = ''
    for i in lines:
        if i[0] == '$': # skip h info
            continue
        elif i[0] == '>' or i[0] == '#':
            if len(sequence) > 0:
                blocks.append([h, sequence])
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

hyphy = HyPhy._THyPhy(os.getcwd(), 1)  # @UndefinedVariable

dump = hyphy.ExecuteBF('MESSAGE_LOGGING = 0;', False)

hyphyAlign.change_settings(hyphy,
                            alphabet = hyphyAlign.nucAlphabet,
                            scoreMatrix = hyphyAlign.nucScoreMatrix,
                            gapOpen = 20,
                            gapOpen2 = 20,
                            gapExtend = 10,
                            gapExtend2 = 10,
                            noTerminalPenalty = 1)


with open('HCV_REF_2012_genome.fasta', 'rU') as handle:
    genomes = convert_fasta(handle)  # keep one per genotype

projects = ProjectConfig.loadDefault()
genes = ['E1', 'E2', 'NS2', 'NS3', 'NS4a', 'NS4b', 'NS5a', 'NS5b', 'core', 'p7']
h77 = dict([(gene, projects.getReference("HCV1A-H77-{}-seed".format(gene)))
            for gene in genes])

with open('hcv_genes.fasta', 'w') as outfile:
    processed_subtypes = set()
    for subtype, genome in genomes:
        if subtype in processed_subtypes:
            continue
        for region, refseq in h77.iteritems():
            print subtype, region
            aquery, aref, ascore = hyphyAlign.pair_align(hyphy, refseq, genome)
            left, right = hyphyAlign.get_boundaries(aref)
            outfile.write('>%s-%s\n%s\n' % (subtype,
                                            region,
                                            aquery[left:right].replace('-', '')))
        processed_subtypes.add(subtype)
