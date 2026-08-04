"""
Use project file to punch out genes from FDA amino acid refs
"""
import os
import HyPhy
import hyphyAlign
import json
from micall.utils.externals import ProjectsFile
from seqUtils import convert_fasta

hyphy = HyPhy._THyPhy (os.getcwd(), 1) # instance of HyPhy
hyphyAlign.change_settings(hyphy)  # default settings

handle = open('fda_hcv_polyprotein.fa', 'r')
fasta = convert_fasta(handle)
handle.close()

with ProjectsFile().path() as projects_file_path:
     handle = open(projects_file_path, 'r')
     proj = json.load(handle)
     handle.close()

h77 = {}
for key in proj['regions'].iterkeys():
    if 'H77' in key and not key.endswith('seed'):
        aa = ''.join(proj['regions'][key]['reference'])
        h77.update({str(key): str(aa)})
        
outfile = open('fda_hcv_coords.fa', 'w')

for h, s in fasta:
    for gene, refseq in h77.iteritems():
        aquery, aref, ascore = hyphyAlign.pair_align(hyphy, refseq, s)
        left, right = hyphyAlign.get_boundaries(aref)
        s2 = aquery[left:right].replace('-', '')
        outfile.write('>%s %s\n%s\n' % (gene.split('-')[-1], h, s2))

outfile.close()

        
