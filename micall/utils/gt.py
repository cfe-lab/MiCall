import sys
from subprocess import Popen
import shlex
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIXML
from tempfile import NamedTemporaryFile

def genotype(fasta,
#        db='hcv.fasta',
        db='/opt/hcv_geno/hcv.fasta',
        sc_t=200, disambig=False):

    with NamedTemporaryFile() as tmp_out:
        xml_out = tmp_out.name
        cline = NcbiblastnCommandline(query=fasta,db=db,outfmt=5,out=xml_out,evalue=0.0001,gapopen=5,gapextend=2,penalty=-3,reward=1,max_target_seqs=5000)
        blastn_proc = Popen(shlex.split(str(cline)))
        if blastn_proc.wait():
            raise Exception
        handle = open(xml_out, 'r')
        samples = {} 
        for i, record in enumerate(NCBIXML.parse(handle)):
            top_scores = list(reversed(sorted([(d.title,d.score) for d in
                record.descriptions], key=lambda x: x[1])))
            
            gts_s = None
            if len(top_scores) == 0:
                gts_s = "not found" 
            else:
                top_score = top_scores[0][1]
                gts_s = top_scores[0][0]

            samples[i] = gts_s 
    return samples 
