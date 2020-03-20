import os
from micall.utils import primer_finder

# def test_hivseqinr():
#     hiv = primer_finder.Hivseqinr()
#     content = hiv.download('.\lol')
#     if os.path.isfile('./lol/R_HXB2.fasta'):
#         os.rmdir('./lol')
#         return True
#     return False

def test_hxb2_db_creation():
    hiv = primer_finder.Hivseqinr('C:/Users/Dan/Documents/CFE/projects/micall/repo/MiCall/lol', '../../hivseqinr/MiCallSummaryDaniel/190501_M05995/190501_M05995.fasta')