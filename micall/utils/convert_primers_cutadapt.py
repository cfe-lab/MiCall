from csv import DictReader
from io import StringIO
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from micall.core.project_config import ProjectConfig


def main():
    source_path = Path(__file__).parent / 'ptrimmer_data'
    target_path = source_path.parent.parent / 'core'
    seqs = {name: [] for name in ('left', 'right')}
    load_hcv(seqs)
    load_corona(source_path, seqs)
    for name, seq_list in seqs.items():
        if name == 'left':
            fwd_seqs = [SeqRecord('X' + seq.seq,
                                  id=seq.id,
                                  description=seq.description)
                        for seq in seq_list]
            rev_seqs = [SeqRecord(seq.reverse_complement().seq + 'X',
                                  id=seq.id,
                                  description=seq.description)
                        for seq in seq_list]
        else:
            fwd_seqs = [SeqRecord(seq.seq + 'X',
                                  id=seq.id,
                                  description=seq.description)
                        for seq in seq_list]
            rev_seqs = [SeqRecord('X' + seq.reverse_complement().seq,
                                  id=seq.id,
                                  description=seq.description)
                        for seq in seq_list]
        SeqIO.write(fwd_seqs, str(target_path/f'primers_{name}.fasta'), 'fasta')
        SeqIO.write(rev_seqs, str(target_path/f'primers_{name}_rev.fasta'), 'fasta')


def load_hcv(seqs):
    hcv_definitions = DictReader(StringIO("""\
protocol,name,direction,length,h77_pos,sequence
HCV WG,oligo dA20,R,20,9418-9437,AAAAAAAAAAAAAAAAAAAA
,Pr3,R,30,8616-8645,GGCGGAATTCCTGGTCATAGCCTCCGTGAA
,1abGENF1bp,F,28,266-293,GGGTCGCGAAAGGCCTTGTGGTACTGCC
,TIM-Pr3,R,30,8616-8645,CAGGAAACAGCTATGACGGCGGAATTCCTGGTCATAGCCTCCGTGAA
,1abGENF2,F,30,286-315,GTACTGCCTGATAGGGTGCTTGCGAGTGCC
,Pr6,R,30,8611-8640,AATTCCTGGTCATAGCCTCCGTGAAGACTC
HCV miDi,Pr1,F,31,8245-8275,TGGGGTTCGCGTATGATACCCGCTGCTTTGA
,Pr2,F,31,8245-8275,TGGGGTTTTCTTACGACACCAGGTGCTTTGA
,oligo dA20-TIM,R,20,9418-9437,CAGGAAACAGCTATGACAAAAAAAAAAAAAAAAAAAA
,Pr4,F,29,8253-8281,CCGTATGATACCCGCTGCTTTGACTCAAC
,Pr5,F,29,8253-8281,TCCTACGACACCAGGTGCTTTGATTCAAC
,TIM,R,,1-0,CAGGAAACAGCTATGAC
"""))
    projects = ProjectConfig.loadDefault()
    h77 = projects.getReference('HCV-1a')
    is_comparing = True
    for row in hcv_definitions:
        name = 'HCV ' + row['name']
        start, end = (int(pos) for pos in row['h77_pos'].split('-'))
        primer = SeqRecord(Seq(row['sequence']), name, description='')
        complement = primer.reverse_complement(id=primer.id, description='')
        direction = row['direction']
        if direction == 'F':
            seqs['left'].append(primer)
        else:
            primer, complement = complement, primer
            seqs['right'].append(primer)
        h77_section = Seq(h77[start - 1:end])
        if is_comparing and primer.seq != h77_section:
            print(name, 'does not match.')
            print(primer.seq)
            print(h77_section)


def load_corona(source_path, seqs):
    with (source_path / 'nCoV-2019.tsv').open() as f:
        rows = DictReader(f, delimiter='\t')
        for row in rows:
            new_seq = SeqRecord(Seq(row['seq']), row['name'], description='')
            if 'LEFT' in row['name']:
                fwd_seqs = seqs['left']
            else:
                assert 'RIGHT' in row['name'], row['name']
                fwd_seqs = seqs['right']
            fwd_seqs.append(new_seq)


main()
