import typing
from contextlib import contextmanager
from csv import DictReader
from difflib import Differ
from io import StringIO
from pathlib import Path

import requests
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from micall.core.project_config import ProjectConfig


def main():
    source_path = Path(__file__).parent / 'ptrimmer_data'
    target_path = source_path.parent.parent / 'data'
    seqs = {name: [] for name in ('left', 'right')}
    load_hivb(seqs)
    for name, seq_list in seqs.items():
        start_seqs = [SeqRecord('X' + seq.seq,
                                id=seq.id.replace(' ', '_'),
                                description=seq.description)
                      for seq in seq_list]
        end_seqs = [SeqRecord(seq.reverse_complement().seq + 'X',
                              id=seq.id.replace(' ', '_'),
                              description=seq.description)
                    for seq in seq_list]
        SeqIO.write(start_seqs, str(target_path/f'primers_hivb_{name}.fasta'), 'fasta')
        SeqIO.write(end_seqs, str(target_path/f'primers_hivb_{name}_end.fasta'), 'fasta')


def load_hivb(seqs: typing.Dict[str, typing.List[SeqRecord]]):
    left_source_path = Path(__file__).parent / 'MiCall project HIVB ForwardPrimers.txt'
    right_source_path = Path(__file__).parent / 'MiCall project HIVB ReversePrimers.txt'
    seqs['left'].extend(SeqIO.parse(left_source_path, 'fasta'))
    seqs['right'].extend(SeqIO.parse(right_source_path, 'fasta'))


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
    differ = Differ()
    for row in hcv_definitions:
        name = 'HCV ' + row['name']
        start, end = (int(pos) for pos in row['h77_pos'].split('-'))
        primer = SeqRecord(Seq(row['sequence']), name, description='')
        complement = primer.reverse_complement(id=primer.id, description='')
        direction = row['direction']
        if direction == 'F':
            seqs['left'].append(primer)
        else:
            seqs['right'].append(primer)
            primer, complement = complement, primer
        h77_section = Seq(h77[start - 1:end])
        if is_comparing and primer.seq != h77_section:
            print(name, 'does not match.')
            diffs = differ.compare([str(primer.seq) + '\n'],
                                   [str(h77_section) + '\n'])
            print(*diffs, sep='')


def load_corona(source_path, seqs):
    with load_artic_rows(source_path) as rows:
        for row in rows:
            new_seq = SeqRecord(Seq(row['seq']), row['name'], description='')
            if 'LEFT' in row['name']:
                fwd_seqs = seqs['left']
            else:
                assert 'RIGHT' in row['name'], row['name']
                fwd_seqs = seqs['right']
            fwd_seqs.append(new_seq)


@contextmanager
def load_artic_file(filename: str, data_path: Path):
    local_path = data_path / filename
    if not local_path.exists():
        url = 'https://github.com/artic-network/artic-ncov2019/raw/master/' \
              'primer_schemes/nCoV-2019/V3/' + filename
        response = requests.get(url)
        response.raise_for_status()
        with local_path.open('wb') as f:
            for chunk in response.iter_content(chunk_size=1000):
                f.write(chunk)
    with local_path.open() as source_file:
        yield source_file


@contextmanager
def load_artic_rows(data_path: Path):
    with load_artic_file('nCoV-2019.tsv', data_path) as source_file:
        yield DictReader(source_file, delimiter='\t')


main()
