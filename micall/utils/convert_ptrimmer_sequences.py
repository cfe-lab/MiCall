from _csv import writer
from collections import defaultdict
from contextlib import contextmanager
from csv import DictReader, reader
from itertools import groupby
from operator import itemgetter
from pathlib import Path

import requests
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from micall.utils.fetch_sequences import fetch_by_accession


def main():
    source_path = Path(__file__).parent.parent / 'core'
    data_path = Path(__file__).parent / 'ptrimmer_data'
    data_path.mkdir(exist_ok=True)

    # {group_num: {direction: [(row, definition_row)]}}
    primers = defaultdict(lambda: defaultdict(list))

    ref_seq = None
    with load_artic_rows(data_path) as rows, \
            load_definition_rows(data_path) as definition_rows:
        for row, definition_row in zip(rows, definition_rows):
            name = row['name']
            assert name == definition_row[3]

            name_fields = name.split('_')
            group_num = int(name_fields[1])
            direction = name_fields[2]
            primers[group_num][direction].append((row, definition_row))
    with (source_path/'ptrimmer_primers.txt').open('w') as out_file:
        row_writer = writer(out_file, delimiter='\t')
        for group_num, group in primers.items():
            left_pairs = group['LEFT']
            right_pairs = group['RIGHT']
            left_pairs.reverse()
            right_pairs.reverse()
            for left_row, left_definition_row in left_pairs:
                for right_row, right_definition_row in right_pairs:
                    if ref_seq is None:
                        accession = left_definition_row[0]
                        ref_seq = load_accession(accession, data_path)

                    left_start = int(left_definition_row[1])
                    left_end = int(left_definition_row[2])
                    left_primer_seq = Seq(ref_seq[left_start:left_end])
                    assert left_primer_seq == left_row['seq'], (left_primer_seq,
                                                                left_row['seq'])

                    right_start = int(right_definition_row[1])
                    right_end = int(right_definition_row[2])
                    right_primer_seq = Seq(ref_seq[right_start:right_end])
                    complement = right_primer_seq.reverse_complement()
                    assert complement == right_row['seq'], (complement,
                                                            right_row['seq'])

                    left_fields = left_row['name'].split('_')
                    right_fields = right_row['name'].split('_')
                    name = '_'.join(left_fields[:2])
                    if max(len(fields) for fields in (left_fields,
                                                      right_fields)) > 3:
                        for fields in (left_fields, right_fields):
                            if len(fields) == 3:
                                name += '_main'
                            else:
                                name += '_' + fields[3]
                    new_row = [left_primer_seq,
                               right_primer_seq,
                               right_start - left_end,
                               name]
                    row_writer.writerow(new_row)


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


@contextmanager
def load_definition_rows(data_path: Path):
    with load_artic_file('nCoV-2019.scheme.bed', data_path) as source_file:
        yield reader(source_file, delimiter='\t')


def load_accession(accession: str, data_path: Path):
    ref_file: Path = data_path / (accession + '.fasta')
    if ref_file.exists():
        ref_seq = str(SeqIO.read(ref_file, 'fasta').seq)
    else:
        ref_seq = fetch_by_accession(accession)
        ref_file.write_text(f'>{accession}\n{ref_seq}')
    return ref_seq


def build_row(direction, row):
    assert direction in row['name'], row['name']
    new_seq = SeqRecord(Seq(row['seq']),
                        row['name'],
                        description='')
    return new_seq


main()
