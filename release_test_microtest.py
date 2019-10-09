from argparse import ArgumentParser
from csv import DictReader
from operator import itemgetter
from pathlib import Path


def read_file(sample_name, file_name):
    file_path = (Path(__file__).parent /
                 'micall/tests/microtest/Results/basespace' / file_name)
    with file_path.open() as f:
        yield from (row
                    for row in DictReader(f)
                    if row['sample'] == sample_name)


def check_v3loop(sample_name, expected_counts):
    nuc_rows = list(read_file(sample_name, 'nuc.csv'))
    assert nuc_rows
    for row in nuc_rows:
        pos = int(row['refseq.nuc.pos'])
        assert row['region'] == 'V3LOOP', pos
    count_rows = list(map(itemgetter('refseq.nuc.pos', 'A', 'C', 'G', 'T'),
                          nuc_rows))
    expected_count_rows = [tuple(line.split())
                           for line in expected_counts.splitlines()]
    assert len(count_rows) >= len(expected_count_rows), len(count_rows)
    for line, expected_line in zip(count_rows, expected_count_rows):
        assert line == expected_line, (line, expected_line)
    return nuc_rows


def check_1234():
    expected_counts = """\
1 0 0 0 10
2 0 0 10 0
3 0 10 0 0
4 10 0 0 0
5 0  9 0 1
6 10 0 0 0
7 10 0 0 0
8 0 0 10 0
9 10 0 0 0
10 0 10 0 0
"""
    check_v3loop('1234A-V3LOOP_S1', expected_counts)


def check_2000():
    expected_counts = """\
1 0 0 0 0
2 0 0 0 0
3 0 0 0 0
4 10 0 0 0
5 0 10 0 0
6 10 0 0 0
7 10 0 0 0
8 0 0 10 0
9 10 0 0 0
10 0 10 0 0
"""
    check_v3loop('2000A-V3LOOP_S2', expected_counts)


def check_2010():
    nuc_rows = check_v3loop('2010A-V3LOOP_S3', '')
    assert nuc_rows
    for row in nuc_rows:
        pos = int(row['refseq.nuc.pos'])
        in_gap = 51 < pos <= 57
        assert row['coverage'] == '10' or in_gap, pos


def check_2020():
    nuc_rows = list(read_file('2020A-GP41_S4', 'nuc.csv'))
    assert nuc_rows
    expected_seeds = ('HIV1-B-KR-KJ140263-seed', 'HIV1-B-FR-KF716496-seed')
    for row in nuc_rows:
        pos = int(row['refseq.nuc.pos'])
        assert row['seed'] in expected_seeds, pos
        assert row['region'] == 'GP41', pos
        if 3 < pos <= 81:
            assert row['coverage'] == '10', pos
        else:
            assert row['coverage'] == '0', pos
        if pos == 30:
            assert row['ins'] == '10', pos


def check_2030():
    nuc_rows = list(read_file('2030A-V3LOOP_S5', 'nuc.csv'))
    assert not nuc_rows


def check_2040():
    conseq_rows = list(read_file('2040A-HLA-B_S6', 'conseq.csv'))
    num_rows = len(conseq_rows)
    assert num_rows == 7, num_rows
    max_row = conseq_rows[0]
    cutoff = max_row['consensus-percent-cutoff']
    assert cutoff == 'MAX', cutoff
    conseq = max_row['sequence']

    # Check for mixture
    assert conseq.startswith('GCTCCCWC'), conseq


def check_2050():
    nuc_rows = list(read_file('2050A-V3LOOP_S7', 'nuc.csv'))
    assert not nuc_rows


def check_2060():
    summary_rows = list(read_file('2060A-V3LOOP_S8', 'g2p_summary.csv'))
    assert len(summary_rows) == 1, summary_rows
    row = summary_rows[0]
    assert row['valid'] == '10', row['valid']
    assert row['X4calls'] == '0', row['X4calls']


def check_2070():
    nuc_rows = list(read_file('2070A-PR_S9', 'nuc.csv'))
    assert nuc_rows
    for row in nuc_rows:
        assert row['seed'] == 'HIV1-B-FR-K03455-seed', row['seed']
        assert row['region'] == 'PR', row['region']
        pos = int(row['refseq.nuc.pos'])
        if 118 <= pos <= 207:
            assert row['coverage'] == '13', pos
            if 133 <= pos <= 135:
                assert row['del'] == '13', pos
            elif 193 <= pos <= 195:
                assert row['del'] == '3', pos
        else:
            assert row['coverage'] == '0', pos


def check_2080():
    check_v3loop('2080A-V3LOOP_S10', '')


def check_2100(is_denovo):
    conseq_rows = list(read_file('2100A-HCV-1337B-V3LOOP-PWND-HIV_S12',
                                 'conseq.csv'))
    regions = set(map(itemgetter('region'), conseq_rows))
    if is_denovo:
        expected_regions = {'2_3-HIV1-B-FR-K03455-seed', '1-HCV-1a'}
    else:
        expected_regions = {'HIV1-CON-XX-Consensus-seed',
                            'HCV-1a',
                            'HIV1-B-FR-K03455-seed'}
    assert regions == expected_regions, regions


def check_amino_row(row,
                    coverage=0,
                    stop_codons=0,
                    low_quality=0,
                    partial=0,
                    dels=0,
                    ins=0,
                    clip=0,
                    v3_overlap=0):
    pos = int(row['refseq.aa.pos'])
    assert row['coverage'] == str(coverage), f'{row["coverage"]} at {pos}'
    assert row['*'] == str(stop_codons), f'{row["*"]} at {pos}'
    assert row['X'] == str(low_quality), f'{row["X"]} at {pos}'
    assert row['partial'] == str(partial), f'{row["partial"]} at {pos}'
    assert row['del'] == str(dels), f'{row["del"]} at {pos}'
    assert row['ins'] == str(ins), f'{row["ins"]} at {pos}'
    assert row['clip'] == str(clip), f'{row["clip"]} at {pos}'
    assert row['v3_overlap'] == str(v3_overlap), f'{row["v3_overlap"]} at {pos}'


def check_2110(is_denovo):
    amino_rows = list(read_file('2110A-V3LOOP_S13', 'amino.csv'))
    assert amino_rows
    for row in amino_rows:
        assert row['region'] == 'V3LOOP', row['region']
        pos = int(row['refseq.aa.pos'])
        if pos == 6:
            check_amino_row(row, coverage=13, stop_codons=1)
        elif pos == 7:
            check_amino_row(row, coverage=11, low_quality=2)
        elif pos == 8:
            check_amino_row(row, coverage=10, partial=3)
        elif pos == 9:
            # G2P doesn't currently handle inserts.
            expected_inserts = 4 if is_denovo else 0
            check_amino_row(row, coverage=13, ins=expected_inserts)
        elif pos == 10:
            check_amino_row(row, coverage=13, dels=4)
        elif pos < 18:
            check_amino_row(row, coverage=13)
        elif pos < 33:
            check_amino_row(row)
        else:
            if is_denovo:
                check_amino_row(row, coverage=10)
            else:
                check_amino_row(row, v3_overlap=10)


def check_2120():
    amino_rows = list(read_file('2120A-PR_S14', 'amino.csv'))
    assert amino_rows
    for row in amino_rows:
        pos = int(row['refseq.aa.pos'])
        if row['region'] == 'PR':
            if pos <= 29:
                check_amino_row(row)
            elif pos <= 46:
                check_amino_row(row, coverage=13)
            elif pos <= 49:
                check_amino_row(row, coverage=13, stop_codons=1)
            elif pos <= 52:
                check_amino_row(row, coverage=11, low_quality=2)
            elif pos <= 54:
                check_amino_row(row, coverage=10, partial=3)
            elif pos <= 55:
                check_amino_row(row, coverage=13)
            elif pos <= 56:
                check_amino_row(row, coverage=13, ins=4)
            elif pos <= 58:
                check_amino_row(row, coverage=13)
            elif pos <= 64:
                check_amino_row(row, coverage=10, clip=3)
            elif pos <= 67:
                check_amino_row(row, coverage=10, dels=1)
            elif pos <= 90:
                check_amino_row(row, coverage=10)
            else:
                check_amino_row(row)
        else:
            assert row['region'] == 'RT', row['region']
            if pos <= 33:
                check_amino_row(row, coverage=10)
            else:
                check_amino_row(row)


def check_2130(is_denovo):
    conseq_rows = list(read_file('2130A-HCV_S15', 'conseq.csv'))
    regions = set(map(itemgetter('region'), conseq_rows))
    expected_regions = {'1_2-HCV-2a'} if is_denovo else {'HCV-2a'}
    assert regions == expected_regions, regions


def check_2130midi(is_denovo):
    conseq_rows = list(read_file('2130AMIDI-MidHCV_S16', 'conseq.csv'))
    regions = set(map(itemgetter('region'), conseq_rows))
    expected_regions = {'1_2-HCV-2a'} if is_denovo else {'HCV-2a'}
    assert regions == expected_regions, regions


def check_2140():
    resistance_rows = list(read_file('2140A-HIV_S17', 'resistance.csv'))
    pr_rows = [row for row in resistance_rows if row['region'] == 'PR']
    assert pr_rows
    for row in pr_rows:
        assert row['level'] != '0', (row['drug_name'], row['level_name'])
        if row['drug'] == 'IDV/r':
            assert row['level'] == '3', row['level']


def check_2160():
    amino_rows = list(read_file('2160A-HCV_S19', 'amino.csv'))
    assert amino_rows
    for row in amino_rows:
        assert row['region'] == 'HCV2-JFH-1-NS5b', row['region']
        pos = int(row['refseq.aa.pos'])
        coverage = int(row['coverage'])
        coverage_message = f'{coverage} coverage at {pos}'
        if pos < 40:
            assert coverage < 10, coverage_message
        elif 80 <= pos < 160:
            assert 10 < coverage, coverage_message
        elif 250 <= pos:
            assert coverage < 10, coverage_message


def check_2160midi():
    amino_rows = list(read_file('2160AMIDI-MidHCV_S20', 'amino.csv'))
    assert amino_rows
    for row in amino_rows:
        assert row['region'] == 'HCV2-JFH-1-NS5b', row['region']
        pos = int(row['refseq.aa.pos'])
        coverage = int(row['coverage'])
        coverage_message = f'{coverage} coverage at {pos}'
        if pos < 350:
            assert coverage < 10, coverage_message
        elif 420 <= pos < 520:
            assert 10 < coverage, coverage_message
        elif 580 <= pos:
            assert coverage < 10, coverage_message


def check_2170():
    amino_rows = list(read_file('2170A-HCV_S21', 'amino.csv'))
    assert amino_rows
    for row in amino_rows:
        pos = int(row['refseq.aa.pos'])
        coverage = int(row['coverage'])
        coverage_message = f'{coverage} coverage at {pos}'
        if row['region'] == 'HCV1A-H77-NS5a':
            if pos < 20:
                assert coverage < 10, coverage_message
            elif 50 <= pos:
                assert 10 < coverage, coverage_message
        elif row['region'] == 'HCV1A-H77-NS5b':
            assert 10 < coverage, coverage_message
        elif row['region'] == 'HCV2-JFH-1-NS5a':
            if pos < 15:
                assert coverage < 10, coverage_message
            elif 80 <= pos:
                assert 10 < coverage, coverage_message
        else:
            assert row['region'] == 'HCV2-JFH-1-NS5b', row['region']
            assert 10 < coverage, coverage_message


def main():
    parser = ArgumentParser(description='Check results for microtest samples.')
    parser.add_argument('--denovo',
                        action='store_true',
                        help='Samples were processed with de novo assembly.')
    args = parser.parse_args()
    is_denovo = args.denovo

    check_1234()
    check_2000()
    check_2010()
    check_2020()
    check_2030()
    check_2040()
    check_2050()
    check_2060()
    check_2070()
    check_2080()
    check_2100(is_denovo)
    check_2110(is_denovo)
    check_2120()
    check_2130(is_denovo)
    check_2130midi(is_denovo)
    check_2140()
    check_2160()
    check_2160midi()
    check_2170()
    print('Passed.')


main()
