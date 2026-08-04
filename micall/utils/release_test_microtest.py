import gzip
import inspect
import os
import shutil
import typing
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
from concurrent.futures.process import ProcessPoolExecutor
from csv import DictReader, DictWriter
from operator import itemgetter
from pathlib import Path
from subprocess import run, PIPE, CalledProcessError
from tempfile import mkdtemp
from traceback import print_exc

import sys

from micall.monitor.find_groups import find_groups, SampleGroup


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


class ResultsFolder:
    def __init__(self, path: Path, is_denovo: bool):
        self.path = path
        self.is_denovo = is_denovo

    def read_file(self, sample_name, file_name) -> typing.Iterator[dict]:
        scratch_path = self.path / 'scratch' / sample_name
        if self.is_denovo:
            folder_names = ('output_denovo', 'output_resistance_denovo')
        else:
            folder_names = ('output', 'output_resistance')
        for folder_name in folder_names:
            file_path = scratch_path / folder_name / file_name
            if file_path.exists():
                with file_path.open() as f:
                    yield from DictReader(f)
                break
        else:
            raise FileNotFoundError(f'{sample_name}/{file_name}')

    def check_v3loop(self, sample_name: str, expected_counts: str):
        nuc_rows = list(self.read_file(sample_name, 'nuc.csv'))
        assert nuc_rows
        v3_loop_rows = []
        for row in nuc_rows:
            if row['region'] == 'V3LOOP':
                v3_loop_rows.append(row)
            else:
                pos = int(row['refseq.nuc.pos'])
                assert row['region'] == 'GP120', (pos, row['region'])
        count_rows = list(map(itemgetter('refseq.nuc.pos', 'A', 'C', 'G', 'T'),
                              v3_loop_rows))
        expected_count_rows = [tuple(line.split())
                               for line in expected_counts.splitlines()]
        assert len(count_rows) >= len(expected_count_rows), len(count_rows)
        for line, expected_line in zip(count_rows, expected_count_rows):
            assert line == expected_line, (line, expected_line)
        return v3_loop_rows

    def check_1234(self):
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
        self.check_v3loop('1234A-V3LOOP_S1', expected_counts)

    def check_2000(self):
        expected_counts = """\
4 10 0 0 0
5 0 10 0 0
6 10 0 0 0
7 10 0 0 0
8 0 0 10 0
9 10 0 0 0
10 0 10 0 0
"""
        self.check_v3loop('2000A-V3LOOP_S2', expected_counts)

    def check_2010(self):
        nuc_rows = self.check_v3loop('2010A-V3LOOP_S3', '')
        assert nuc_rows
        for row in nuc_rows:
            pos = int(row['refseq.nuc.pos'])
            in_gap = 51 < pos <= 57
            assert row['coverage'] == '10' or in_gap, pos

    def check_2020(self):
        nuc_rows = list(self.read_file('2020A-GP41_S4', 'nuc.csv'))
        assert nuc_rows
        expected_seeds = ('HIV1-B-KR-KJ140263-seed', 'HIV1-B-FR-KF716496-seed')
        for row in nuc_rows:
            pos = int(row['refseq.nuc.pos'])
            assert row['seed'] in expected_seeds, pos
            assert row['region'] == 'GP41', pos
            if 3 < pos <= 111:
                assert row['coverage'] == '10', pos
            else:
                assert row['coverage'] == '0', pos
            if pos == 30:
                assert row['ins'] == '10', pos

    def check_2030(self):
        nuc_rows = list(self.read_file('2030A-V3LOOP_S5', 'nuc.csv'))
        assert not nuc_rows

    def check_2040(self):
        conseq_rows = list(self.read_file('2040A-HLA-B_S6', 'conseq.csv'))
        num_rows = len(conseq_rows)
        assert num_rows == 7, num_rows
        max_row = conseq_rows[5]
        cutoff = float(max_row['consensus-percent-cutoff'])
        assert cutoff == 0.2, cutoff
        conseq = max_row['sequence']

        # Check for mixture
        assert conseq.startswith('GCTCCCWC'), conseq

    def check_2050(self):
        nuc_rows = list(self.read_file('2050A-V3LOOP_S7', 'nuc.csv'))
        assert not nuc_rows

    def check_2060(self):
        summary_rows = list(self.read_file('2060A-V3LOOP_S8', 'g2p_summary.csv'))
        assert len(summary_rows) == 1, summary_rows
        row = summary_rows[0]
        assert row['valid'] == '10', row['valid']
        assert row['X4calls'] == '0', row['X4calls']

    def check_2070(self):
        nuc_rows = list(self.read_file('2070A-PR_S9', 'nuc.csv'))
        assert nuc_rows
        for row in nuc_rows:
            assert row['seed'] == 'HIV1-B-FR-K03455-seed', row['seed']
            assert row['region'] == 'PR', row['region']
            pos = int(row['refseq.nuc.pos'])
            if 118 <= pos <= 240:
                assert row['coverage'] == '15', pos
                if 133 <= pos <= 135:
                    assert row['del'] == '15', pos
                elif 193 <= pos <= 195:
                    assert row['del'] == '3', pos
            else:
                assert row['coverage'] == '0', pos

    def check_2080(self):
        self.check_v3loop('2080A-V3LOOP_S10', '')

    def check_2100(self):
        conseq_rows = list(self.read_file('2100A-HCV-1337B-V3LOOP-PWND-HIV_S12',
                                          'conseq.csv'))
        regions = set(map(itemgetter('region'), conseq_rows))
        if self.is_denovo:
            expected_regions = {'HIV1-CON-XX-Consensus-seed',
                                '1-HIV1-B-FR-K03455-seed',
                                '2-HCV-1a'}
        else:
            expected_regions = {'HIV1-CON-XX-Consensus-seed',
                                'HCV-1a',
                                'HIV1-B-FR-K03455-seed'}
        assert regions == expected_regions, regions

    def check_2110(self):
        amino_rows = list(self.read_file('2110A-V3LOOP_S13', 'amino.csv'))
        assert amino_rows
        for row in amino_rows:
            pos = int(row['refseq.aa.pos'])
            if row['region'] == 'GP120':
                if pos == 304:
                    check_amino_row(row, coverage=23, ins=4)
                continue
            assert row['region'] == 'V3LOOP', row['region']
            if pos == 6:
                check_amino_row(row, coverage=23, stop_codons=1)
            elif pos == 7:
                check_amino_row(row, coverage=21, low_quality=2)
            elif pos == 8:
                check_amino_row(row, coverage=20, partial=3)
            elif pos == 9:
                check_amino_row(row, coverage=23)
            elif pos == 10:
                check_amino_row(row, coverage=23, dels=4)
            elif pos < 18:
                check_amino_row(row, coverage=23)
            else:
                check_amino_row(row)

    def check_2120(self):
        amino_rows = list(self.read_file('2120A-PR_S14', 'amino.csv'))
        assert amino_rows
        found_regions = set(map(itemgetter('region'), amino_rows))
        assert found_regions == {'PR'}, found_regions
        for row in amino_rows:
            pos = int(row['refseq.aa.pos'])
            if pos <= 29:
                check_amino_row(row)
            elif pos <= 46:
                check_amino_row(row, coverage=23)
            elif pos <= 49:
                check_amino_row(row, coverage=23, stop_codons=1)
            elif pos <= 52:
                check_amino_row(row, coverage=21, low_quality=2)
            elif pos <= 54:
                check_amino_row(row, coverage=20, partial=3)
            elif pos <= 55:
                check_amino_row(row, coverage=23)
            elif pos <= 56:
                check_amino_row(row, coverage=23, ins=4)
            elif pos <= 58:
                check_amino_row(row, coverage=23)
            elif pos <= 64:
                check_amino_row(row, coverage=20, clip=3)
            elif pos <= 67:
                check_amino_row(row, coverage=20, dels=1)
            elif pos <= 90:
                check_amino_row(row, coverage=20)
            else:
                check_amino_row(row)

    def check_2130(self):
        conseq_rows = list(self.read_file('2130A-HCV_S15', 'conseq.csv'))
        regions = set(map(itemgetter('region'), conseq_rows))
        expected_regions = ({'1-HCV-2a'}
                            if self.is_denovo
                            else {'HCV-2a'})
        assert regions == expected_regions, regions

    def check_2130midi(self):
        conseq_rows = list(self.read_file('2130AMIDI-MidHCV_S16', 'conseq.csv'))
        regions = set(map(itemgetter('region'), conseq_rows))
        expected_regions = ({'1-HCV-2a'}
                            if self.is_denovo
                            else {'HCV-2a'})
        assert regions == expected_regions, regions

    def check_2140(self):
        resistance_rows = list(self.read_file('2140A-HIV_S17', 'resistance.csv'))
        pr_rows = [row for row in resistance_rows if row['region'] == 'PR']
        assert pr_rows
        for row in pr_rows:
            assert row['level'] != '0', (row['drug_name'], row['level_name'])
            if row['drug'] == 'IDV/r':
                assert row['level'] == '3', row['level']

    def check_2160(self):
        amino_rows = list(self.read_file('2160A-HCV_S19', 'amino.csv'))
        assert amino_rows
        for row in amino_rows:
            assert row['region'] == 'HCV2-JFH-1-NS5b', row['region']
            pos = int(row['refseq.aa.pos'])
            coverage = int(row['coverage'])
            coverage_message = f'{coverage} coverage at {pos}'
            if pos < 30:
                assert coverage < 10, coverage_message
            elif 70 < pos < 150:
                assert 10 < coverage, coverage_message
            elif 230 < pos:
                assert coverage < 10, coverage_message

    def check_2160midi(self):
        amino_rows = list(self.read_file('2160AMIDI-MidHCV_S20', 'amino.csv'))
        assert amino_rows
        for row in amino_rows:
            assert row['region'] == 'HCV2-JFH-1-NS5b', row['region']
            pos = int(row['refseq.aa.pos'])
            coverage = int(row['coverage'])
            coverage_message = f'{coverage} coverage at {pos}'
            if pos < 370:
                assert coverage < 10, coverage_message
            elif 400 <= pos < 520:
                assert 10 <= coverage, coverage_message
            elif 540 <= pos:
                assert coverage < 10, coverage_message

    def check_2170(self):
        amino_rows = list(self.read_file('2170A-HCV_S21', 'amino.csv'))
        assert amino_rows
        for row in amino_rows:
            pos = int(row['refseq.aa.pos'])
            coverage = int(row['coverage'])
            coverage_message = f'{coverage} coverage at {pos}'
            if row['region'] == 'HCV1A-H77-NS5a':
                if pos < 15:
                    assert coverage < 10, coverage_message
                elif 50 <= pos:
                    assert 10 < coverage, coverage_message
            elif row['region'] == 'HCV1A-H77-NS5b':
                if pos <= 570:
                    assert 10 < coverage, coverage_message
            elif row['region'] == 'HCV2-JFH-1-NS5a':
                if pos < 15:
                    assert coverage < 10, coverage_message
                elif 100 <= pos:
                    assert 10 < coverage, coverage_message
            else:
                assert row['region'] == 'HCV2-JFH-1-NS5b', row['region']
                if pos < 540:
                    assert 10 < coverage, coverage_message

    def check_2180(self):
        amino_rows = list(self.read_file('2180A-HIV_S22', 'amino.csv'))
        assert amino_rows
        for row in amino_rows:
            pos = int(row['refseq.aa.pos'])
            coverage = int(row['coverage'])
            coverage_message = f'{coverage} coverage at {pos}'
            if row['region'] == 'GP120':
                assert row['seed'] == 'HIV1-B-FR-K03455-seed', row['seed']
                if pos < 50:
                    assert coverage < 10, coverage_message
                elif 230 < pos < 330:
                    assert 10 < coverage, coverage_message
                elif 380 < pos:
                    assert coverage < 10, coverage_message
            elif row['region'] == 'V3LOOP':
                assert row['seed'] == 'HIV1-CON-XX-Consensus-seed', row['seed']
                assert 10 < coverage, coverage_message
            else:
                # Last 2 codons of vpu have coverage.
                assert row['region'] == 'HIV1B-vpu', row['region']
                pos = row['query.nuc.pos']
                assert pos in ('6305', '6308'), pos

    def check_2190(self):
        mutation_rows = self.read_file('2190A-SARSCOV2_S23', 'nuc_mutations.csv')
        mutations = [''.join(fields)
                     for fields in map(itemgetter('wt', 'refseq_nuc_pos', 'var'),
                                       mutation_rows)]
        assert mutations == ['T13199C', 'T23C', 'T23C'], mutations

    def check_2200(self):
        amino_rows = list(self.read_file('2200A-SARSCOV2_S24', 'amino.csv'))
        assert amino_rows
        for row in amino_rows:
            pos = int(row['refseq.aa.pos'])
            coverage = int(row['coverage'])
            coverage_message = f'{coverage} coverage at {pos}'
            assert row['region'] in ('SARS-CoV-2-nsp1', 'SARS-CoV-2-ORF1ab'), row
            if 27 <= pos <= 102:
                assert coverage == 100, coverage_message
            else:
                assert coverage == 0, coverage_message

    def check_2210(self):
        conseq_rows = [row
                       for row in self.read_file('2210A-NFLHIVDNA_S25',
                                                 'conseq_all.csv')
                       if not row['region']]  # Whole-genome conseq only.
        assert len(conseq_rows) == 1, len(conseq_rows)
        conseq = conseq_rows[0]['sequence']
        assert 585 <= len(conseq), len(conseq)


def gzip_compress(source_path: Path, target_path: Path):
    with source_path.open('rb') as source:
        with gzip.open(target_path, 'wb') as target:
            shutil.copyfileobj(source, target)


def create_sample_scratch(fastq_file):
    sample_name = '_'.join(str(fastq_file.name).split('_')[:2])
    scratch_path: Path = fastq_file.parent / 'scratch' / sample_name
    scratch_path.mkdir(parents=True, exist_ok=True)
    return scratch_path


class SampleRunner:
    def __init__(self, image_path: Path, is_docker: bool = False):
        self.image_path = image_path
        self.is_docker = is_docker
        self.is_denovo = False
        self.bad_cycles_path = None

    def process_quality(self, quality_file: Path):
        scratch_path: Path = quality_file.parent / 'scratch' / 'quality'
        scratch_path.mkdir(parents=True)
        input_path = scratch_path / 'input'
        output_path = scratch_path / 'output'
        input_path.mkdir()
        output_path.mkdir()
        quality_input = input_path / quality_file.name
        self.bad_cycles_path = output_path / 'bad_cycles.csv'
        shutil.copy(str(quality_file), str(quality_input))
        run(self.build_command([quality_input],
                               [self.bad_cycles_path],
                               app_name='filter_quality'),
            stdout=PIPE,
            stderr=PIPE,
            check=True)

    def process_sample(self, fastq_file: Path):
        fastq_file2 = fastq_file.parent / (fastq_file.name.replace('_R1_',
                                                                   '_R2_'))
        scratch_path = create_sample_scratch(fastq_file)
        sample_name = scratch_path.name
        if self.is_denovo:
            input_path = scratch_path / 'input_denovo'
            output_path = scratch_path / 'output_denovo'
        else:
            input_path = scratch_path / 'input'
            output_path = scratch_path / 'output'
        input_path.mkdir()
        output_path.mkdir()
        fastq_input1 = input_path / fastq_file.name
        fastq_input2 = input_path / fastq_file2.name
        gzip_compress(fastq_file, fastq_input1)
        gzip_compress(fastq_file2, fastq_input2)
        sample_info_path = input_path / 'sample_info.csv'
        with sample_info_path.open('w') as f:
            writer = DictWriter(f, ['sample', 'project'])
            writer.writeheader()
            sections = sample_name.split('_')
            fields = sections[0].split('-')
            project_code = fields[-1]
            writer.writerow(dict(sample=sample_name, project=project_code))

        if self.is_denovo:
            output_names = [
                'g2p.csv',
                'g2p_summary.csv',
                'remap_counts.csv',
                'remap_conseq.csv',
                'unmapped1.fastq',
                'unmapped2.fastq',
                'conseq_ins.csv',
                'failed.csv',
                'cascade.csv',
                'nuc.csv',
                'amino.csv',
                'insertions.csv',
                'conseq.csv',
                'conseq_all.csv',
                'concordance.csv',
                'concordance_seed.csv',
                'failed_align.csv',
                'coverage_scores.csv',
                'coverage_maps.tar',
                'aligned.csv',
                'g2p_aligned.csv',
                'genome_coverage.csv',
                'genome_coverage.svg',
                'genome_concordance.svg',
                'unstitched_cascade.csv',
                'unstitched_conseq.csv',
                'unstitched_contigs.csv',
                'contigs.csv',
                'stitcher_plot_svg',
                'read_entropy.csv',
                'conseq_region.csv',
            ]

        else:
            output_names = [
                'g2p.csv',
                'g2p_summary.csv',
                'remap_counts.csv',
                'remap_conseq.csv',
                'unmapped1.fastq',
                'unmapped2.fastq',
                'conseq_ins.csv',
                'failed.csv',
                'cascade.csv',
                'nuc.csv',
                'amino.csv',
                'insertions.csv',
                'conseq.csv',
                'conseq_all.csv',
                'concordance.csv',
                'concordance_seed.csv',
                'failed_align.csv',
                'coverage_scores.csv',
                'coverage_maps.tar',
                'aligned.csv',
                'g2p_aligned.csv',
                'genome_coverage.csv',
                'genome_coverage.svg',
                'genome_concordance.svg',
            ]

        output_paths = [output_path/name for name in output_names]
        app_name = 'denovo' if self.is_denovo else None
        run_with_retries(self.build_command([sample_info_path,
                                             fastq_input1,
                                             fastq_input2,
                                             self.bad_cycles_path],
                                            output_paths,
                                            app_name))

        for path in output_paths:

            if path == (output_path/"conseq_ins.csv"):
                # This file is special. See https://github.com/cfe-lab/MiCall/issues/1085
                path = output_path/"scratch"/"conseq_ins.csv"

            assert os.path.exists(path), f"Expected output file {path!r} to be created."

        return sample_name

    def process_resistance(self, sample_group: SampleGroup):
        main_fastq_path = sample_group.names[0]
        midi_fastq_path = sample_group.names[1]
        main_scratch = create_sample_scratch(main_fastq_path)
        if self.is_denovo:
            output_path: Path = main_scratch / 'output_denovo'
            input_path: Path = main_scratch / 'input_resistance_denovo'
        else:
            output_path: Path = main_scratch / 'output'
            input_path: Path = main_scratch / 'input_resistance'
        main_amino_path = output_path / 'amino.csv'
        main_nuc_path = output_path / 'nuc.csv'
        input_path.mkdir()
        main_amino_input = input_path / 'main_amino.csv'
        main_nuc_input = input_path / 'main_nuc.csv'
        shutil.copy(str(main_amino_path), str(main_amino_input))
        shutil.copy(str(main_nuc_path), str(main_nuc_input))
        if midi_fastq_path is None:
            midi_amino_input = main_amino_input
        else:
            midi_scratch = create_sample_scratch(midi_fastq_path)
            if self.is_denovo:
                midi_amino_path = midi_scratch / 'output_denovo' / 'amino.csv'
            else:
                midi_amino_path = midi_scratch / 'output' / 'amino.csv'
            midi_amino_input = input_path / 'midi_amino.csv'
            shutil.copy(str(midi_amino_path), str(midi_amino_input))
        if self.is_denovo:
            output_path2 = output_path.parent / 'output_resistance_denovo'
        else:
            output_path2 = output_path.parent / 'output_resistance'
        output_path2.mkdir()
        output_names = ['resistance.csv',
                        'mutations.csv',
                        'nuc_mutations.csv',
                        'resistance_fail.csv',
                        'resistance.pdf',
                        'resistance_consensus.csv']
        output_paths = [output_path2/name for name in output_names]
        command_args = self.build_command([main_amino_input,
                                           midi_amino_input,
                                           main_nuc_input],
                                          output_paths,
                                          app_name='resistance')
        # print(*command_args)
        run_with_retries(command_args)
        return sample_group.enum

    def build_command(self, inputs, outputs, app_name=None):
        input_path = inputs[0].parent
        output_path = outputs[0].parent

        if self.is_docker:
            command = ['docker',
                       'run',
                       '--rm',
                       '--read-only',
                       '--volume', '{}:/mnt/input'.format(input_path),
                       '--volume', '{}:/mnt/output'.format(output_path),
                       '--volume', '{}:/tmp'.format(output_path / 'tmp'),
                       '--entrypoint', 'micall',
                       '--', str(self.image_path),
                       ]

            app_arguments = {
                None: ['micall_kive'],
                'filter_quality': ['filter_quality'],
                'resistance': ['micall_kive_resistance'],
                'denovo': ['micall_kive', '--denovo'],
            }[app_name]

            command.extend(app_arguments)

        else:
            command = ['singularity',
                       'run',
                       '--contain',
                       '--cleanenv',
                       '-B',
                       '{}:/mnt/input,{}:/mnt/output'.format(input_path,
                                                             output_path)]
            if app_name:
                command.append('--app')
                command.append(app_name)
            command.append(str(self.image_path))

        for arguments, guest_path in zip((inputs, outputs),
                                         ('/mnt/input', '/mnt/output')):
            for argument in arguments:
                command.append(os.path.join(guest_path, argument.name))

        return command


def find_full_groups(fastq_files, sandbox_path):
    groups = list(find_groups([p.name for p in fastq_files],
                              sandbox_path / 'SampleSheet.csv'))
    full_groups = []
    for group in groups:
        full_names = tuple(name and (sandbox_path / name)
                           for name in group.names)
        full_groups.append(SampleGroup(group.enum,
                                       full_names,
                                       group.project_codes))
    return full_groups


def run_with_retries(command_args: typing.List[str], retries=2):
    while retries:
        try:
            run(command_args, stdout=PIPE, stderr=PIPE, check=True)
            break
        except CalledProcessError as ex:
            if ex.returncode != 255:
                raise
            # Singularity 2.5 seems to be flaky, just retry.
            retries -= 1


def main():
    parser = ArgumentParser(description='Validate with small test samples.',
                            formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('--sandbox',
                        help='Folder to copy microtest samples into.')
    parser.add_argument('--sample',
                        help='Prefix of sample name to run.')
    parser.add_argument('--docker', action='store_true',
                        help='If to use docker instead of Singularity.')
    # noinspection PyTypeChecker
    parser.add_argument('image',
                        type=Path,
                        help='Singularity image to run tests in.')
    args = parser.parse_args()
    source_path: Path = Path(__file__).parent.parent / 'tests' / 'microtest'
    if args.sandbox is None:
        sandbox_path = source_path
        shutil.rmtree(source_path / 'scratch', ignore_errors=True)
        all_sandboxes = []
    else:
        sandbox_path = Path(mkdtemp(prefix='microtest_', dir=args.sandbox))
        all_sandboxes = Path(args.sandbox).glob('microtest_*')
        print(f'Copying to {sandbox_path}.')
        source_files = list(source_path.glob('*.fastq'))
        source_files.append(source_path / 'quality.csv')
        source_files.append(source_path / 'SampleSheet.csv')
        for source_file in source_files:
            target_file: Path = sandbox_path / source_file.name
            shutil.copy(str(source_file), str(target_file))
    runner = SampleRunner(args.image, is_docker=args.docker)
    with ProcessPoolExecutor() as pool:
        search_pattern = '*_R1_*.fastq'
        if args.sample:
            search_pattern = args.sample + search_pattern
        fastq_files = sorted(sandbox_path.glob(search_pattern))
        sample_groups = find_full_groups(fastq_files, sandbox_path)
        try:
            runner.process_quality(sandbox_path / 'quality.csv')
            for sample in pool.map(runner.process_sample, fastq_files):
                print(f'Processed {sample}.')
            for sample_enum in pool.map(runner.process_resistance, sample_groups):
                print(f'Processed resistance for {sample_enum}.')
            runner.is_denovo = True
            for sample in pool.map(runner.process_sample, fastq_files):
                print(f'Processed denovo {sample}.')
            for sample_enum in pool.map(runner.process_resistance, sample_groups):
                print(f'Processed resistance for denovo {sample_enum}.')
        except CalledProcessError as ex:
            print(ex.stderr.decode('utf8'))
            raise

    results_folder = ResultsFolder(sandbox_path, is_denovo=False)
    folder_methods = inspect.getmembers(results_folder, inspect.ismethod)
    check_pattern = 'check_'
    if args.sample:
        check_pattern += args.sample
    check_methods = [(name, method)
                     for name, method in folder_methods
                     if name.startswith(check_pattern)]
    error_count = 0
    for name, method in check_methods:
        if not inspect.signature(method).parameters:
            # Doesn't require any parameters, so it's a top level check method.
            print(name)
            results_folder.is_denovo = False
            try:
                method()
            except AssertionError:
                print_exc(file=sys.stdout)
                error_count += 1

            print(name, '(denovo)')
            results_folder.is_denovo = True
            try:
                method()
            except AssertionError:
                print_exc(file=sys.stdout)
                error_count += 1
    print('Finished checking.')
    assert error_count == 0, error_count

    # If they all passed, clean up.
    if all_sandboxes:
        print('Cleaning up sandboxes:')
        for sandbox in all_sandboxes:
            print(' ', sandbox)
            shutil.rmtree(sandbox)
    print('Done.')

if __name__ == '__main__':
    main()
