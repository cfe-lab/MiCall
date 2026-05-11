import argparse
import csv
import os
import shutil
import tarfile
import tempfile
from pathlib import Path
from typing import Sequence, TextIO

DOWNLOADED_RESULTS = [
    'remap_counts_csv',
    'conseq_csv',
    'conseq_all_csv',
    'conseq_region_csv',
    'concordance_csv',
    'concordance_seed_csv',
    'insertions_csv',
    'failed_csv',
    'nuc_csv',
    'amino_csv',
    'failed_align_csv',
    'g2p_csv',
    'g2p_summary_csv',
    'coverage_scores_csv',
    'coverage_maps_tar',
    'cascade_csv',
    'mixed_counts_csv',
    'mixed_amino_csv',
    'mixed_amino_merged_csv',
    'resistance_csv',
    'mutations_csv',
    'nuc_mutations_csv',
    'resistance_fail_csv',
    'resistance_consensus_csv',
    'wg_fasta',
    'mid_fasta',
    'unstitched_cascade_csv',
    'unstitched_conseq_csv',
    'unstitched_contigs_csv',
    'contigs_csv',
    'stitcher_plot_svg',
    'alignment_svg',
    'alignment_png',
    'assembly_fasta',
    'genome_coverage_csv',
    'genome_coverage_svg',
    'genome_concordance_svg',
    'read_entropy_csv',
    'outcome_summary_csv',
    'conseqs_primers_csv',
    'contigs_primers_csv',
    'table_precursor_csv',
    'proviral_landscape_csv',
    'hivseqinr_results_tar',
    'detailed_results_tar',
]

REQUIRED_MANIFEST_COLUMNS = {'index', 'sample', 'output_name'}


def get_output_filename(output_name: str) -> str:
    return '.'.join(output_name.rsplit('_', 1))


def extract_csv(source: TextIO, target: TextIO, sample_name: str, source_count: int) -> int:
    reader = csv.DictReader(source)
    fieldnames = reader.fieldnames
    if fieldnames is None:
        return 0
    fieldnames = list(fieldnames)
    has_sample = 'sample' in fieldnames
    if not has_sample:
        fieldnames.insert(0, 'sample')
    writer = csv.DictWriter(target, fieldnames, lineterminator='\n')
    if source_count == 0:
        writer.writeheader()
    for row in reader:
        if not has_sample:
            row['sample'] = sample_name
        writer.writerow(row)
    return 1


def extract_fasta(source: TextIO, target: TextIO, sample_name: str) -> int:
    for line in source:
        if line.startswith('>'):
            target.write(f'>{sample_name},{line[1:]}')
        else:
            target.write(line)
    return 1


def remove_empty_directory(path: Path) -> None:
    if path.exists() and not any(path.iterdir()):
        path.rmdir()


def extract_coverage_maps(sample_names: list[str], scratch_path: Path, results_path: Path) -> None:
    coverage_path = results_path / 'coverage_maps'
    coverage_path.mkdir(parents=True, exist_ok=True)
    for sample_name in sample_names:
        source_path = scratch_path / sample_name / 'coverage_maps.tar'
        try:
            with tarfile.open(source_path) as source_tar:
                for source_info in source_tar:
                    filename = os.path.basename(source_info.name)
                    target_path = coverage_path / f'{sample_name}.{filename}'
                    source = source_tar.extractfile(source_info)
                    if source is None:
                        continue
                    with source, open(target_path, 'wb') as target:
                        shutil.copyfileobj(source, target)
        except FileNotFoundError:
            continue
    remove_empty_directory(coverage_path)


def extract_archive(sample_names: list[str], scratch_path: Path, results_path: Path, output_name: str) -> None:
    archive_name = output_name[:-4]
    output_path = results_path / archive_name
    output_path.mkdir(parents=True, exist_ok=True)
    for sample_name in sample_names:
        source_path = scratch_path / sample_name / f'{archive_name}.tar'
        try:
            with tarfile.open(source_path) as source_tar:
                sample_target_path = output_path / sample_name
                sample_target_path.mkdir(parents=True, exist_ok=True)
                for source_info in source_tar:
                    filename = os.path.basename(source_info.name)
                    target_path = sample_target_path / filename
                    source = source_tar.extractfile(source_info)
                    if source is None:
                        continue
                    with source, open(target_path, 'wb') as target:
                        shutil.copyfileobj(source, target)
                remove_empty_directory(sample_target_path)
        except FileNotFoundError:
            continue
    remove_empty_directory(output_path)


def move_alignment_plot(sample_names: list[str], extension: str, scratch_path: Path, results_path: Path) -> None:
    alignment_path = results_path / 'alignment'
    alignment_path.mkdir(parents=True, exist_ok=True)
    for sample_name in sample_names:
        source_path = scratch_path / sample_name / f'alignment{extension}'
        target_path = alignment_path / f'{sample_name}_alignment{extension}'
        if source_path.exists():
            source_path.rename(target_path)
    remove_empty_directory(alignment_path)


def move_genome_coverage(sample_names: list[str], scratch_path: Path, results_path: Path) -> None:
    plots_path = results_path / 'genome_coverage'
    plots_path.mkdir(parents=True, exist_ok=True)
    for sample_name in sample_names:
        source_path = scratch_path / sample_name / 'genome_coverage.svg'
        target_path = plots_path / f'{sample_name}_genome_coverage.svg'
        if source_path.exists():
            source_path.rename(target_path)
        concordance_path = scratch_path / sample_name / 'genome_concordance.svg'
        target_concordance_path = plots_path / f'{sample_name}_genome_concordance.svg'
        if concordance_path.exists():
            concordance_path.rename(target_concordance_path)
    remove_empty_directory(plots_path)


def move_stitcher_plot(sample_names: list[str], scratch_path: Path, results_path: Path) -> None:
    plots_path = results_path / 'stitcher_plots'
    plots_path.mkdir(parents=True, exist_ok=True)
    for sample_name in sample_names:
        source_path = scratch_path / sample_name / 'stitcher_plot.svg'
        target_path = plots_path / f'{sample_name}_stitcher_plot.svg'
        if source_path.exists():
            source_path.rename(target_path)
    remove_empty_directory(plots_path)


def copy_outputs(sample_names: list[str], scratch_path: Path, results_path: Path) -> None:
    results_path.mkdir(parents=True, exist_ok=True)
    for output_name in DOWNLOADED_RESULTS:
        if output_name == 'coverage_maps_tar':
            extract_coverage_maps(sample_names, scratch_path, results_path)
            continue
        if output_name.endswith('_tar'):
            extract_archive(sample_names, scratch_path, results_path, output_name)
            continue
        if output_name == 'alignment_svg':
            move_alignment_plot(sample_names, '.svg', scratch_path, results_path)
            continue
        if output_name == 'alignment_png':
            move_alignment_plot(sample_names, '.png', scratch_path, results_path)
            continue
        if output_name == 'genome_coverage_svg':
            move_genome_coverage(sample_names, scratch_path, results_path)
            continue
        if output_name == 'stitcher_plot_svg':
            move_stitcher_plot(sample_names, scratch_path, results_path)
            continue

        source_count = 0
        filename = get_output_filename(output_name)
        target_path = results_path / filename
        with target_path.open('w', encoding='utf-8', newline='') as target:
            for sample_name in sample_names:
                source_path = scratch_path / sample_name / filename
                if not source_path.exists():
                    continue
                with source_path.open(encoding='utf-8', newline='') as source:
                    if output_name.endswith('_fasta'):
                        source_count += extract_fasta(source, target, sample_name)
                    else:
                        source_count += extract_csv(source, target, sample_name, source_count)
        if source_count == 0 and target_path.exists():
            target_path.unlink()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument('--run_outputs', nargs='*', default=[], type=Path)
    parser.add_argument('metadata_csv', type=Path)
    parser.add_argument('collated_results_tar', type=Path)
    return parser.parse_args()


def stage_inputs_by_sample(run_outputs: Sequence[Path], metadata_csv: Path, scratch_path: Path) -> list[str]:
    with metadata_csv.open(newline='', encoding='utf-8') as manifest_file:
        reader = csv.DictReader(manifest_file)
        if reader.fieldnames is None:
            raise ValueError('Metadata manifest is missing a header row.')
        missing_columns = REQUIRED_MANIFEST_COLUMNS - set(reader.fieldnames)
        if missing_columns:
            missing = ', '.join(sorted(missing_columns))
            raise ValueError(f'Metadata manifest is missing required columns: {missing}.')
        rows = list(reader)

    sample_names: set[str] = set()
    for row_number, row in enumerate(rows, start=2):
        index_text = row.get('index', '')
        try:
            file_index = int(index_text)
        except (TypeError, ValueError) as ex:
            raise ValueError(
                f'Metadata manifest row {row_number} has invalid index {index_text!r}.') from ex

        sample_name = (row.get('sample') or '').strip()
        output_name = (row.get('output_name') or '').strip()
        if file_index < 0 or file_index >= len(run_outputs):
            raise ValueError(
                f'Invalid run_outputs index {file_index} in metadata manifest row {row_number}.')
        if not sample_name:
            raise ValueError(f'Metadata manifest row {row_number} has empty sample value.')
        if (Path(sample_name).name != sample_name or
            sample_name in ('.', '..') or
            '/' in sample_name or
            '\\' in sample_name):
            raise ValueError(
                f'Metadata manifest row {row_number} has invalid sample name {sample_name!r}.')
        if output_name not in DOWNLOADED_RESULTS:
            raise ValueError(
                f'Metadata manifest row {row_number} has unknown output_name {output_name!r}.')

        source_path = run_outputs[file_index]
        if not source_path.is_file():
            raise FileNotFoundError(
                f'Metadata manifest row {row_number} references missing run output {source_path}.')
        sample_path = scratch_path / sample_name
        sample_path.mkdir(parents=True, exist_ok=True)
        target_path = sample_path / get_output_filename(output_name)
        if target_path.exists():
            raise ValueError(
                f'Metadata manifest row {row_number} duplicates output {output_name!r} '
                f'for sample {sample_name!r}.')
        shutil.copyfile(source_path, target_path)
        sample_names.add(sample_name)

    return sorted(sample_names)


def main() -> None:
    args = parse_args()
    args.collated_results_tar.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory() as tmp_text:
        tmp_path = Path(tmp_text)
        scratch_path = tmp_path / 'scratch'
        collated_path = tmp_path / 'collated'
        scratch_path.mkdir(parents=True, exist_ok=True)
        collated_path.mkdir(parents=True, exist_ok=True)

        sample_names = stage_inputs_by_sample(args.run_outputs,
                              args.metadata_csv,
                              scratch_path)
        copy_outputs(sample_names, scratch_path, collated_path)

        with tarfile.open(args.collated_results_tar, 'w') as output_tar:
            for file_path in sorted(collated_path.rglob('*')):
                if file_path.is_file():
                    output_tar.add(file_path, file_path.relative_to(collated_path))


if __name__ == '__main__':
    main()
