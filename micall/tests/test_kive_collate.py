import sys
import tarfile
from pathlib import Path

import pytest

from micall.utils import kive_collate


def test_parse_args_with_optional_multiple_and_separator(monkeypatch, tmp_path):
    metadata_path = tmp_path / 'metadata.csv'
    metadata_path.write_text('index,sample,output_name\n')
    output_path = tmp_path / 'out.tar'
    monkeypatch.setattr(
        sys,
        'argv',
        ['kive_collate', '--inputs', 'a.csv', 'b.csv', str(metadata_path), '--', str(output_path)])

    args = kive_collate.parse_args()

    assert args.inputs == [Path('a.csv'), Path('b.csv'), metadata_path]
    assert args.output == output_path
    assert not args.verbose
    assert not args.debug
    assert not args.quiet


def test_parse_args_with_debug_flag(monkeypatch, tmp_path):
    metadata_path = tmp_path / 'metadata.csv'
    metadata_path.write_text('index,sample,output_name\n')
    output_path = tmp_path / 'out.tar'

    monkeypatch.setattr(
        sys,
        'argv',
        ['kive_collate', '--debug', '--inputs', 'a.csv', str(metadata_path), '--', str(output_path)])

    args = kive_collate.parse_args()

    assert args.debug
    assert not args.verbose
    assert not args.quiet


def test_main_collates_csv_and_fasta_from_multiple_samples(monkeypatch, tmp_path):
    sample1_cascade = tmp_path / 'sample1_cascade.csv'
    sample1_fasta = tmp_path / 'sample1_wg.fasta'
    sample2_cascade = tmp_path / 'sample2_cascade.csv'
    sample2_fasta = tmp_path / 'sample2_wg.fasta'
    sample1_cascade.write_text('x,y\n1,2\n')
    sample2_cascade.write_text('x,y\n3,4\n')
    sample1_fasta.write_text('>seed\nACTG\n')
    sample2_fasta.write_text('>seed\nACTG\n')
    metadata_path = tmp_path / 'metadata.csv'
    metadata_path.write_text(
        'index,sample,output_name\n'
        '0,E11111,cascade_csv\n'
        '1,E11111,wg_fasta\n'
        '2,E22222,cascade_csv\n'
        '3,E22222,wg_fasta\n')
    output_path = tmp_path / 'collated.tar'

    monkeypatch.setattr(
        sys,
        'argv',
        ['kive_collate', '--inputs',
         str(sample1_cascade),
         str(sample1_fasta),
         str(sample2_cascade),
         str(sample2_fasta),
         str(metadata_path),
         '--',
         str(output_path)])

    kive_collate.main()

    extract_path = tmp_path / 'extract'
    extract_path.mkdir()
    with tarfile.open(output_path) as output_tar:
        output_tar.extractall(extract_path, filter='data')

    cascade_text = (extract_path / 'cascade.csv').read_text()
    assert cascade_text == (
        'sample,x,y\n'
        'E11111,1,2\n'
        'E22222,3,4\n'
    )

    fasta_text = (extract_path / 'wg.fasta').read_text()
    assert fasta_text == (
        '>E11111,seed\n'
        'ACTG\n'
        '>E22222,seed\n'
        'ACTG\n'
    )


def test_stage_inputs_by_sample_rejects_invalid_index(tmp_path):
    metadata_path = tmp_path / 'metadata.csv'
    metadata_path.write_text('index,sample,output_name\nabc,E11111,cascade_csv\n')
    run_outputs = [tmp_path / 'cascade.csv']
    run_outputs[0].write_text('x,y\n1,2\n')

    with pytest.raises(ValueError, match='invalid index'):
        kive_collate.stage_inputs_by_sample(run_outputs, metadata_path, tmp_path / 'scratch')


def test_stage_inputs_by_sample_rejects_invalid_sample_name(tmp_path):
    metadata_path = tmp_path / 'metadata.csv'
    metadata_path.write_text('index,sample,output_name\n0,../escape,cascade_csv\n')
    run_outputs = [tmp_path / 'cascade.csv']
    run_outputs[0].write_text('x,y\n1,2\n')

    with pytest.raises(ValueError, match='invalid sample name'):
        kive_collate.stage_inputs_by_sample(run_outputs, metadata_path, tmp_path / 'scratch')


def test_stage_inputs_by_sample_rejects_missing_required_columns(tmp_path):
    metadata_path = tmp_path / 'metadata.csv'
    metadata_path.write_text('index,sample\n0,E11111\n')
    run_outputs = [tmp_path / 'cascade.csv']
    run_outputs[0].write_text('x,y\n1,2\n')

    with pytest.raises(ValueError, match='missing required columns'):
        kive_collate.stage_inputs_by_sample(run_outputs, metadata_path, tmp_path / 'scratch')


def test_stage_inputs_by_sample_rejects_duplicate_output_for_sample(tmp_path):
    metadata_path = tmp_path / 'metadata.csv'
    metadata_path.write_text(
        'index,sample,output_name\n'
        '0,E11111,cascade_csv\n'
        '1,E11111,cascade_csv\n'
    )
    run_outputs = [tmp_path / 'cascade1.csv', tmp_path / 'cascade2.csv']
    run_outputs[0].write_text('x,y\n1,2\n')
    run_outputs[1].write_text('x,y\n3,4\n')

    with pytest.raises(ValueError, match='duplicates output'):
        kive_collate.stage_inputs_by_sample(run_outputs, metadata_path, tmp_path / 'scratch')


def test_parse_args_with_explicit_output_separator(monkeypatch, tmp_path):
    metadata_path = tmp_path / 'metadata.csv'
    metadata_path.write_text('index,sample,output_name\n')
    output_path = tmp_path / 'out.tar'
    monkeypatch.setattr(
        sys,
        'argv',
        ['kive_collate', '--inputs', 'a.csv', 'b.csv', str(metadata_path), '--', str(output_path)])

    args = kive_collate.parse_args()

    assert args.inputs == [Path('a.csv'), Path('b.csv'), metadata_path]
    assert args.output == output_path
