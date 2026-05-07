import sys
import tarfile
from pathlib import Path

from micall.utils import kive_collate


def _make_sample_tar(root: Path, sample_name: str, cascade_text: str) -> Path:
    sample_dir = root / sample_name
    sample_dir.mkdir(parents=True)
    (sample_dir / 'cascade.csv').write_text(cascade_text)
    (sample_dir / 'wg.fasta').write_text('>seed\nACTG\n')

    tar_path = root / f'{sample_name}.tar'
    with tarfile.open(tar_path, 'w') as sample_tar:
        sample_tar.add(sample_dir, arcname=sample_name)
    return tar_path


def test_parse_args_with_optional_multiple_and_separator(monkeypatch, tmp_path):
    output_path = tmp_path / 'out.tar'
    monkeypatch.setattr(
        sys,
        'argv',
        ['kive_collate', '--sample_results_tars', 'a.tar', 'b.tar', '--', str(output_path)])

    args = kive_collate.parse_args()

    assert args.sample_results_tars == ['a.tar', 'b.tar']
    assert args.collated_results_tar == str(output_path)


def test_main_collates_csv_and_fasta_from_multiple_samples(monkeypatch, tmp_path):
    sample1_tar = _make_sample_tar(tmp_path, 'E11111', 'x,y\n1,2\n')
    sample2_tar = _make_sample_tar(tmp_path, 'E22222', 'x,y\n3,4\n')
    output_path = tmp_path / 'collated.tar'

    monkeypatch.setattr(
        sys,
        'argv',
        ['kive_collate', '--sample_results_tars', str(sample1_tar), str(sample2_tar), str(output_path)])

    kive_collate.main()

    extract_path = tmp_path / 'extract'
    extract_path.mkdir()
    with tarfile.open(output_path) as output_tar:
        output_tar.extractall(extract_path)

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
