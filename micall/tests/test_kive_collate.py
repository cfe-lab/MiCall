import sys
import tarfile

from micall.utils import kive_collate


def test_parse_args_with_optional_multiple_and_separator(monkeypatch, tmp_path):
    metadata_path = tmp_path / 'metadata.csv'
    metadata_path.write_text('index,sample,output_name\n')
    output_path = tmp_path / 'out.tar'
    monkeypatch.setattr(
        sys,
        'argv',
        ['kive_collate', '--run_outputs', 'a.csv', 'b.csv', '--', str(metadata_path), str(output_path)])

    args = kive_collate.parse_args()

    assert args.run_outputs == ['a.csv', 'b.csv']
    assert args.metadata_csv == str(metadata_path)
    assert args.collated_results_tar == str(output_path)


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
        ['kive_collate', '--run_outputs',
         str(sample1_cascade),
         str(sample1_fasta),
         str(sample2_cascade),
         str(sample2_fasta),
         '--',
         str(metadata_path),
         str(output_path)])

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
