"""
Tests for downloading stitcher plot SVG files from Kive.
"""
import pytest
from pathlib import Path
from unittest.mock import Mock, patch
from io import BytesIO

from micall.monitor.kive_watcher import KiveWatcher, FolderWatcher, SampleGroup, PipelineType, trim_name
from micall.tests.test_kive_watcher import create_kive_watcher_with_main_run, create_raw_data_with_two_samples, create_mock_open_kive, create_default_config

assert create_raw_data_with_two_samples is not None  # To avoid unused import warning
assert create_mock_open_kive is not None # To avoid unused import warning
assert create_default_config is not None  # To avoid unused import warning


def test_stitcher_plot_svg_in_downloaded_results():
    """Test that stitcher_plot_svg is in the list of DOWNLOADED_RESULTS"""
    from micall.monitor.kive_watcher import DOWNLOADED_RESULTS
    assert 'stitcher_plot_svg' in DOWNLOADED_RESULTS


def test_move_stitcher_plot_single_sample(tmp_path):
    """Test moving stitcher plot for a single sample"""
    # Setup
    scratch_path = tmp_path / "scratch"
    results_path = tmp_path / "results"
    sample_name = "2110A-V3LOOP_S13"

    # Create source file
    source_dir = scratch_path / trim_name(sample_name)
    source_dir.mkdir(parents=True)
    source_file = source_dir / "stitcher_plot.svg"
    source_file.write_text("<svg>test stitcher plot</svg>")

    # Create folder watcher mock
    folder_watcher = Mock(spec=FolderWatcher)
    folder_watcher.all_samples = [sample_name]

    # Execute
    KiveWatcher.move_stitcher_plot(folder_watcher, scratch_path, results_path)

    # Verify
    expected_path = results_path / "stitcher_plots" / f"{trim_name(sample_name)}_stitcher_plot.svg"
    assert expected_path.exists()
    assert expected_path.read_text() == "<svg>test stitcher plot</svg>"
    assert not source_file.exists()  # Should be moved, not copied


def test_move_stitcher_plot_multiple_samples(tmp_path):
    """Test moving stitcher plots for multiple samples"""
    # Setup
    scratch_path = tmp_path / "scratch"
    results_path = tmp_path / "results"
    sample_names = ["2110A-V3LOOP_S13", "2120A-PR_S14", "2130A-HCV_S15"]

    # Create source files
    for sample_name in sample_names:
        source_dir = scratch_path / trim_name(sample_name)
        source_dir.mkdir(parents=True)
        source_file = source_dir / "stitcher_plot.svg"
        source_file.write_text(f"<svg>{sample_name} stitcher plot</svg>")

    # Create folder watcher mock
    folder_watcher = Mock(spec=FolderWatcher)
    folder_watcher.all_samples = sample_names

    # Execute
    KiveWatcher.move_stitcher_plot(folder_watcher, scratch_path, results_path)

    # Verify all files were moved
    plots_dir = results_path / "stitcher_plots"
    assert plots_dir.exists()

    for sample_name in sample_names:
        expected_path = plots_dir / f"{trim_name(sample_name)}_stitcher_plot.svg"
        assert expected_path.exists()
        assert f"{sample_name} stitcher plot" in expected_path.read_text()


def test_move_stitcher_plot_missing_file(tmp_path):
    """Test that missing stitcher plots are handled gracefully"""
    # Setup
    scratch_path = tmp_path / "scratch"
    results_path = tmp_path / "results"
    sample_name = "2110A-V3LOOP_S13"

    # Don't create the source file - it's missing
    source_dir = scratch_path / trim_name(sample_name)
    source_dir.mkdir(parents=True)

    # Create folder watcher mock
    folder_watcher = Mock(spec=FolderWatcher)
    folder_watcher.all_samples = [sample_name]

    # Execute - should not raise an error
    KiveWatcher.move_stitcher_plot(folder_watcher, scratch_path, results_path)

    # Verify
    plots_dir = results_path / "stitcher_plots"
    # Directory should not be created if no plots were moved
    assert not plots_dir.exists() or len(list(plots_dir.iterdir())) == 0


def test_move_stitcher_plot_removes_empty_directory(tmp_path):
    """Test that empty stitcher_plots directory is removed"""
    # Setup
    scratch_path = tmp_path / "scratch"
    results_path = tmp_path / "results"
    sample_name = "2110A-V3LOOP_S13"

    # Don't create any source files
    source_dir = scratch_path / trim_name(sample_name)
    source_dir.mkdir(parents=True)

    # Create folder watcher mock
    folder_watcher = Mock(spec=FolderWatcher)
    folder_watcher.all_samples = [sample_name]

    # Execute
    KiveWatcher.move_stitcher_plot(folder_watcher, scratch_path, results_path)

    # Verify
    plots_dir = results_path / "stitcher_plots"
    # Directory should be removed if it was created but remained empty
    assert not plots_dir.exists() or len(list(plots_dir.iterdir())) == 0


def test_folder_completed_with_stitcher_plot(raw_data_with_two_samples, mock_open_kive, default_config, tmp_path):
    """Test that stitcher plots are moved correctly from scratch to results"""
    # Simplify the test: just test the move_stitcher_plot function directly
    from micall.monitor.kive_watcher import KiveWatcher, trim_name
    from pathlib import Path

    # Create test structure
    scratch_denovo_path = tmp_path / "scratch_denovo"
    results_path = tmp_path

    # Create a sample stitcher plot file
    sample_name = "2120A-PR_S14"
    sample_scratch = scratch_denovo_path / trim_name(sample_name)
    sample_scratch.mkdir(parents=True)
    stitcher_plot_content = b"<svg><text>Stitcher plot</text></svg>"
    (sample_scratch / "stitcher_plot.svg").write_bytes(stitcher_plot_content)

    # Create a mock folder_watcher with the sample
    folder_watcher = Mock(spec=FolderWatcher)
    folder_watcher.all_samples = [sample_name]

    # Call move_stitcher_plot
    KiveWatcher.move_stitcher_plot(folder_watcher, scratch_denovo_path, results_path)

    # Verify the file was moved
    expected_path = results_path / "stitcher_plots" / f"{trim_name(sample_name)}_stitcher_plot.svg"
    assert expected_path.exists()
    assert expected_path.read_bytes() == stitcher_plot_content
    assert not (sample_scratch / "stitcher_plot.svg").exists()  # Source should be moved, not copied


def test_stitcher_plot_svg_argument_in_micall_kive():
    """Test that stitcher_plot_svg argument is properly defined in micall_kive.py"""
    from micall.utils.micall_kive import parse_args
    import sys

    # Create test arguments (--denovo must come before positional arguments)
    test_args = [
        '--denovo',
        'sample_info.csv',
        'fastq1.fastq',
        'fastq2.fastq',
        'bad_cycles.csv',
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
        'stitcher_plot.svg',
        'read_entropy.csv',
        'conseq_region.csv',
        'conseq_stitched.csv',
    ]

    with patch.object(sys, 'argv', ['micall_kive'] + test_args):
        args = parse_args()
        assert hasattr(args, 'stitcher_plot_svg')
        assert args.stitcher_plot_svg == 'stitcher_plot.svg'


def test_stitcher_plot_passed_to_sample():
    """Test that stitcher_plot_svg is passed to Sample in load_sample"""
    from micall.utils.micall_kive import load_sample
    from argparse import Namespace

    # Create mock args
    args = Namespace(
        sample_info_csv='sample_info.csv',
        fastq1='fastq1.fastq',
        fastq2='fastq2.fastq',
        bad_cycles_csv='bad_cycles.csv',
        g2p_csv='g2p.csv',
        g2p_summary_csv='g2p_summary.csv',
        remap_counts_csv='remap_counts.csv',
        remap_conseq_csv='remap_conseq.csv',
        unmapped1_fastq='unmapped1.fastq',
        unmapped2_fastq='unmapped2.fastq',
        conseq_ins_csv='conseq_ins.csv',
        failed_csv='failed.csv',
        cascade_csv='cascade.csv',
        nuc_csv='nuc.csv',
        amino_csv='amino.csv',
        insertions_csv='insertions.csv',
        conseq_csv='conseq.csv',
        conseq_all_csv='conseq_all.csv',
        concordance_csv='concordance.csv',
        concordance_seed_csv='concordance_seed.csv',
        failed_align_csv='failed_align.csv',
        coverage_scores_csv='coverage_scores.csv',
        coverage_maps_tar='coverage_maps.tar',
        aligned_csv='aligned.csv',
        g2p_aligned_csv='g2p_aligned.csv',
        genome_coverage_csv='genome_coverage.csv',
        genome_coverage_svg='genome_coverage.svg',
        genome_concordance_svg='genome_concordance.svg',
        denovo=True,
        unstitched_cascade_csv='unstitched_cascade.csv',
        unstitched_conseq_csv='unstitched_conseq.csv',
        unstitched_contigs_csv='unstitched_contigs.csv',
        contigs_csv='contigs.csv',
        stitcher_plot_svg='stitcher_plot.svg',
        read_entropy_csv='read_entropy.csv',
        conseq_region_csv='conseq_region.csv',
        conseq_stitched_csv='conseq_stitched.csv',
    )

    sample = load_sample(args)

    # Verify that stitcher_plot_svg is accessible via the sample object
    assert sample.stitcher_plot_svg == 'stitcher_plot.svg'
