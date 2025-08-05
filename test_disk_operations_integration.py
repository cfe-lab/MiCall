#!/usr/bin/env python3
"""
Additional integration tests for disk_operations to verify correctness
and ensure no breaking changes to the kive_watcher integration.
"""

import tempfile
from pathlib import Path
import errno
from unittest.mock import patch

import pytest

from micall.monitor import disk_operations


class TestKiveWatcherIntegration:
    """Test scenarios that mirror how kive_watcher uses disk operations."""

    def test_extract_coverage_maps_scenario(self):
        """Test the exact pattern used in extract_coverage_maps."""
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            results_path = temp_path / "results"
            coverage_path = results_path / "coverage_maps"
            
            # This mirrors the exact pattern from kive_watcher.extract_coverage_maps
            disk_operations.mkdir_p(coverage_path, exist_ok=True)
            assert coverage_path.exists()
            
            # Simulate creating some files (would normally come from tar extraction)
            test_file = coverage_path / "sample1.coverage.svg"
            disk_operations.write_text(test_file, "test content")
            
            # The remove_empty_directory should NOT remove the directory because it has files
            disk_operations.remove_empty_directory(coverage_path)
            assert coverage_path.exists()  # Should still exist because it has files
            assert test_file.exists()
            
            # Now remove the file and try again
            disk_operations.unlink(test_file)
            disk_operations.remove_empty_directory(coverage_path)
            assert not coverage_path.exists()  # Should be removed now that it's empty

    def test_extract_archive_scenario(self):
        """Test the pattern used in extract_archive method."""
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            results_path = temp_path / "results"
            output_path = results_path / "detailed_results"
            sample_target_path = output_path / "sample1"
            
            # Create the directory structure
            disk_operations.mkdir_p(sample_target_path, exist_ok=True)
            assert sample_target_path.exists()
            
            # Test empty directory removal
            disk_operations.remove_empty_directory(sample_target_path)
            assert not sample_target_path.exists()
            
            # Test with files - create again and add content
            disk_operations.mkdir_p(sample_target_path, exist_ok=True)
            test_file = sample_target_path / "results.txt"
            disk_operations.write_text(test_file, "important data")
            
            # Should not remove when there are files
            disk_operations.remove_empty_directory(sample_target_path)
            assert sample_target_path.exists()
            assert test_file.exists()
            
            # Remove files and parent directory
            disk_operations.unlink(test_file)
            disk_operations.remove_empty_directory(sample_target_path)
            assert not sample_target_path.exists()
            
            # Test parent directory cleanup
            disk_operations.remove_empty_directory(output_path)
            assert not output_path.exists()

    def test_alignment_plot_scenario(self):
        """Test the pattern used in move_alignment_plot method."""
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            results_path = temp_path / "results"
            alignment_path = results_path / "alignment"
            
            # Create alignment directory
            disk_operations.mkdir_p(alignment_path, exist_ok=True)
            
            # Test case where no alignment files are moved (empty directory)
            disk_operations.remove_empty_directory(alignment_path)
            assert not alignment_path.exists()
            
            # Test case where some files exist
            disk_operations.mkdir_p(alignment_path, exist_ok=True)
            alignment_file = alignment_path / "sample1_alignment.svg"
            disk_operations.write_text(alignment_file, "<svg>test</svg>")
            
            # Should not remove directory with files
            disk_operations.remove_empty_directory(alignment_path)
            assert alignment_path.exists()

    def test_network_drive_retry_simulation(self):
        """Test retry behavior with network drive issues."""
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            test_dir = temp_path / "test_dir"
            
            # Create directory
            disk_operations.mkdir_p(test_dir, exist_ok=True)
            
            # Mock intermittent permission errors at the Path.rmdir level
            call_count = 0
            original_rmdir = Path.rmdir
            
            def mock_rmdir(self):
                nonlocal call_count
                if str(self) == str(test_dir):
                    call_count += 1
                    if call_count <= 2:  # Fail first 2 attempts
                        raise OSError(errno.EACCES, "Permission denied")
                return original_rmdir(self)
            
            with patch.object(Path, 'rmdir', mock_rmdir):
                # Should succeed after retries
                disk_operations.remove_empty_directory(test_dir)
                assert not test_dir.exists()
                assert call_count == 3  # Should have retried and succeeded on 3rd attempt

    def test_remove_empty_directory_non_empty_case(self):
        """Test that ENOTEMPTY is handled silently as expected."""
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            parent_dir = temp_path / "parent"
            child_dir = parent_dir / "child"
            
            # Create nested structure
            disk_operations.mkdir_p(child_dir, exist_ok=True)
            
            # Add file to child
            test_file = child_dir / "important.txt"
            disk_operations.write_text(test_file, "data")
            
            # Try to remove parent - should silently handle ENOTEMPTY
            disk_operations.remove_empty_directory(parent_dir)
            assert parent_dir.exists()  # Should still exist
            assert child_dir.exists()   # Child should still exist
            assert test_file.exists()   # File should still exist

    def test_disk_file_operation_integration(self):
        """Test the disk_file_operation context manager pattern used in kive_watcher."""
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            target_file = temp_path / "test_output.csv"
            
            # Test write operation (used in copy_outputs)
            with disk_operations.disk_file_operation(target_file, 'w') as target:
                target.write("sample,value\n")
                target.write("sample1,123\n")
            
            assert target_file.exists()
            
            # Test read operation (used in copy_outputs)
            with disk_operations.disk_file_operation(target_file, 'r') as source:
                content = source.read()
                assert "sample,value" in content
                assert "sample1,123" in content

    def test_error_handling_write_text(self):
        """Test error handling for write_text operation."""
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            error_file = temp_path / "errorprocessing"
            
            # Test normal operation
            disk_operations.write_text(error_file, "Finding error metrics failed.\n")
            assert error_file.exists()
            content = error_file.read_text()
            assert content == "Finding error metrics failed.\n"

    def test_concurrent_directory_operations(self):
        """Test handling of concurrent directory operations."""
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            
            # Create multiple directories as might happen in concurrent processing
            dirs = []
            for i in range(5):
                dir_path = temp_path / f"sample_{i}"
                disk_operations.mkdir_p(dir_path, exist_ok=True)
                dirs.append(dir_path)
            
            # Remove them all
            for dir_path in dirs:
                disk_operations.remove_empty_directory(dir_path)
                assert not dir_path.exists()

    def test_rmtree_with_ignore_errors(self):
        """Test rmtree operation with ignore_errors=True as used in kive_watcher."""
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            results_path = temp_path / "Results" / "version_1.0"
            
            # Create complex directory structure
            disk_operations.mkdir_p(results_path / "sub1" / "sub2", exist_ok=True)
            disk_operations.write_text(results_path / "file1.txt", "content")
            disk_operations.write_text(results_path / "sub1" / "file2.txt", "content")
            
            # Test rmtree with ignore_errors=True
            disk_operations.rmtree(results_path, ignore_errors=True)
            assert not results_path.exists()


def test_stress_test_remove_empty_directory():
    """Stress test remove_empty_directory with many operations."""
    with tempfile.TemporaryDirectory() as temp_dir:
        temp_path = Path(temp_dir)
        
        # Create and remove many directories
        for i in range(100):
            test_dir = temp_path / f"dir_{i}"
            disk_operations.mkdir_p(test_dir, exist_ok=True)
            
            if i % 2 == 0:
                # Even numbered dirs - leave empty
                disk_operations.remove_empty_directory(test_dir)
                assert not test_dir.exists()
            else:
                # Odd numbered dirs - add file then remove
                test_file = test_dir / "temp.txt"
                disk_operations.write_text(test_file, "temp")
                
                # Should not remove (has file)
                disk_operations.remove_empty_directory(test_dir)
                assert test_dir.exists()
                
                # Remove file then directory
                disk_operations.unlink(test_file)
                disk_operations.remove_empty_directory(test_dir)
                assert not test_dir.exists()
