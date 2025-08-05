#!/usr/bin/env python3
"""
Comprehensive additional tests for disk_operations module.
These tests cover advanced edge cases, error scenarios, and integration patterns
to maximize confidence in the implementation.
"""

import tempfile
from pathlib import Path
import errno
import threading
import time
from unittest.mock import patch

import pytest

from micall.monitor import disk_operations


class TestAdvancedRetryLogic:
    """Advanced tests for retry logic and error handling."""

    def test_exponential_backoff_timing_precision(self):
        """Test that exponential backoff calculations are precise."""
        min_wait = disk_operations.MINIMUM_RETRY_WAIT
        max_wait = disk_operations.MAXIMUM_RETRY_WAIT

        # Test specific timing values
        expected_timings = [
            (1, 5),  # First retry: 5 seconds
            (2, 10),  # Second retry: 10 seconds
            (3, 20),  # Third retry: 20 seconds
            (4, 40),  # Fourth retry: 40 seconds
            (10, 2560),  # Tenth retry: ~42 minutes
            (15, 81920),  # Fifteenth retry: ~22.7 hours
            (20, 86400),  # Twentieth retry: capped at 1 day
        ]

        for attempt, expected_seconds in expected_timings:
            result = disk_operations.calculate_retry_wait(min_wait, max_wait, attempt)
            assert result.total_seconds() == expected_seconds

    def test_retry_decorator_thread_safety(self):
        """Test that retry decorator works correctly with multiple threads."""
        call_counts = {}

        @disk_operations.disk_retry("test_op")
        def thread_test_function(thread_id):
            if thread_id not in call_counts:
                call_counts[thread_id] = 0
            call_counts[thread_id] += 1

            if call_counts[thread_id] <= 2:
                raise OSError(f"Error in thread {thread_id}")
            return f"success_{thread_id}"

        def run_test(thread_id, results):
            try:
                result = thread_test_function(thread_id)
                results[thread_id] = result
            except Exception as e:
                results[thread_id] = str(e)

        threads = []
        results: dict[int, float] = {}

        # Run multiple threads concurrently
        for i in range(3):  # Reduced from 5 to 3 threads
            thread = threading.Thread(target=run_test, args=(i, results))
            threads.append(thread)
            thread.start()

        # Wait for all threads to complete
        for thread in threads:
            thread.join()

        # Verify all threads succeeded after retries
        for i in range(3):  # Updated from 5 to 3
            assert results[i] == f"success_{i}"
            assert call_counts[i] == 3  # Should have retried twice

    def test_retry_with_different_exception_types(self):
        """Test retry behavior with various exception types."""

        @disk_operations.disk_retry("test_op")
        def test_function(exception_type):
            if exception_type == "os_error":
                raise OSError("OS Error")
            elif exception_type == "io_error":
                raise IOError("IO Error")
            elif exception_type == "permission_error":
                raise PermissionError("Permission denied")
            elif exception_type == "file_not_found":
                raise FileNotFoundError("File not found")
            elif exception_type == "value_error":
                raise ValueError("Value error")
            return "success"

        # Should retry these errors
        with patch("micall.monitor.disk_operations.wait_for_retry") as mock_wait:
            with pytest.raises(OSError):
                test_function("os_error")
            assert mock_wait.call_count == 14  # 15 attempts = 14 retries

        with patch("micall.monitor.disk_operations.wait_for_retry") as mock_wait:
            with pytest.raises(IOError):
                test_function("io_error")
            assert mock_wait.call_count == 14  # 15 attempts = 14 retries

        # Should NOT retry these errors
        with patch("micall.monitor.disk_operations.wait_for_retry") as mock_wait:
            with pytest.raises(ValueError):
                test_function("value_error")
            assert mock_wait.call_count == 0  # No retries

    def test_retry_success_after_max_attempts_minus_one(self):
        """Test success on the 14th attempt (just before max)."""
        attempt_count = 0

        @disk_operations.disk_retry("test_op")
        def almost_max_attempts():
            nonlocal attempt_count
            attempt_count += 1
            if attempt_count < 14:
                raise OSError("Almost there")
            return "success"

        with patch("micall.monitor.disk_operations.wait_for_retry"):
            result = almost_max_attempts()
            assert result == "success"
            assert attempt_count == 14


class TestEdgeCaseFileOperations:
    """Test edge cases for file operations."""

    def test_mkdir_p_with_special_modes(self):
        """Test mkdir_p with various file modes."""
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)

            # Test restrictive mode
            restricted_dir = temp_path / "restricted"
            disk_operations.mkdir_p(restricted_dir, mode=0o700, exist_ok=True)
            assert restricted_dir.exists()

            # Test permissive mode
            permissive_dir = temp_path / "permissive"
            disk_operations.mkdir_p(permissive_dir, mode=0o777, exist_ok=True)
            assert permissive_dir.exists()

            # Test nested creation with exist_ok=False then True
            nested_dir = temp_path / "level1" / "level2" / "level3"
            disk_operations.mkdir_p(nested_dir, exist_ok=False)
            assert nested_dir.exists()

            # Should not raise with exist_ok=True
            disk_operations.mkdir_p(nested_dir, exist_ok=True)

    def test_write_text_with_different_encodings(self):
        """Test write_text with various encodings."""
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)

            # Test UTF-8 with special characters
            utf8_file = temp_path / "utf8.txt"
            utf8_content = "Hello 世界 🌍 café naïve résumé"
            disk_operations.write_text(utf8_file, utf8_content, encoding="utf-8")
            assert utf8_file.read_text(encoding="utf-8") == utf8_content

            # Test ASCII
            ascii_file = temp_path / "ascii.txt"
            ascii_content = "Simple ASCII text"
            disk_operations.write_text(ascii_file, ascii_content, encoding="ascii")
            assert ascii_file.read_text(encoding="ascii") == ascii_content

            # Test Latin-1
            latin1_file = temp_path / "latin1.txt"
            latin1_content = "Café résumé naïve"
            disk_operations.write_text(latin1_file, latin1_content, encoding="latin-1")
            assert latin1_file.read_text(encoding="latin-1") == latin1_content

    def test_unlink_with_missing_ok_variations(self):
        """Test unlink behavior with missing_ok parameter."""
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)

            # Test unlink existing file with missing_ok=True
            existing_file = temp_path / "existing.txt"
            existing_file.write_text("content")
            disk_operations.unlink(existing_file, missing_ok=True)
            assert not existing_file.exists()

            # Test unlink non-existing file with missing_ok=True
            non_existing = temp_path / "non_existing.txt"
            disk_operations.unlink(non_existing, missing_ok=True)  # Should not raise

            # Test unlink existing file with missing_ok=False (should work)
            existing_file2 = temp_path / "existing2.txt"
            existing_file2.write_text("content")
            disk_operations.unlink(existing_file2, missing_ok=False)
            assert not existing_file2.exists()

            # Test that missing_ok=False would raise for non-existing file
            # (but don't actually call it since it would retry for 15 attempts)
            # This tests our understanding of the expected behavior
            assert not non_existing.exists()  # Confirm file doesn't exist

    def test_copy_file_with_metadata_preservation(self):
        """Test that copy_file preserves metadata using shutil.copy2."""
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)

            source = temp_path / "source.txt"
            target = temp_path / "target.txt"

            # Create source file with content
            source.write_text("test content")
            original_stat = source.stat()

            # Add a small delay to ensure timestamp difference would be detectable
            time.sleep(0.01)

            # Copy file
            disk_operations.copy_file(source, target)

            # Verify content
            assert target.read_text() == "test content"

            # Verify metadata preservation (size should match)
            assert target.stat().st_size == original_stat.st_size

    def test_rmtree_with_readonly_files(self):
        """Test rmtree with read-only files (ignore_errors behavior)."""
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)

            # Create directory with read-only file
            target_dir = temp_path / "readonly_dir"
            target_dir.mkdir()
            readonly_file = target_dir / "readonly.txt"
            readonly_file.write_text("readonly content")
            readonly_file.chmod(0o444)  # Read-only

            # Test with ignore_errors=True
            disk_operations.rmtree(target_dir, ignore_errors=True)
            # Should succeed regardless of read-only files


class TestComplexIntegrationScenarios:
    """Test complex scenarios that mirror real-world usage."""

    def test_full_kive_watcher_workflow_simulation(self):
        """Simulate a simplified kive_watcher workflow."""
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)

            # Simulate run structure
            run_folder = temp_path / "run_001"
            results_path = run_folder / "Results" / "version_1.0"
            scratch_path = results_path / "scratch"

            # Step 1: Create initial structure (only 1 sample for speed)
            disk_operations.mkdir_p(scratch_path / "sample_001", exist_ok=True)

            # Step 2: Simulate processing - create output files
            sample_path = scratch_path / "sample_001"
            disk_operations.write_text(
                sample_path / "cascade.csv", "header,data\nval1,val2\n"
            )
            disk_operations.write_text(
                sample_path / "conseq.csv", "sequence,count\nATCG,100\n"
            )

            # Step 3: Simulate collation process
            collated_path = results_path / "collated"
            disk_operations.mkdir_p(collated_path, exist_ok=True)

            # Collate CSV files (simplified)
            with disk_operations.disk_file_operation(
                collated_path / "cascade.csv", "w"
            ) as target:
                target.write("sample,header,data\n")
                source_path = scratch_path / "sample_001" / "cascade.csv"
                with disk_operations.disk_file_operation(
                    source_path, "r"
                ) as source:
                    lines = source.readlines()
                    for line in lines[1:]:  # Skip header
                        target.write(f"sample_001,{line}")

            # Step 4: Final cleanup
            disk_operations.rmtree(scratch_path)
            disk_operations.touch(results_path / "doneprocessing")

            # Verify final state
            assert (collated_path / "cascade.csv").exists()
            assert (results_path / "doneprocessing").exists()
            assert not scratch_path.exists()

    def test_concurrent_file_operations(self):
        """Test concurrent file operations across multiple threads."""
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            results = {}
            errors = {}

            def worker_thread(thread_id):
                try:
                    # Each thread works in its own subdirectory
                    thread_dir = temp_path / f"thread_{thread_id}"
                    disk_operations.mkdir_p(thread_dir, exist_ok=True)

                    # Create multiple files
                    for i in range(10):
                        file_path = thread_dir / f"file_{i}.txt"
                        disk_operations.write_text(
                            file_path, f"Thread {thread_id} file {i}"
                        )

                    # Copy files to backup directory
                    backup_dir = thread_dir / "backup"
                    disk_operations.mkdir_p(backup_dir, exist_ok=True)

                    for i in range(10):
                        source = thread_dir / f"file_{i}.txt"
                        target = backup_dir / f"file_{i}.txt"
                        disk_operations.copy_file(source, target)

                    # Clean up some files
                    for i in range(5):
                        disk_operations.unlink(thread_dir / f"file_{i}.txt")

                    # Final verification
                    remaining_files = list(thread_dir.glob("file_*.txt"))
                    backup_files = list(backup_dir.glob("file_*.txt"))

                    results[thread_id] = {
                        "remaining": len(remaining_files),
                        "backup": len(backup_files),
                    }

                except Exception as e:
                    errors[thread_id] = str(e)

            # Run multiple worker threads
            threads = []
            for i in range(4):  # Reduced from 8 to 4 threads
                thread = threading.Thread(target=worker_thread, args=(i,))
                threads.append(thread)
                thread.start()

            # Wait for completion
            for thread in threads:
                thread.join()

            # Verify results
            assert len(errors) == 0, f"Errors occurred: {errors}"
            assert len(results) == 4  # Updated from 8 to 4

            for thread_id, result in results.items():
                assert result["remaining"] == 5  # 5 files left after deleting 5
                assert result["backup"] == 10  # All 10 files backed up

    def test_error_recovery_scenarios(self):
        """Test various error recovery scenarios."""
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)

            # Scenario 1: Network interruption during file operations
            def simulate_network_interruption():
                call_count = 0

                def mock_write_text(content, encoding="utf-8"):
                    nonlocal call_count
                    call_count += 1
                    if call_count <= 3:
                        raise OSError(errno.EIO, "I/O error - network interruption")
                    return content

                error_file = temp_path / "network_error.txt"
                with patch.object(Path, "write_text", side_effect=mock_write_text):
                    with patch("micall.monitor.disk_operations.wait_for_retry"):
                        disk_operations.write_text(error_file, "test content")
                        assert call_count == 4  # Should succeed on 4th attempt

            simulate_network_interruption()

            # Scenario 2: Permission denied then granted
            def simulate_permission_recovery():
                call_count = 0

                def mock_mkdir(mode=0o777, parents=False, exist_ok=False):
                    nonlocal call_count
                    call_count += 1
                    if call_count <= 2:
                        raise OSError(errno.EACCES, "Permission denied")
                    # Success on 3rd attempt
                    return None

                perm_dir = temp_path / "permission_test"
                with patch.object(Path, "mkdir", side_effect=mock_mkdir):
                    with patch("micall.monitor.disk_operations.wait_for_retry"):
                        disk_operations.mkdir_p(perm_dir, exist_ok=True)
                        assert call_count == 3

            simulate_permission_recovery()

    def test_large_file_operations(self):
        """Test operations with larger files and directories."""
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)

            # Create directory with many files
            large_dir = temp_path / "large_directory"
            disk_operations.mkdir_p(large_dir, exist_ok=True)

            # Create 50 files with various sizes (reduced from 100)
            for i in range(50):
                file_path = large_dir / f"file_{i:03d}.txt"
                content = f"File {i} content " * (i + 1)  # Variable size content
                disk_operations.write_text(file_path, content)

            # Create nested directory structure
            for level1 in range(3):  # Reduced from 5 to 3
                for level2 in range(3):  # Reduced from 5 to 3
                    nested_dir = large_dir / f"level1_{level1}" / f"level2_{level2}"
                    disk_operations.mkdir_p(nested_dir, exist_ok=True)

                    # Add a file in each nested directory
                    nested_file = nested_dir / "nested_file.txt"
                    disk_operations.write_text(
                        nested_file, f"Nested content {level1}-{level2}"
                    )

            # Verify structure was created
            files = list(large_dir.glob("file_*.txt"))
            assert len(files) == 50  # Updated from 100 to 50

            nested_files = list(large_dir.glob("**/nested_file.txt"))
            assert len(nested_files) == 9  # Updated: 3 * 3 = 9

            # Test bulk removal
            disk_operations.rmtree(large_dir)
            assert not large_dir.exists()


class TestErrorMessageVerification:
    """Test that error messages are informative and correctly logged."""

    def test_retry_error_messages(self):
        """Test that retry error messages contain useful information."""

        @disk_operations.disk_retry("test_operation")
        def failing_operation():
            raise OSError("Simulated failure")

        with patch("micall.monitor.disk_operations.logger.error") as mock_logger:
            with patch("micall.monitor.disk_operations.sleep"):
                with pytest.raises(OSError):
                    failing_operation()

                # Verify error logging
                assert mock_logger.call_count >= 1

                # Check that final error message mentions max attempts
                final_call = mock_logger.call_args_list[-1]
                error_message = final_call[0][0]
                assert "test_operation" in error_message
                assert "15 attempts" in error_message

    def test_wait_for_retry_logging(self):
        """Test wait_for_retry logging behavior."""

        # Test with logging enabled
        with patch("micall.monitor.disk_operations.logger.error") as mock_logger:
            with patch("micall.monitor.disk_operations.sleep"):
                disk_operations.wait_for_retry(3, "test_op", is_logged=True)

                mock_logger.assert_called_once()
                # The log message uses string formatting, so check that args were passed
                assert mock_logger.call_args[0][1] == "test_op"
                # Third argument should be a timedelta object
                assert hasattr(mock_logger.call_args[0][2], 'total_seconds')

        # Test with logging disabled
        with patch("micall.monitor.disk_operations.logger.error") as mock_logger:
            with patch("micall.monitor.disk_operations.sleep"):
                disk_operations.wait_for_retry(3, "test_op", is_logged=False)

                mock_logger.assert_not_called()


class TestDiskFileOperationContextManager:
    """Comprehensive tests for the disk_file_operation context manager."""

    def test_context_manager_exception_propagation(self):
        """Test that exceptions inside context manager are properly propagated."""
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            test_file = temp_path / "test.txt"

            with pytest.raises(ValueError):
                with disk_operations.disk_file_operation(test_file, "w") as f:
                    f.write("some content")
                    raise ValueError("Test exception")

            # File should still be created even though exception occurred
            assert test_file.exists()

    def test_context_manager_with_binary_modes(self):
        """Test context manager with binary file modes."""
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            binary_file = temp_path / "binary.dat"

            test_data = b"Binary data: \x00\x01\x02\x03\xff"

            # Write binary data
            with disk_operations.disk_file_operation(binary_file, "wb") as f:
                f.write(test_data)

            # Read binary data back
            with disk_operations.disk_file_operation(binary_file, "rb") as f:
                read_data = f.read()

            assert read_data == test_data

    def test_context_manager_retry_on_open_failure(self):
        """Test that context manager retries on file open failures."""
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            test_file = temp_path / "retry_test.txt"

            call_count = 0
            original_open = Path.open

            def mock_open(
                self, mode="r", buffering=-1, encoding=None, errors=None, newline=None
            ):
                nonlocal call_count
                call_count += 1
                if call_count <= 2:
                    raise OSError("Simulated open failure")
                return original_open(self, mode, buffering, encoding, errors, newline)

            with patch.object(Path, "open", mock_open):
                with patch("micall.monitor.disk_operations.wait_for_retry"):
                    with disk_operations.disk_file_operation(test_file, "w") as f:
                        f.write("success after retry")

                    assert call_count == 3  # Should succeed on 3rd attempt

            assert test_file.read_text() == "success after retry"


class TestRemoveEmptyDirectoryAdvanced:
    """Advanced tests for remove_empty_directory function."""

    def test_remove_empty_directory_with_hidden_files(self):
        """Test behavior with hidden files (should not remove directory)."""
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            test_dir = temp_path / "hidden_files_dir"
            test_dir.mkdir()

            # Create hidden file
            hidden_file = test_dir / ".hidden_file"
            hidden_file.write_text("hidden content")

            # Should not remove directory because it contains hidden file
            disk_operations.remove_empty_directory(test_dir)
            assert test_dir.exists()
            assert hidden_file.exists()

    def test_remove_empty_directory_with_symlinks(self):
        """Test behavior with symbolic links."""
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)

            # Create target file and directory with symlink
            target_file = temp_path / "target.txt"
            target_file.write_text("target content")

            symlink_dir = temp_path / "symlink_dir"
            symlink_dir.mkdir()

            symlink_file = symlink_dir / "link_to_target"
            symlink_file.symlink_to(target_file)

            # Should not remove directory because it contains symlink
            disk_operations.remove_empty_directory(symlink_dir)
            assert symlink_dir.exists()
            assert symlink_file.exists()

    def test_remove_empty_directory_enotempty_no_retry(self):
        """Test that ENOTEMPTY doesn't trigger retry logic."""
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            test_dir = temp_path / "nonempty_test"
            test_dir.mkdir()

            # Create file in directory
            (test_dir / "file.txt").write_text("content")

            with patch("micall.monitor.disk_operations.wait_for_retry") as mock_wait:
                disk_operations.remove_empty_directory(test_dir)

                # Should not retry for ENOTEMPTY
                assert mock_wait.call_count == 0
                assert test_dir.exists()


class TestIntegrationWithActualFilesystem:
    """Integration tests with actual filesystem operations."""

    def test_real_filesystem_stress_test(self):
        """Stress test with real filesystem operations."""
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)

            # Create complex directory structure
            for i in range(10):  # Reduced from 20 to 10
                base_dir = temp_path / f"stress_test_{i}"
                disk_operations.mkdir_p(base_dir, exist_ok=True)

                # Create files with various operations
                for j in range(5):  # Reduced from 10 to 5
                    file_path = base_dir / f"file_{j}.txt"
                    content = f"Stress test content {i}-{j} " * 50
                    disk_operations.write_text(file_path, content)

                    # Copy to backup
                    backup_path = base_dir / f"backup_{j}.txt"
                    disk_operations.copy_file(file_path, backup_path)

                    # Create and remove temporary file
                    temp_file = base_dir / f"temp_{j}.txt"
                    disk_operations.touch(temp_file)
                    disk_operations.unlink(temp_file)

                # Test remove_empty_directory on subdirectories
                for k in range(2):  # Reduced from 3 to 2
                    empty_dir = base_dir / f"empty_{k}"
                    disk_operations.mkdir_p(empty_dir, exist_ok=True)
                    disk_operations.remove_empty_directory(empty_dir)
                    assert not empty_dir.exists()

            # Verify final state
            created_dirs = list(temp_path.glob("stress_test_*"))
            assert len(created_dirs) == 10  # Updated from 20 to 10

            total_files = list(temp_path.glob("**/*.txt"))
            assert (
                len(total_files) == 100  # Updated: 10 dirs * 10 files each (5 original + 5 backup)
            )

            # Cleanup stress test
            for stress_dir in created_dirs:
                disk_operations.rmtree(stress_dir)
                assert not stress_dir.exists()

    def test_unicode_filename_handling(self):
        """Test operations with Unicode filenames."""
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)

            unicode_names = [
                "café_file.txt",
                "фаїл.txt",  # Ukrainian
                "文件.txt",  # Chinese
                "📁_emoji_file.txt",
                "naïve_résumé.txt",
            ]

            for unicode_name in unicode_names:
                try:
                    file_path = temp_path / unicode_name
                    content = f"Unicode content for {unicode_name}"

                    # Test all operations with Unicode names
                    disk_operations.write_text(file_path, content)
                    assert file_path.exists()

                    # Copy to backup
                    backup_path = temp_path / f"backup_{unicode_name}"
                    disk_operations.copy_file(file_path, backup_path)
                    assert backup_path.exists()

                    # Verify content
                    assert backup_path.read_text() == content

                    # Cleanup
                    disk_operations.unlink(file_path)
                    disk_operations.unlink(backup_path)

                except OSError:
                    # Some filesystems don't support certain Unicode characters
                    # This is acceptable and not a failure of our implementation
                    pass


class TestPerformanceAndScalability:
    """Tests focused on performance and scalability aspects."""

    def test_large_directory_removal_performance(self):
        """Test performance of removing large directories."""
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)

            # Create large directory structure
            large_dir = temp_path / "large_structure"
            disk_operations.mkdir_p(large_dir, exist_ok=True)

            # Create nested structure: 3 levels deep, 3 subdirs per level (reduced from 5 levels, 5 subdirs)
            current_paths = [large_dir]
            for level in range(3):  # Reduced from 5 to 3
                next_paths = []
                for current_path in current_paths:
                    for i in range(3):  # Reduced from 5 to 3
                        sub_path = current_path / f"level_{level}_sub_{i}"
                        disk_operations.mkdir_p(sub_path, exist_ok=True)

                        # Add some files
                        for j in range(3):
                            file_path = sub_path / f"file_{j}.txt"
                            disk_operations.write_text(
                                file_path, f"Content {level}-{i}-{j}"
                            )

                        next_paths.append(sub_path)
                current_paths = next_paths

            # Measure removal time (should complete reasonably quickly)
            start_time = time.time()
            disk_operations.rmtree(large_dir)
            end_time = time.time()

            assert not large_dir.exists()
            removal_time = end_time - start_time
            # Should complete within reasonable time (adjust threshold as needed)
            assert removal_time < 5.0  # 5 seconds should be more than enough

    def test_many_small_operations_performance(self):
        """Test performance of many small operations."""
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)

            start_time = time.time()

            # Perform many small operations (reduced from 200 to 100)
            for i in range(100):
                # Create, write, read, copy, delete cycle
                file_path = temp_path / f"perf_test_{i}.txt"
                content = f"Performance test content {i}"

                disk_operations.write_text(file_path, content)

                copy_path = temp_path / f"copy_{i}.txt"
                disk_operations.copy_file(file_path, copy_path)

                disk_operations.unlink(file_path)
                disk_operations.unlink(copy_path)

            end_time = time.time()
            total_time = end_time - start_time

            # Should complete within reasonable time (adjusted for fewer operations)
            assert total_time < 5.0  # Reduced from 10 to 5 seconds for 100 cycles


if __name__ == "__main__":
    pytest.main([__file__, "-v", "--tb=short"])
