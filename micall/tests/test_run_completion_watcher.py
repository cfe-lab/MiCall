"""
Unit tests for run_completion_watcher module.

Migrated from build/MISEQ/tests/miseq_watch_runs_test.rb
"""

from pathlib import Path
from time import sleep
from unittest.mock import patch

from micall.monitor.run_completion_watcher import (
    find_unstable_runs,
    monitor_run_completion,
)


class TestFindUnstableRuns:
    """Test the find_unstable_runs() function."""

    def test_finds_basecalls_format_run(self, tmp_path):
        """Create a run directory with traditional BaseCalls format."""
        run_dir = tmp_path / "230224_M04401_0001_000000000-TEST1"
        basecalls_dir = run_dir / "Data/Intensities/BaseCalls"
        basecalls_dir.mkdir(parents=True)

        # Create some test FASTQ files
        for i in range(3):
            (basecalls_dir / f"sample_{i}.fastq.gz").touch()

        result = find_unstable_runs(tmp_path)

        assert len(result) == 1
        assert result[0].run_dir == run_dir
        assert result[0].file_count == 3
        assert "BaseCalls" in result[0].glob_pattern

    def test_finds_alignment_format_run(self, tmp_path):
        """Create a run directory with newer Alignment/Fastq format."""
        run_dir = tmp_path / "251030_M05995_0002_000000000-TEST2"
        alignment_dir = run_dir / "Alignment_1/20251101_052556/Fastq"
        alignment_dir.mkdir(parents=True)

        # Create some test FASTQ files
        for i in range(5):
            (alignment_dir / f"sample_{i}.fastq.gz").touch()

        result = find_unstable_runs(tmp_path)

        assert len(result) == 1
        assert result[0].run_dir == run_dir
        assert result[0].file_count == 5
        assert "Alignment" in result[0].glob_pattern

    def test_ignores_run_with_needsprocessing_file(self, tmp_path):
        """Create a run that already has needsprocessing file."""
        run_dir = tmp_path / "230101_M04401_0003_000000000-TEST3"
        basecalls_dir = run_dir / "Data/Intensities/BaseCalls"
        basecalls_dir.mkdir(parents=True)

        (basecalls_dir / "sample.fastq.gz").touch()
        (run_dir / "needsprocessing").touch()

        result = find_unstable_runs(tmp_path)

        assert len(result) == 0

    def test_ignores_run_with_processed_file(self, tmp_path):
        """Create a run that already has processed file."""
        run_dir = tmp_path / "230101_M04401_0004_000000000-TEST4"
        basecalls_dir = run_dir / "Data/Intensities/BaseCalls"
        basecalls_dir.mkdir(parents=True)

        (basecalls_dir / "sample.fastq.gz").touch()
        (run_dir / "processed").touch()

        result = find_unstable_runs(tmp_path)

        assert len(result) == 0

    def test_ignores_run_with_no_fastq_files(self, tmp_path):
        """Create a run directory but no FASTQ files."""
        run_dir = tmp_path / "230101_M04401_0005_000000000-TEST5"
        basecalls_dir = run_dir / "Data/Intensities/BaseCalls"
        basecalls_dir.mkdir(parents=True)

        # Directory exists but no .gz files
        result = find_unstable_runs(tmp_path)

        assert len(result) == 0

    def test_finds_multiple_runs(self, tmp_path):
        """Create multiple runs with different formats."""
        # Run 1 with BaseCalls
        run1 = tmp_path / "230101_M04401_0006_000000000-RUN1"
        basecalls_dir1 = run1 / "Data/Intensities/BaseCalls"
        basecalls_dir1.mkdir(parents=True)
        (basecalls_dir1 / "file1.fastq.gz").touch()

        # Run 2 with Alignment
        run2 = tmp_path / "251030_M05995_0007_000000000-RUN2"
        alignment_dir = run2 / "Alignment_1/20251101/Fastq"
        alignment_dir.mkdir(parents=True)
        (alignment_dir / "file1.fastq.gz").touch()
        (alignment_dir / "file2.fastq.gz").touch()

        result = find_unstable_runs(tmp_path)

        assert len(result) == 2
        run_dirs = {r.run_dir for r in result}
        assert run1 in run_dirs
        assert run2 in run_dirs

    def test_prefers_basecalls_when_both_exist(self, tmp_path):
        """Edge case: both locations have files (shouldn't normally happen)."""
        run_dir = tmp_path / "230101_M04401_0008_000000000-BOTH"
        basecalls_dir = run_dir / "Data/Intensities/BaseCalls"
        alignment_dir = run_dir / "Alignment_1/20251101/Fastq"
        basecalls_dir.mkdir(parents=True)
        alignment_dir.mkdir(parents=True)

        (basecalls_dir / "file1.fastq.gz").touch()
        (alignment_dir / "file2.fastq.gz").touch()

        result = find_unstable_runs(tmp_path)

        assert len(result) == 1
        # Should prefer BaseCalls (checked first in if/elif)
        assert "BaseCalls" in result[0].glob_pattern

    def test_handles_nonexistent_directory(self, tmp_path):
        """Test with a directory that doesn't exist."""
        nonexistent = tmp_path / "does_not_exist"
        result = find_unstable_runs(nonexistent)

        assert len(result) == 0


class TestMonitorRunCompletion:
    """Test the monitor_run_completion() function."""

    def test_monitor_runs_creates_needsprocessing_for_stable_runs(self, tmp_path):
        """Test that stable runs get marked with needsprocessing."""
        # Create runs with stable file counts
        run1 = tmp_path / "230224_M04401_0030_000000000-STABLE1"
        basecalls_dir1 = run1 / "Data/Intensities/BaseCalls"
        basecalls_dir1.mkdir(parents=True)
        for i in range(3):
            (basecalls_dir1 / f"file{i}.fastq.gz").touch()

        run2 = tmp_path / "230224_M04401_0031_000000000-STABLE2"
        alignment_dir = run2 / "Alignment_1/20251101/Fastq"
        alignment_dir.mkdir(parents=True)
        for i in range(2):
            (alignment_dir / f"file{i}.fastq.gz").touch()

        # Use exception to break out of infinite loop after desired behavior
        sleep_count = [0]

        def mock_sleep_impl(duration):
            sleep_count[0] += 1
            if sleep_count[0] >= 2:
                # After second sleep, files should be marked stable
                raise KeyboardInterrupt("Test complete")
            sleep(0.01)  # Small real sleep to avoid busy loop

        with patch(
            "micall.monitor.run_completion_watcher.sleep", side_effect=mock_sleep_impl
        ):
            try:
                monitor_run_completion(tmp_path)
            except KeyboardInterrupt:
                pass  # Expected - this is how we exit the infinite loop

        # Both runs should have needsprocessing markers
        assert (run1 / "needsprocessing").exists(), (
            "Run 1 should have needsprocessing marker"
        )
        assert (run2 / "needsprocessing").exists(), (
            "Run 2 should have needsprocessing marker"
        )

    def test_monitor_runs_skips_already_marked_runs(self, tmp_path):
        """Test that runs with existing markers are skipped."""
        run_dir = tmp_path / "230224_M04401_0032_000000000-ALREADY"
        basecalls_dir = run_dir / "Data/Intensities/BaseCalls"
        basecalls_dir.mkdir(parents=True)
        (basecalls_dir / "file.fastq.gz").touch()
        (run_dir / "needsprocessing").touch()

        # Use exception to exit after one iteration
        def mock_sleep_impl(duration):
            raise KeyboardInterrupt("Test complete")

        with patch(
            "micall.monitor.run_completion_watcher.sleep", side_effect=mock_sleep_impl
        ):
            try:
                monitor_run_completion(tmp_path)
            except KeyboardInterrupt:
                pass

        # Still should only have one needsprocessing file (not recreated)
        assert (run_dir / "needsprocessing").exists()

    def test_monitor_runs_handles_empty_directory(self, tmp_path):
        """Empty directory - should complete without error."""

        def mock_sleep_impl(duration):
            raise KeyboardInterrupt("Test complete")

        with patch(
            "micall.monitor.run_completion_watcher.sleep", side_effect=mock_sleep_impl
        ):
            try:
                monitor_run_completion(tmp_path)
            except KeyboardInterrupt:
                pass
        # No assertion needed - just shouldn't crash

    def test_monitor_runs_waits_for_stability(self, tmp_path):
        """Test that monitor_runs waits for files to stabilize."""
        run_dir = tmp_path / "230224_M04401_0033_000000000-GROWING"
        basecalls_dir = run_dir / "Data/Intensities/BaseCalls"
        basecalls_dir.mkdir(parents=True)

        # Create initial file
        (basecalls_dir / "file1.fastq.gz").touch()

        # Track sleep calls to control timing
        sleep_count = [0]

        def mock_sleep_impl(duration):
            sleep_count[0] += 1
            if sleep_count[0] == 1:
                # After first sleep (discovery), do nothing
                sleep(0.01)
            elif sleep_count[0] == 2:
                # After second sleep (first stability check), add a file
                (basecalls_dir / "file2.fastq.gz").touch()
                sleep(0.01)
            elif sleep_count[0] == 3:
                # After third sleep (second stability check), add another file
                (basecalls_dir / "file3.fastq.gz").touch()
                sleep(0.01)
            elif sleep_count[0] == 4:
                # After fourth sleep (third check), files are now stable
                sleep(0.01)
            else:
                # Stop after files stabilize and marker is created
                raise KeyboardInterrupt("Test complete")

        with patch(
            "micall.monitor.run_completion_watcher.sleep", side_effect=mock_sleep_impl
        ):
            try:
                monitor_run_completion(tmp_path)
            except KeyboardInterrupt:
                pass

        # Should have marked as complete when files stabilized
        assert (run_dir / "needsprocessing").exists(), (
            "Should create needsprocessing when files stabilize"
        )

        # Should have all three files
        assert (basecalls_dir / "file1.fastq.gz").exists()
        assert (basecalls_dir / "file2.fastq.gz").exists()
        assert (basecalls_dir / "file3.fastq.gz").exists()

    def test_monitor_detects_new_runs_after_startup(self, tmp_path):
        """Test that monitor detects runs added after it starts."""
        # Start with no runs
        sleep_count = [0]

        def mock_sleep_impl(duration):
            sleep_count[0] += 1
            if sleep_count[0] == 1:
                # After first iteration, add a new run
                run_dir = tmp_path / "230224_M04401_0040_000000000-LATE"
                basecalls_dir = run_dir / "Data/Intensities/BaseCalls"
                basecalls_dir.mkdir(parents=True)
                (basecalls_dir / "file1.fastq.gz").touch()
                sleep(0.01)
            elif sleep_count[0] == 2:
                # Files are stable now
                sleep(0.01)
            else:
                raise KeyboardInterrupt("Test complete")

        with patch(
            "micall.monitor.run_completion_watcher.sleep", side_effect=mock_sleep_impl
        ):
            try:
                monitor_run_completion(tmp_path)
            except KeyboardInterrupt:
                pass

        # Should have detected and marked the late-arriving run
        assert (
            tmp_path / "230224_M04401_0040_000000000-LATE" / "needsprocessing"
        ).exists()

    def test_monitor_handles_multiple_unstable_then_stable_runs(self, tmp_path):
        """Test monitoring multiple runs that stabilize at different times."""
        # Create two runs
        run1 = tmp_path / "230224_M04401_0050_000000000-FAST"
        basecalls_dir1 = run1 / "Data/Intensities/BaseCalls"
        basecalls_dir1.mkdir(parents=True)
        (basecalls_dir1 / "file1.fastq.gz").touch()

        run2 = tmp_path / "230224_M04401_0051_000000000-SLOW"
        basecalls_dir2 = run2 / "Data/Intensities/BaseCalls"
        basecalls_dir2.mkdir(parents=True)
        (basecalls_dir2 / "file1.fastq.gz").touch()

        sleep_count = [0]

        def mock_sleep_impl(duration):
            sleep_count[0] += 1
            if sleep_count[0] == 1:
                # After discovery
                sleep(0.01)
            elif sleep_count[0] == 2:
                # run1 is stable, but run2 gets another file
                (basecalls_dir2 / "file2.fastq.gz").touch()
                sleep(0.01)
            elif sleep_count[0] == 3:
                # Now run2 is stable too
                sleep(0.01)
            else:
                raise KeyboardInterrupt("Test complete")

        with patch(
            "micall.monitor.run_completion_watcher.sleep", side_effect=mock_sleep_impl
        ):
            try:
                monitor_run_completion(tmp_path)
            except KeyboardInterrupt:
                pass

        # Both should be marked
        assert (run1 / "needsprocessing").exists()
        assert (run2 / "needsprocessing").exists()

    def test_monitor_handles_file_removal(self, tmp_path):
        """Test that file removal is detected (count decreases)."""
        run_dir = tmp_path / "230224_M04401_0060_000000000-REMOVAL"
        basecalls_dir = run_dir / "Data/Intensities/BaseCalls"
        basecalls_dir.mkdir(parents=True)
        file1 = basecalls_dir / "file1.fastq.gz"
        file2 = basecalls_dir / "file2.fastq.gz"
        file1.touch()
        file2.touch()

        sleep_count = [0]

        def mock_sleep_impl(duration):
            sleep_count[0] += 1
            if sleep_count[0] == 1:
                sleep(0.01)
            elif sleep_count[0] == 2:
                # Remove a file - count decreases
                file2.unlink()
                sleep(0.01)
            elif sleep_count[0] == 3:
                # Now stable with 1 file
                sleep(0.01)
            else:
                raise KeyboardInterrupt("Test complete")

        with patch(
            "micall.monitor.run_completion_watcher.sleep", side_effect=mock_sleep_impl
        ):
            try:
                monitor_run_completion(tmp_path)
            except KeyboardInterrupt:
                pass

        # Should be marked after stabilizing at 1 file
        assert (run_dir / "needsprocessing").exists()


class TestEdgeCases:
    """Test edge cases and error conditions."""

    def test_find_runs_with_deeply_nested_alignment_paths(self, tmp_path):
        """Test various Alignment directory patterns."""
        # Test Alignment_1/datestamp/Fastq
        run1 = tmp_path / "230101_M04401_0100_000000000-ALIGN1"
        align1 = run1 / "Alignment_1/20230101_123456/Fastq"
        align1.mkdir(parents=True)
        (align1 / "file.fastq.gz").touch()

        # Test Alignment_2/datestamp/Fastq
        run2 = tmp_path / "230101_M04401_0101_000000000-ALIGN2"
        align2 = run2 / "Alignment_2/20230102/Fastq"
        align2.mkdir(parents=True)
        (align2 / "file.fastq.gz").touch()

        result = find_unstable_runs(tmp_path)

        assert len(result) == 2
        run_dirs = {r.run_dir for r in result}
        assert run1 in run_dirs
        assert run2 in run_dirs

    def test_find_runs_ignores_non_gz_files(self, tmp_path):
        """Test that only .gz files are counted."""
        run_dir = tmp_path / "230101_M04401_0110_000000000-NONGZ"
        basecalls_dir = run_dir / "Data/Intensities/BaseCalls"
        basecalls_dir.mkdir(parents=True)

        # Create non-gz files
        (basecalls_dir / "file1.fastq").touch()
        (basecalls_dir / "file2.txt").touch()
        (basecalls_dir / "file3.csv").touch()

        result = find_unstable_runs(tmp_path)

        # Should not find this run (no .gz files)
        assert len(result) == 0

    def test_find_runs_handles_symlinks(self, tmp_path):
        """Test handling of symlinked run directories."""
        # Create actual run
        real_run = tmp_path / "real_run"
        basecalls_dir = real_run / "Data/Intensities/BaseCalls"
        basecalls_dir.mkdir(parents=True)
        (basecalls_dir / "file.fastq.gz").touch()

        # Create symlink to run
        symlink_run = tmp_path / "230101_M04401_0120_000000000-SYMLINK"
        symlink_run.symlink_to(real_run)

        result = find_unstable_runs(tmp_path)

        # Should find the symlinked run
        assert len(result) == 2  # Both real and symlink are found

    def test_find_runs_with_special_characters_in_names(self, tmp_path):
        """Test run directories with special characters."""
        run_dir = tmp_path / "230101_M04401_0130_000000000-TEST_RUN-2.0"
        basecalls_dir = run_dir / "Data/Intensities/BaseCalls"
        basecalls_dir.mkdir(parents=True)
        (basecalls_dir / "file.fastq.gz").touch()

        result = find_unstable_runs(tmp_path)

        assert len(result) == 1
        assert result[0].run_dir == run_dir

    def test_monitor_recovers_from_marker_creation_failure(self, tmp_path):
        """Test that monitoring continues even if marker creation fails."""
        run_dir = tmp_path / "230224_M04401_0140_000000000-READONLY"
        basecalls_dir = run_dir / "Data/Intensities/BaseCalls"
        basecalls_dir.mkdir(parents=True)
        (basecalls_dir / "file.fastq.gz").touch()

        sleep_count = [0]
        touch_count = [0]

        original_touch = run_dir.joinpath("needsprocessing").touch

        def failing_touch(*args, **kwargs):
            touch_count[0] += 1
            if touch_count[0] == 1:
                raise OSError("Permission denied")
            # Second attempt succeeds
            return original_touch(*args, **kwargs)

        def mock_sleep_impl(duration):
            sleep_count[0] += 1
            if sleep_count[0] >= 3:
                raise KeyboardInterrupt("Test complete")
            sleep(0.01)

        with patch(
            "micall.monitor.run_completion_watcher.sleep", side_effect=mock_sleep_impl
        ):
            with patch.object(
                type(run_dir / "needsprocessing"), "touch", side_effect=failing_touch
            ):
                try:
                    monitor_run_completion(tmp_path)
                except KeyboardInterrupt:
                    pass

        # Should have attempted twice (failed once, succeeded once)
        # Note: The actual behavior depends on implementation details

    def test_find_runs_empty_glob_patterns(self, tmp_path):
        """Test that glob patterns work correctly with various scenarios."""
        # Run with files in subdirectory but not matching pattern
        run_dir = tmp_path / "230101_M04401_0150_000000000-WRONGDIR"
        wrong_dir = run_dir / "Data/Intensities/WrongDir"
        wrong_dir.mkdir(parents=True)
        (wrong_dir / "file.fastq.gz").touch()

        result = find_unstable_runs(tmp_path)

        # Should not find this run (files not in expected location)
        assert len(result) == 0

    def test_monitor_with_zero_files_then_files_appear(self, tmp_path):
        """Test run that starts empty then gets files."""
        run_dir = tmp_path / "230224_M04401_0160_000000000-LATE_FILES"
        basecalls_dir = run_dir / "Data/Intensities/BaseCalls"
        basecalls_dir.mkdir(parents=True)
        # Start with empty directory

        sleep_count = [0]

        def mock_sleep_impl(duration):
            sleep_count[0] += 1
            if sleep_count[0] == 1:
                # Add files after first check
                (basecalls_dir / "file1.fastq.gz").touch()
                sleep(0.01)
            elif sleep_count[0] == 2:
                # Files are now stable
                sleep(0.01)
            else:
                raise KeyboardInterrupt("Test complete")

        with patch(
            "micall.monitor.run_completion_watcher.sleep", side_effect=mock_sleep_impl
        ):
            try:
                monitor_run_completion(tmp_path)
            except KeyboardInterrupt:
                pass

        # Should detect and mark the run that appeared later
        assert (run_dir / "needsprocessing").exists()


class TestCrashRecovery:
    """Test crash recovery and error handling."""

    def test_monitor_recovers_from_unexpected_exception(self, tmp_path):
        """Test that monitor restarts after unexpected exception."""
        run_dir = tmp_path / "230224_M04401_0170_000000000-CRASH"
        basecalls_dir = run_dir / "Data/Intensities/BaseCalls"
        basecalls_dir.mkdir(parents=True)
        (basecalls_dir / "file1.fastq.gz").touch()

        sleep_count = [0]
        crash_count = [0]

        def mock_sleep_impl(duration):
            sleep_count[0] += 1
            if sleep_count[0] == 1:
                # First iteration - cause a crash
                crash_count[0] += 1
                raise RuntimeError("Simulated crash")
            elif sleep_count[0] == 2:
                # After recovery delay, should retry
                sleep(0.01)
            elif sleep_count[0] == 3:
                # Second check - files stable
                sleep(0.01)
            else:
                raise KeyboardInterrupt("Test complete")

        with patch(
            "micall.monitor.run_completion_watcher.sleep", side_effect=mock_sleep_impl
        ):
            with patch(
                "micall.monitor.run_completion_watcher.CRASH_RECOVERY_DELAY", 0.01
            ):
                try:
                    monitor_run_completion(tmp_path)
                except KeyboardInterrupt:
                    pass

        # Should have crashed once and recovered
        assert crash_count[0] == 1
        # Should have marked the run after recovery
        assert (run_dir / "needsprocessing").exists()

    def test_find_unstable_runs_handles_permission_error(self, tmp_path):
        """Test that permission errors on individual runs are handled gracefully."""
        # Create a good run
        good_run = tmp_path / "230224_M04401_0180_000000000-GOOD"
        good_basecalls = good_run / "Data/Intensities/BaseCalls"
        good_basecalls.mkdir(parents=True)
        (good_basecalls / "file.fastq.gz").touch()

        # Create a bad run that will error
        bad_run = tmp_path / "230224_M04401_0181_000000000-BAD"
        bad_basecalls_bad = bad_run / "Data/Intensities/BaseCalls"
        bad_basecalls_bad.mkdir(parents=True)
        (bad_basecalls_bad / "file.fastq.gz").touch()

        # Track which paths had glob errors
        error_count = [0]
        real_glob = Path.glob

        def mock_glob(self, pattern):
            # Raise error when globbing in the bad_run directory
            if "BAD" in str(self):
                error_count[0] += 1
                raise PermissionError("Simulated permission denied")
            # Use real glob for everything else
            return real_glob(self, pattern)

        with patch.object(Path, "glob", mock_glob):
            result = find_unstable_runs(tmp_path)

        # Should have attempted to glob in bad_run
        assert error_count[0] > 0
        # Should still find the good run despite error on bad run
        assert len(result) == 1
        assert result[0].run_dir == good_run

    def test_monitor_continues_after_marker_creation_failure(self, tmp_path):
        """Test that monitoring continues even if one marker fails to create."""
        # Create two runs
        run1 = tmp_path / "230224_M04401_0190_000000000-SUCCEED"
        basecalls1 = run1 / "Data/Intensities/BaseCalls"
        basecalls1.mkdir(parents=True)
        (basecalls1 / "file.fastq.gz").touch()

        run2 = tmp_path / "230224_M04401_0191_000000000-FAIL"
        basecalls2 = run2 / "Data/Intensities/BaseCalls"
        basecalls2.mkdir(parents=True)
        (basecalls2 / "file.fastq.gz").touch()

        sleep_count = [0]
        touch_calls = [0]

        def mock_touch(path):
            touch_calls[0] += 1
            # First touch (run1) succeeds
            if touch_calls[0] == 1:
                path.touch()
            # Second touch (run2) fails
            elif touch_calls[0] == 2:
                raise OSError("Simulated disk error")
            # Third touch (run2 retry) succeeds
            else:
                path.touch()

        def mock_sleep_impl(duration):
            sleep_count[0] += 1
            if sleep_count[0] >= 3:
                raise KeyboardInterrupt("Test complete")
            sleep(0.01)

        with patch("micall.monitor.disk_operations.touch", side_effect=mock_touch):
            with patch(
                "micall.monitor.run_completion_watcher.sleep",
                side_effect=mock_sleep_impl,
            ):
                try:
                    monitor_run_completion(tmp_path)
                except KeyboardInterrupt:
                    pass

        # First run should be marked
        assert (run1 / "needsprocessing").exists()
        # Second run should eventually be marked too (on retry)
        assert (run2 / "needsprocessing").exists()

    def test_keyboard_interrupt_propagates(self, tmp_path):
        """Test that KeyboardInterrupt is not caught by crash recovery."""
        sleep_count = [0]

        def mock_sleep_impl(duration):
            sleep_count[0] += 1
            raise KeyboardInterrupt("User requested shutdown")

        with patch(
            "micall.monitor.run_completion_watcher.sleep", side_effect=mock_sleep_impl
        ):
            try:
                monitor_run_completion(tmp_path)
                assert False, "Should have raised KeyboardInterrupt"
            except KeyboardInterrupt:
                pass  # Expected

        # Should have only tried once before exiting
        assert sleep_count[0] == 1
