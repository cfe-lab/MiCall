"""
Test for time-based retry logging behavior.

Tests that retry operatio                with patch('micall.monitor.kive_watcher.logger') as mock_logger:
            with patch('time.time', return_value=self.current_time.timestamp()):
                # Should not log during first hour
                network_wait_for_retry(1, self.start_time)
                mock_logger.error.assert_not_called()only log after one hour has passed since first failure.
"""

import unittest
from unittest.mock import patch
from datetime import datetime, timedelta

from micall.monitor.disk_operations import wait_for_retry as disk_wait_for_retry
from micall.monitor.kive_watcher import wait_for_retry as network_wait_for_retry


class TestRetryLogging(unittest.TestCase):
    """Test that retry operations use time-based logging."""

    def setUp(self):
        """Set up test fixtures."""
        self.logger_patcher = patch("micall.monitor.disk_operations.logger")
        self.mock_logger = self.logger_patcher.start()

        self.sleep_patcher = patch("time.sleep")
        self.mock_sleep = self.sleep_patcher.start()

        # Also patch sleep in the specific modules
        self.disk_sleep_patcher = patch("micall.monitor.disk_operations.sleep")
        self.disk_sleep_patcher.start()

        # Use a fixed time for testing
        self.start_time = datetime(2024, 1, 1, 12, 0, 0)

    def tearDown(self):
        """Clean up test fixtures."""
        self.logger_patcher.stop()
        self.sleep_patcher.stop()
        self.disk_sleep_patcher.stop()

    def test_disk_retry_logs_after_one_hour(self):
        """Test that disk retries only log after one hour."""
        # Time just started (within first hour)
        current_time = self.start_time + timedelta(minutes=30)

        with patch("micall.monitor.disk_operations.datetime") as mock_datetime:
            mock_datetime.now.return_value = current_time

            # Should not log during first hour
            disk_wait_for_retry(1, "test_operation", self.start_time)
            self.mock_logger.error.assert_not_called()

        # Time after one hour
        current_time = self.start_time + timedelta(hours=1, minutes=5)

        with patch("micall.monitor.disk_operations.datetime") as mock_datetime:
            mock_datetime.now.return_value = current_time

            # Should log on second retry after 1 hour
            disk_wait_for_retry(2, "test_operation", self.start_time)
            self.mock_logger.error.assert_called_once()

    def test_network_retry_logs_after_one_hour(self):
        """Test that network retries only log after one hour."""
        network_logger_patcher = patch("micall.monitor.kive_watcher.logger")
        network_sleep_patcher = patch("micall.monitor.kive_watcher.sleep")
        mock_network_logger = network_logger_patcher.start()
        network_sleep_patcher.start()

        try:
            # Time just started (within first hour)
            current_time = self.start_time + timedelta(minutes=45)

            with patch("micall.monitor.kive_watcher.datetime") as mock_datetime:
                mock_datetime.now.return_value = current_time

                # Should not log during first hour
                network_wait_for_retry(1, self.start_time)
                mock_network_logger.warning.assert_not_called()

            # Time after one hour
            current_time = self.start_time + timedelta(hours=1, minutes=10)

            with patch("micall.monitor.kive_watcher.datetime") as mock_datetime:
                mock_datetime.now.return_value = current_time

                # Should log after one hour
                network_wait_for_retry(2, self.start_time)
                mock_network_logger.warning.assert_called_once()

        finally:
            network_logger_patcher.stop()
            network_sleep_patcher.stop()

    def test_disk_retry_with_old_start_time_always_logs(self):
        """Test that disk retries with very old start_time always log."""
        # Use a start time from 2 hours ago
        old_start_time = self.start_time - timedelta(hours=2)
        disk_wait_for_retry(1, "test_operation", old_start_time)
        self.mock_logger.error.assert_called_once()

    def test_network_retry_with_old_start_time_always_logs(self):
        """Test that network retries with very old start_time always log."""
        network_logger_patcher = patch("micall.monitor.kive_watcher.logger")
        network_sleep_patcher = patch("micall.monitor.kive_watcher.sleep")
        mock_network_logger = network_logger_patcher.start()
        network_sleep_patcher.start()

        try:
            # Use a start time from 2 hours ago
            old_start_time = self.start_time - timedelta(hours=2)
            network_wait_for_retry(1, old_start_time)
            mock_network_logger.warning.assert_called_once()

        finally:
            network_logger_patcher.stop()
            network_sleep_patcher.stop()

    def test_retry_timing_behavior(self):
        """Test that retries check elapsed time correctly."""
        # Test disk operations - should not log within first hour
        current_time = self.start_time + timedelta(minutes=30)
        with patch("micall.monitor.disk_operations.datetime") as mock_datetime:
            mock_datetime.now.return_value = current_time
            disk_wait_for_retry(1, "test_operation", self.start_time)
            self.mock_logger.error.assert_not_called()

        # Test after one hour - should log
        current_time = self.start_time + timedelta(hours=2)
        with patch("micall.monitor.disk_operations.datetime") as mock_datetime:
            mock_datetime.now.return_value = current_time
            disk_wait_for_retry(2, "test_operation", self.start_time)
            self.mock_logger.error.assert_called_once()

        # Reset mock for network operations
        self.mock_logger.reset_mock()

        # Test network operations - should not log within first hour
        with (
            patch("micall.monitor.kive_watcher.logger") as mock_network_logger,
            patch("micall.monitor.kive_watcher.sleep"),
        ):
            current_time = self.start_time + timedelta(minutes=30)
            with patch("micall.monitor.kive_watcher.datetime") as mock_datetime:
                mock_datetime.now.return_value = current_time
                network_wait_for_retry(1, self.start_time)
                mock_network_logger.warning.assert_not_called()

            with patch("micall.monitor.kive_watcher.datetime") as mock_datetime:
                current_time = self.start_time + timedelta(hours=2)
                mock_datetime.now.return_value = current_time
                network_wait_for_retry(2, self.start_time)
                mock_network_logger.warning.assert_called_once()
