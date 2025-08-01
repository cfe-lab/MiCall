"""
Tests for disk_operations module.

Tests the low-level disk operation wrappers with network-style retry logic.
"""
import tempfile
from pathlib import Path
from unittest.mock import patch, Mock

import pytest

from micall.monitor.disk_operations import (
    calculate_retry_wait, wait_for_retry, disk_retry,
    mkdir_p, rmtree, move, copy_file, write_text, touch, unlink, rmdir, rename,
    copy_fileobj, disk_file_operation, remove_empty_directory,
    MINIMUM_RETRY_WAIT, MAXIMUM_RETRY_WAIT
)


class TestRetryLogic:
    """Test the exponential backoff and retry mechanisms."""

    def test_calculate_retry_wait_first_attempt(self):
        """First attempt should use minimum wait time."""
        result = calculate_retry_wait(MINIMUM_RETRY_WAIT, MAXIMUM_RETRY_WAIT, 1)
        assert result == MINIMUM_RETRY_WAIT

    def test_calculate_retry_wait_exponential_growth(self):
        """Wait time should grow exponentially."""
        wait1 = calculate_retry_wait(MINIMUM_RETRY_WAIT, MAXIMUM_RETRY_WAIT, 1)
        wait2 = calculate_retry_wait(MINIMUM_RETRY_WAIT, MAXIMUM_RETRY_WAIT, 2)
        wait3 = calculate_retry_wait(MINIMUM_RETRY_WAIT, MAXIMUM_RETRY_WAIT, 3)
        
        assert wait2.total_seconds() == wait1.total_seconds() * 2
        assert wait3.total_seconds() == wait1.total_seconds() * 4

    def test_calculate_retry_wait_max_cap(self):
        """Wait time should be capped at maximum."""
        result = calculate_retry_wait(MINIMUM_RETRY_WAIT, MAXIMUM_RETRY_WAIT, 30)
        assert result == MAXIMUM_RETRY_WAIT

    @patch('micall.monitor.disk_operations.sleep')
    @patch('micall.monitor.disk_operations.logger')
    def test_wait_for_retry_logs_error(self, mock_logger, mock_sleep):
        """wait_for_retry should log error and sleep."""
        wait_for_retry(3, "test operation")
        
        mock_logger.error.assert_called_once()
        mock_sleep.assert_called_once()
        args, _ = mock_logger.error.call_args
        # Check that the log message contains the operation name and "waiting"
        assert "test operation" in args[1]
        assert "waiting" in args[0]  # The format string contains "waiting"

    @patch('micall.monitor.disk_operations.sleep')
    @patch('micall.monitor.disk_operations.logger')
    def test_wait_for_retry_no_logging(self, mock_logger, mock_sleep):
        """wait_for_retry with is_logged=False should not log."""
        wait_for_retry(2, "test operation", is_logged=False)
        
        mock_logger.error.assert_not_called()
        mock_sleep.assert_called_once()


class TestDiskRetryDecorator:
    """Test the @disk_retry decorator."""

    def test_decorator_success_first_try(self):
        """Function succeeding on first attempt should not retry."""
        mock_func = Mock(return_value="success")
        decorated = disk_retry("test_op")(mock_func)
        
        result = decorated("arg1", kwarg1="value1")
        
        assert result == "success"
        mock_func.assert_called_once_with("arg1", kwarg1="value1")

    @patch('micall.monitor.disk_operations.wait_for_retry')
    def test_decorator_retry_on_os_error(self, mock_wait):
        """OSError should trigger retry mechanism."""
        mock_func = Mock(side_effect=[OSError("Disk error"), "success"])
        decorated = disk_retry("test_op")(mock_func)
        
        result = decorated("arg1")
        
        assert result == "success"
        assert mock_func.call_count == 2
        mock_wait.assert_called_once_with(1, "test_op", False)

    @patch('micall.monitor.disk_operations.wait_for_retry')
    def test_decorator_retry_on_io_error(self, mock_wait):
        """IOError should trigger retry mechanism."""
        mock_func = Mock(side_effect=[IOError("I/O error"), "success"])
        decorated = disk_retry("test_op")(mock_func)
        
        result = decorated()
        
        assert result == "success"
        assert mock_func.call_count == 2
        mock_wait.assert_called_once_with(1, "test_op", False)

    @patch('micall.monitor.disk_operations.wait_for_retry')
    def test_decorator_max_attempts_exceeded(self, mock_wait):
        """Should raise exception after max attempts."""
        error = OSError("Persistent disk error")
        mock_func = Mock(side_effect=error)
        decorated = disk_retry("test_op")(mock_func)
        
        with pytest.raises(OSError, match="Persistent disk error"):
            decorated()
        
        assert mock_func.call_count == 15  # max attempts
        assert mock_wait.call_count == 14  # one less than max attempts

    def test_decorator_non_retryable_error(self):
        """Non-disk errors should not trigger retry."""
        error = ValueError("Not a disk error")
        mock_func = Mock(side_effect=error)
        decorated = disk_retry("test_op")(mock_func)
        
        with pytest.raises(ValueError, match="Not a disk error"):
            decorated()
        
        mock_func.assert_called_once()


class TestFileSystemOperations:
    """Test individual file system operation wrappers."""

    @patch('pathlib.Path.mkdir')
    def test_mkdir_p_success(self, mock_mkdir):
        """mkdir_p should call Path.mkdir with correct arguments."""
        path = Path("/test/dir")
        
        mkdir_p(path, mode=0o755, parents=True, exist_ok=True)
        
        mock_mkdir.assert_called_once_with(mode=0o755, parents=True, exist_ok=True)

    @patch('pathlib.Path.mkdir')
    @patch('micall.monitor.disk_operations.wait_for_retry')
    def test_mkdir_p_retry_on_failure(self, mock_wait, mock_mkdir):
        """mkdir_p should retry on OSError."""
        path = Path("/test/dir")
        mock_mkdir.side_effect = [OSError("Permission denied"), None]
        
        mkdir_p(path)
        
        assert mock_mkdir.call_count == 2
        mock_wait.assert_called_once_with(1, "mkdir", False)

    @patch('shutil.rmtree')
    @patch('pathlib.Path.exists')
    def test_rmtree_success(self, mock_exists, mock_rmtree):
        """rmtree should call shutil.rmtree with correct arguments."""
        path = Path("/test/dir")
        mock_exists.return_value = True
        
        rmtree(path, ignore_errors=True)
        
        mock_rmtree.assert_called_once_with(path, ignore_errors=True)

    @patch('shutil.rmtree')
    @patch('pathlib.Path.exists')
    @patch('micall.monitor.disk_operations.wait_for_retry')
    def test_rmtree_retry_on_failure(self, mock_wait, mock_exists, mock_rmtree):
        """rmtree should retry on OSError."""
        path = Path("/test/dir")
        mock_exists.return_value = True
        mock_rmtree.side_effect = [OSError("Directory not empty"), None]
        
        rmtree(path)
        
        assert mock_rmtree.call_count == 2
        mock_wait.assert_called_once_with(1, "rmtree", False)

    @patch('shutil.move')
    def test_move_success(self, mock_move):
        """move should call shutil.move."""
        src = Path("/src/file")
        dst = Path("/dst/file")
        
        move(src, dst)
        
        mock_move.assert_called_once_with(str(src), str(dst))

    @patch('shutil.copy2')
    def test_copy_file_success(self, mock_copy2):
        """copy_file should call shutil.copy2."""
        src = Path("/src/file")
        dst = Path("/dst/file")
        
        copy_file(src, dst)
        
        mock_copy2.assert_called_once_with(str(src), str(dst))

    @patch('pathlib.Path.write_text')
    def test_write_text_success(self, mock_write_text):
        """write_text should call Path.write_text."""
        path = Path("/test/file.txt")
        content = "test content"
        
        write_text(path, content, encoding="utf-8")
        
        mock_write_text.assert_called_once_with(content, encoding="utf-8")

    @patch('pathlib.Path.touch')
    def test_touch_success(self, mock_touch):
        """touch should call Path.touch."""
        path = Path("/test/file")
        
        touch(path)
        
        mock_touch.assert_called_once()

    @patch('pathlib.Path.unlink')
    def test_unlink_success(self, mock_unlink):
        """unlink should call Path.unlink."""
        path = Path("/test/file")
        
        unlink(path, missing_ok=True)
        
        mock_unlink.assert_called_once_with(missing_ok=True)

    @patch('pathlib.Path.rmdir')
    def test_rmdir_success(self, mock_rmdir):
        """rmdir should call Path.rmdir."""
        path = Path("/test/dir")
        
        rmdir(path)
        
        mock_rmdir.assert_called_once()

    @patch('pathlib.Path.rename')
    def test_rename_success(self, mock_rename):
        """rename should call Path.rename."""
        src = Path("/src/file")
        dst = Path("/dst/file")
        
        rename(src, dst)
        
        mock_rename.assert_called_once_with(dst)

    @patch('shutil.copyfileobj')
    def test_copy_fileobj_success(self, mock_copyfileobj):
        """copy_fileobj should call shutil.copyfileobj."""
        src_file = Mock()
        dst_file = Mock()
        
        copy_fileobj(src_file, dst_file)
        
        mock_copyfileobj.assert_called_once_with(src_file, dst_file)


class TestDiskFileOperation:
    """Test the disk_file_operation context manager."""

    @patch('pathlib.Path.open')
    @patch('pathlib.Path.is_file')
    @patch('pathlib.Path.exists')
    def test_context_manager_success(self, mock_exists, mock_is_file, mock_open):
        """Context manager should open and close file correctly."""
        mock_file = Mock()
        mock_exists.return_value = True
        mock_is_file.return_value = True
        mock_open.return_value = mock_file
        path = Path("/test/file.txt")
        
        with disk_file_operation(path, 'r') as f:
            assert f is mock_file
        
        mock_open.assert_called_once_with('r')

    @patch('pathlib.Path.open')
    @patch('pathlib.Path.is_file')
    @patch('pathlib.Path.exists')
    @patch('micall.monitor.disk_operations.wait_for_retry')
    def test_context_manager_retry_on_failure(self, mock_wait, mock_exists, mock_is_file, mock_open):
        """Context manager should retry on OSError."""
        mock_file = Mock()
        mock_exists.return_value = True
        mock_is_file.return_value = True
        mock_open.side_effect = [
            OSError("File not found"),
            mock_file
        ]
        path = Path("/test/file.txt")
        
        with disk_file_operation(path, 'r'):
            pass
        
        assert mock_open.call_count == 2
        mock_wait.assert_called_once_with(1, "open(r)", False)


class TestRemoveEmptyDirectory:
    """Test the remove_empty_directory function."""

    @patch('pathlib.Path.rmdir')
    def test_remove_empty_directory_success(self, mock_rmdir):
        """Should remove directory if it's empty."""
        path = Path("/test/empty_dir")
        
        remove_empty_directory(path)
        
        mock_rmdir.assert_called_once()

    @patch('pathlib.Path.rmdir')
    def test_remove_empty_directory_not_empty(self, mock_rmdir):
        """Should not raise error if directory is not empty."""
        import errno
        mock_rmdir.side_effect = OSError(errno.ENOTEMPTY, "Directory not empty")
        path = Path("/test/not_empty_dir")
        
        remove_empty_directory(path)  # Should not raise
        
        mock_rmdir.assert_called_once()

    @patch('pathlib.Path.rmdir')
    def test_remove_empty_directory_error_checking(self, mock_rmdir):
        """Should re-raise non-ENOTEMPTY errors."""
        import errno
        mock_rmdir.side_effect = OSError(errno.EACCES, "Permission denied")
        path = Path("/test/no_access_dir")
        
        with pytest.raises(OSError):
            remove_empty_directory(path)
        
        mock_rmdir.assert_called_once()

    @patch('pathlib.Path.rmdir')
    @patch('micall.monitor.disk_operations.wait_for_retry')
    def test_remove_empty_directory_retry_on_rmdir_failure(self, mock_wait, mock_rmdir):
        """Should retry rmdir on failure."""
        mock_rmdir.side_effect = [OSError("Directory busy"), None]
        path = Path("/test/busy_dir")
        
        remove_empty_directory(path)
        
        assert mock_rmdir.call_count == 2
        mock_wait.assert_called_once_with(1, "remove_empty_directory", False)


class TestIntegrationWithRealFileSystem:
    """Integration tests using real file system operations."""

    def test_mkdir_p_real_filesystem(self):
        """Test mkdir_p with real filesystem."""
        with tempfile.TemporaryDirectory() as tmpdir:
            test_path = Path(tmpdir) / "test" / "nested" / "directory"
            
            # Should succeed without errors
            mkdir_p(test_path, parents=True, exist_ok=True)
            
            assert test_path.exists()
            assert test_path.is_dir()

    def test_write_text_and_unlink_real_filesystem(self):
        """Test write_text and unlink with real filesystem."""
        with tempfile.TemporaryDirectory() as tmpdir:
            test_file = Path(tmpdir) / "test.txt"
            content = "Hello, World!"
            
            # Write file
            write_text(test_file, content)
            assert test_file.exists()
            assert test_file.read_text() == content
            
            # Remove file
            unlink(test_file)
            assert not test_file.exists()

    def test_touch_real_filesystem(self):
        """Test touch with real filesystem."""
        with tempfile.TemporaryDirectory() as tmpdir:
            test_file = Path(tmpdir) / "touched.txt"
            
            touch(test_file)
            
            assert test_file.exists()
            assert test_file.stat().st_size == 0

    def test_copy_file_real_filesystem(self):
        """Test copy_file with real filesystem."""
        with tempfile.TemporaryDirectory() as tmpdir:
            src_file = Path(tmpdir) / "source.txt"
            dst_file = Path(tmpdir) / "destination.txt"
            content = "Copy this content"
            
            src_file.write_text(content)
            copy_file(src_file, dst_file)
            
            assert dst_file.exists()
            assert dst_file.read_text() == content

    def test_move_real_filesystem(self):
        """Test move with real filesystem."""
        with tempfile.TemporaryDirectory() as tmpdir:
            src_file = Path(tmpdir) / "source.txt"
            dst_file = Path(tmpdir) / "moved.txt"
            content = "Move this content"
            
            src_file.write_text(content)
            move(src_file, dst_file)
            
            assert not src_file.exists()
            assert dst_file.exists()
            assert dst_file.read_text() == content

    def test_rename_real_filesystem(self):
        """Test rename with real filesystem."""
        with tempfile.TemporaryDirectory() as tmpdir:
            src_file = Path(tmpdir) / "original.txt"
            dst_file = Path(tmpdir) / "renamed.txt"
            content = "Rename this file"
            
            src_file.write_text(content)
            rename(src_file, dst_file)
            
            assert not src_file.exists()
            assert dst_file.exists()
            assert dst_file.read_text() == content

    def test_rmtree_real_filesystem(self):
        """Test rmtree with real filesystem."""
        with tempfile.TemporaryDirectory() as tmpdir:
            test_dir = Path(tmpdir) / "test_directory"
            test_dir.mkdir()
            (test_dir / "file1.txt").write_text("content1")
            (test_dir / "file2.txt").write_text("content2")
            
            rmtree(test_dir)
            
            assert not test_dir.exists()

    def test_disk_file_operation_real_filesystem(self):
        """Test disk_file_operation context manager with real filesystem."""
        with tempfile.TemporaryDirectory() as tmpdir:
            test_file = Path(tmpdir) / "context_test.txt"
            content = "Context manager test"
            
            # Write using context manager
            with disk_file_operation(test_file, 'w') as f:
                f.write(content)
            
            # Read using context manager
            with disk_file_operation(test_file, 'r') as f:
                read_content = f.read()
            
            assert read_content == content

    def test_remove_empty_directory_real_filesystem(self):
        """Test remove_empty_directory with real filesystem."""
        with tempfile.TemporaryDirectory() as tmpdir:
            empty_dir = Path(tmpdir) / "empty"
            non_empty_dir = Path(tmpdir) / "non_empty"
            
            empty_dir.mkdir()
            non_empty_dir.mkdir()
            (non_empty_dir / "file.txt").write_text("content")
            
            # Remove empty directory
            remove_empty_directory(empty_dir)
            assert not empty_dir.exists()
            
            # Don't remove non-empty directory
            remove_empty_directory(non_empty_dir)
            assert non_empty_dir.exists()


class TestErrorScenarios:
    """Test various error scenarios and edge cases."""

    @patch('pathlib.Path.mkdir')
    def test_mkdir_p_permission_error_chain(self, mock_mkdir):
        """Test multiple permission errors in sequence."""
        path = Path("/protected/directory")
        mock_mkdir.side_effect = [
            PermissionError("Access denied"),
            PermissionError("Still denied"),
            None  # Success on third try
        ]
        
        with patch('micall.monitor.disk_operations.wait_for_retry') as mock_wait:
            mkdir_p(path)
        
        assert mock_mkdir.call_count == 3
        assert mock_wait.call_count == 2

    @patch('shutil.rmtree')
    @patch('pathlib.Path.exists')
    def test_rmtree_file_in_use_error(self, mock_exists, mock_rmtree):
        """Test handling of file-in-use errors."""
        path = Path("/busy/directory")
        mock_exists.return_value = True
        mock_rmtree.side_effect = [
            OSError("The process cannot access the file because it is being used by another process"),
            None
        ]
        
        with patch('micall.monitor.disk_operations.wait_for_retry') as mock_wait:
            rmtree(path)
        
        assert mock_rmtree.call_count == 2
        mock_wait.assert_called_once_with(1, "rmtree", False)

    @patch('pathlib.Path.open')
    @patch('pathlib.Path.is_file')
    @patch('pathlib.Path.exists')
    def test_disk_file_operation_intermittent_errors(self, mock_exists, mock_is_file, mock_open):
        """Test handling of intermittent file access errors."""
        path = Path("/intermittent/file.txt")
        mock_exists.return_value = True
        mock_is_file.return_value = True
        mock_open.side_effect = [
            IOError("Sharing violation"),
            IOError("File locked"),
            Mock()  # Success
        ]
        
        with patch('micall.monitor.disk_operations.wait_for_retry') as mock_wait:
            with disk_file_operation(path, 'r'):
                pass
        
        assert mock_open.call_count == 3
        assert mock_wait.call_count == 2
