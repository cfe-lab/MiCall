"""
File operations with exponential backoff retry logic.

This module provides a consistent interface for file operations that can fail
due to network drives, permissions, or resource contention. All operations
use the same exponential backoff pattern as the network operations.
"""

import errno
import functools
import logging
import os
import shutil
import tarfile
import time
from datetime import timedelta
from itertools import count
from pathlib import Path
from typing import Callable, Any, Optional, Union

logger = logging.getLogger(__name__)

# Use the same constants as network operations for consistency
MINIMUM_RETRY_WAIT = timedelta(seconds=1)  # Shorter for disk ops
MAXIMUM_RETRY_WAIT = timedelta(seconds=120)  # Much shorter than network (2 minutes)
MAX_RETRIES = 5


class FileOperationError(Exception):
    """Base exception for file operation failures."""
    pass


class DiskSpaceError(FileOperationError):
    """Raised when disk space is insufficient."""
    pass


class PermissionError(FileOperationError):
    """Raised when file permissions prevent operation."""
    pass


class CorruptedFileError(FileOperationError):
    """Raised when file is corrupted or malformed."""
    pass


def calculate_retry_wait(min_wait: timedelta, max_wait: timedelta, attempt_count: int) -> timedelta:
    """Calculate exponential backoff delay."""
    min_seconds = int(min_wait.total_seconds())
    seconds = min_seconds * (2 ** (attempt_count - 1))
    seconds = min(seconds, max_wait.total_seconds())
    return timedelta(seconds=seconds)


def wait_for_retry(attempt_count: int, is_logged: bool = True, operation_name: str = "file operation"):
    """Wait with exponential backoff before retrying."""
    delay = calculate_retry_wait(MINIMUM_RETRY_WAIT, MAXIMUM_RETRY_WAIT, attempt_count)
    if is_logged:
        logger.warning('Retrying %s in %s (attempt %d)', operation_name, delay, attempt_count)
    time.sleep(delay.total_seconds())


def with_file_retry(max_retries: int = MAX_RETRIES, 
                    operation_name: Optional[str] = None):
    """
    Decorator to add exponential backoff retry logic to file operations.
    
    Args:
        max_retries: Maximum number of retry attempts
        operation_name: Human-readable name for logging (auto-detected if None)
    """
    def decorator(func: Callable) -> Callable:
        @functools.wraps(func)
        def wrapper(*args, **kwargs) -> Any:
            op_name = operation_name or f"{func.__name__}"

            for attempt_count in count(1):
                try:
                    return func(*args, **kwargs)
                except OSError as e:
                    # Handle specific OS errors with appropriate exceptions
                    if e.errno == errno.ENOSPC:
                        raise DiskSpaceError(f"Disk full during {op_name}: {e}")
                    elif e.errno == errno.EACCES:
                        raise PermissionError(f"Permission denied for {op_name}: {e}")
                    elif e.errno == errno.ENOENT:
                        # File not found - don't retry, it's probably not transient
                        raise FileNotFoundError(f"File not found during {op_name}: {e}")
                    elif attempt_count >= max_retries:
                        raise FileOperationError(f"Failed {op_name} after {max_retries} attempts: {e}")
                    else:
                        # Transient error - retry
                        logger.warning(f"OS error during {op_name} (attempt {attempt_count}): {e}")
                        wait_for_retry(attempt_count, operation_name=op_name)
                except (IOError, FileNotFoundError) as e:
                    if attempt_count >= max_retries:
                        raise FileOperationError(f"Failed {op_name} after {max_retries} attempts: {e}")
                    else:
                        logger.warning(f"IO error during {op_name} (attempt {attempt_count}): {e}")
                        wait_for_retry(attempt_count, operation_name=op_name)
                except Exception as e:
                    # Unexpected error - don't retry unknown exceptions
                    logger.error(f"Unexpected error during {op_name}: {e}")
                    raise
        
        return wrapper
    return decorator


# File operation wrappers with retry logic

@with_file_retry(operation_name="create directory")
def safe_mkdir(path: Union[Path, str], parents: bool = False, exist_ok: bool = True) -> None:
    """Create directory with retry logic."""
    Path(path).mkdir(parents=parents, exist_ok=exist_ok)


@with_file_retry(operation_name="remove directory tree")
def safe_rmtree(path: Union[Path, str], ignore_errors: bool = False) -> None:
    """Remove directory tree with retry logic."""
    shutil.rmtree(path, ignore_errors=ignore_errors)


@with_file_retry(operation_name="remove file")
def safe_unlink(path: Union[Path, str]) -> None:
    """Remove file with retry logic."""
    Path(path).unlink()


@with_file_retry(operation_name="move/rename file")
def safe_rename(src: Union[Path, str], dst: Union[Path, str]) -> None:
    """Move/rename file with retry logic."""
    os.rename(str(src), str(dst))


@with_file_retry(operation_name="copy file")
def safe_copyfile(src: Union[Path, str], dst: Union[Path, str]) -> None:
    """Copy file with retry logic."""
    shutil.copyfile(src, dst)


@with_file_retry(operation_name="copy file object")
def safe_copyfileobj(fsrc, fdst, length: int = 16*1024) -> None:
    """Copy file object with retry logic."""
    shutil.copyfileobj(fsrc, fdst, length)


@with_file_retry(operation_name="open file")
def safe_open(path: Union[Path, str], mode: str = 'r', **kwargs):
    """Open file with retry logic."""
    return open(path, mode, **kwargs)


@with_file_retry(operation_name="write text file")
def safe_write_text(path: Union[Path, str], content: str, **kwargs) -> None:
    """Write text to file with retry logic."""
    Path(path).write_text(content, **kwargs)


@with_file_retry(operation_name="extract tar file")
def safe_tar_extract(tar_path: Union[Path, str], extract_path: Union[Path, str]) -> None:
    """Extract tar file with retry logic and validation."""
    try:
        with tarfile.open(tar_path, 'r') as tar:
            # Validate tar file integrity and security
            for member in tar.getmembers():
                if not member.isfile():
                    continue
                # Check for path traversal attacks
                if os.path.isabs(member.name) or ".." in member.name:
                    raise CorruptedFileError(f"Unsafe path in tar: {member.name}")
            tar.extractall(extract_path)
    except tarfile.ReadError as e:
        raise CorruptedFileError(f"Corrupted tar file {tar_path}: {e}")


@with_file_retry(operation_name="open tar file")
def safe_tar_open(path: Union[Path, str], mode: str = 'r'):
    """Open tar file with retry logic and corruption detection."""
    try:
        # Use positional arguments to avoid type checker issues
        return tarfile.open(str(path))
    except tarfile.ReadError as e:
        raise CorruptedFileError(f"Cannot open tar file {path}: {e}")


class SafeFileContext:
    """Context manager for safe file operations with automatic cleanup."""
    
    def __init__(self, path: Union[Path, str], mode: str = 'r', **kwargs):
        self.path = Path(path)
        self.mode = mode
        self.kwargs = kwargs
        self.file = None
        
    def __enter__(self):
        self.file = safe_open(self.path, self.mode, **self.kwargs)
        return self.file
        
    def __exit__(self, exc_type, exc_val, exc_tb):
        if self.file:
            try:
                self.file.close()
            except Exception as e:
                logger.warning(f"Error closing file {self.path}: {e}")


# Convenience functions for common patterns

def safe_file_operation(operation: Callable, *args, max_retries: int = MAX_RETRIES, 
                       operation_name: str = "file operation", **kwargs) -> Any:
    """
    Generic wrapper for any file operation with retry logic.
    
    Args:
        operation: Function to execute
        *args: Arguments to pass to the function
        max_retries: Maximum retry attempts
        operation_name: Name for logging
        **kwargs: Keyword arguments to pass to the function
    """
    @with_file_retry(max_retries=max_retries, operation_name=operation_name)
    def wrapper():
        return operation(*args, **kwargs)
    
    return wrapper()
