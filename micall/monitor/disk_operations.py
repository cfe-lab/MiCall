"""
Low-level disk operations with network-style retry logic.

This module provides basic disk operation wrappers that use the same retry patterns
as the Kive network operations. Treats disk as a network resource with exponential
backoff for transient failures.
"""

import logging
import shutil
from pathlib import Path
from datetime import datetime, timedelta
from time import sleep

logger = logging.getLogger(__name__)

# Use same retry constants as network operations
MINIMUM_RETRY_WAIT = timedelta(seconds=5)
MAXIMUM_RETRY_WAIT = timedelta(days=1)


def calculate_retry_wait(min_wait, max_wait, attempt_count):
    """Calculate exponential backoff wait time (same as kive_watcher)."""
    min_seconds = int(min_wait.total_seconds())
    seconds = min_seconds * (2 ** (attempt_count - 1))
    seconds = min(seconds, max_wait.total_seconds())
    return timedelta(seconds=seconds)


def calculate_cumulative_wait_time(attempt_count):
    """Calculate total time that will have elapsed after all previous attempts."""
    total_time = timedelta()
    for i in range(1, attempt_count):
        total_time += calculate_retry_wait(MINIMUM_RETRY_WAIT, MAXIMUM_RETRY_WAIT, i)
    return total_time


def wait_for_retry(attempt_count, operation_name, start_time):
    """Wait with exponential backoff, logging errors after 1 hour, info messages before."""
    delay = calculate_retry_wait(MINIMUM_RETRY_WAIT, MAXIMUM_RETRY_WAIT, attempt_count)

    # Determine if we should log based on elapsed time
    elapsed = datetime.now() - start_time
    should_log_error = elapsed >= timedelta(hours=1)

    if should_log_error:
        logger.error(
            "Disk operation %s failed, waiting %s before retrying.",
            operation_name,
            delay,
            exc_info=True,
        )
    else:
        logger.info(
            "Disk operation %s failed, waiting %s before retrying.",
            operation_name,
            delay,
            exc_info=True,
        )
    sleep(delay.total_seconds())


def disk_retry(operation_name="disk operation"):
    """Decorator that adds network-style retry logic to disk operations."""

    def decorator(func):
        def wrapper(*args, **kwargs):
            attempt_count = 0
            max_attempts = 15  # Same as network operations
            start_time = None

            while attempt_count < max_attempts:
                attempt_count += 1
                try:
                    return func(*args, **kwargs)
                except (OSError, IOError) as ex:
                    if attempt_count >= max_attempts:
                        logger.error(
                            f"Disk operation {operation_name} failed after {max_attempts} attempts: {ex}"
                        )
                        raise

                    # Record start time on first failure
                    if start_time is None:
                        start_time = datetime.now()

                    wait_for_retry(attempt_count, operation_name, start_time)
                except:
                    # Non-disk errors should not be retried
                    raise

            return func(*args, **kwargs)  # This shouldn't be reached

        return wrapper

    return decorator


# Basic low-level disk operation wrappers


@disk_retry("mkdir")
def mkdir_p(path: Path, mode=0o777, parents=True, exist_ok=False):
    """Create directory with parents, network-aware retry."""
    path.mkdir(mode=mode, parents=parents, exist_ok=exist_ok)


@disk_retry("rmtree")
def rmtree(path: Path, ignore_errors=False):
    """Remove directory tree, network-aware retry."""
    if path.exists():
        shutil.rmtree(path, ignore_errors=ignore_errors)


@disk_retry("move")
def move(src: Path, dst: Path):
    """Move file or directory, network-aware retry."""
    shutil.move(str(src), str(dst))


@disk_retry("copy")
def copy_file(src: Path, dst: Path):
    """Copy file, network-aware retry."""
    shutil.copy2(str(src), str(dst))


@disk_retry("write_text")
def write_text(path: Path, content: str, encoding="utf-8"):
    """Write text to file, network-aware retry."""
    path.write_text(content, encoding=encoding)


@disk_retry("touch")
def touch(path: Path):
    """Touch file (create if not exists), network-aware retry."""
    path.touch()


@disk_retry("unlink")
def unlink(path: Path, missing_ok=True):
    """Remove file, network-aware retry."""
    path.unlink(missing_ok=missing_ok)


@disk_retry("rmdir")
def rmdir(path: Path):
    """Remove empty directory, network-aware retry."""
    try:
        path.rmdir()
    except FileNotFoundError:
        pass


@disk_retry("rename")
def rename(src: Path, dst: Path):
    """Rename/move file or directory, network-aware retry."""
    src.rename(dst)


@disk_retry("remove_empty_directory")
def remove_empty_directory(path: Path):
    """Remove empty directory, network-aware retry."""
    import errno

    try:
        path.rmdir()
    except OSError as ex:
        if ex.errno == errno.ENOTEMPTY:
            # Directory not empty is expected, just ignore
            pass
        else:
            # All other errors (including permission errors) should be re-raised
            # and will be subject to retry logic
            raise


@disk_retry("copyfileobj")
def copy_fileobj(src_file, dst_file):
    """Copy file object contents, network-aware retry."""
    shutil.copyfileobj(src_file, dst_file)


# Contextual operation wrapper for file operations
class disk_file_operation:
    """Context manager for file operations with retry logic."""

    def __init__(self, path: Path, mode: str, operation_name=None):
        self.path = path
        self.mode = mode
        self.operation_name = operation_name or f"open({mode})"
        self.file_handle = None

    def __enter__(self):
        @disk_retry(self.operation_name)
        def open_file():
            return self.path.open(self.mode)

        self.file_handle = open_file()
        return self.file_handle

    def __exit__(self, exc_type, exc_val, exc_tb):
        if self.file_handle:
            self.file_handle.close()
