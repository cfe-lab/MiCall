"""
Dynamic scoping for work_dir parameter.

This module provides a way to manage the work_dir parameter throughout
the call chain using Python's contextvars for thread-safe dynamic scoping.

Example usage:
    with WorkDir.using(Path("/tmp/work")):
        # Any function called here can access the work_dir
        denovo(fastq1, fastq2, fasta)
"""

from pathlib import Path
from contextvars import ContextVar
from contextlib import contextmanager
from dataclasses import dataclass


@dataclass
class WorkDir:
    """Manages work_dir with dynamic scoping."""

    work_dir: Path

    @staticmethod
    def _get() -> "WorkDir":
        """
        Get the current WorkDir.

        Returns:
            The current work_dir holder.

        Raises:
            LookupError: If no work_dir has been set.
        """
        return _work_dir.get()

    @staticmethod
    def get() -> Path:
        """
        Get the current work_dir path.

        Returns:
            The current work directory path.

        Raises:
            LookupError: If no work_dir has been set.
        """
        return WorkDir._get().work_dir

    @staticmethod
    @contextmanager
    def using(path: Path):
        """
        Set the work_dir for the duration of this context.

        The work_dir is restored to its previous state when exiting.

        Args:
            path: The work directory path to use in this context.

        Example:
            with WorkDir.using(Path("/tmp/work")):
                # work_dir is now /tmp/work
                do_something()
        """
        ctx = WorkDir(work_dir=path)
        token = _work_dir.set(ctx)
        try:
            yield ctx
        finally:
            _work_dir.reset(token)


# Global context variable for work_dir
_work_dir: ContextVar[WorkDir] = ContextVar("work_dir")
