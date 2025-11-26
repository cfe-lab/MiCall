"""
Dynamic scoping for stderr parameter.

This module provides a way to manage the stderr parameter throughout
the call chain using Python's contextvars for thread-safe dynamic scoping.

Example usage:
    with open("/tmp/stderr.log", "w") as stderr_file:
        with Stderr.using(stderr_file):
            # Any function called here can access the stderr file object
            prelim_map(fastq1, fastq2, output)
"""

from typing import TextIO, Optional
from contextvars import ContextVar
from contextlib import contextmanager
from dataclasses import dataclass
import sys


@dataclass
class Stderr:
    """Manages stderr with dynamic scoping."""

    stderr: Optional[TextIO]

    @staticmethod
    def _get() -> "Stderr":
        """
        Get the current Stderr.

        Returns:
            The current stderr holder.

        Raises:
            LookupError: If no stderr has been set.
        """
        return _stderr.get()

    @staticmethod
    def get() -> TextIO:
        """
        Get the current stderr file object.

        Returns:
            The current stderr file object, or sys.stderr if not set.
        """
        try:
            stderr_obj = Stderr._get().stderr
            return stderr_obj if stderr_obj is not None else sys.stderr
        except LookupError:
            return sys.stderr

    @staticmethod
    @contextmanager
    def using(stderr: Optional[TextIO]):
        """
        Set the stderr for the duration of this context.

        The stderr is restored to its previous state when exiting.

        Args:
            stderr: The stderr file object to use in this context, or None for sys.stderr.

        Example:
            with open("/tmp/stderr.log", "w") as f:
                with Stderr.using(f):
                    # stderr is now the file object
                    do_something()
        """
        ctx = Stderr(stderr=stderr)
        token = _stderr.set(ctx)
        try:
            yield ctx
        finally:
            _stderr.reset(token)


# Global context variable for stderr
_stderr: ContextVar[Stderr] = ContextVar("stderr")
