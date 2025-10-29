"""
Dynamic scoping for remap callback parameter.

This module provides a way to manage the remap callback parameter throughout
the call chain using Python's contextvars for thread-safe dynamic scoping.

Example usage:
    def my_callback(message, progress, max_progress):
        print(f"{message}: {progress}/{max_progress}")

    with RemapCallback.using(my_callback):
        # Any function called here can access the callback
        remap(fastq1, fastq2, ...)
"""

from typing import Optional, Callable
from contextvars import ContextVar
from contextlib import contextmanager
from dataclasses import dataclass


# Type for remap callback function
RemapCallbackFunction = Optional[Callable[[str, int, int], None]]


_remap_callback: ContextVar["RemapCallback"] = ContextVar("remap_callback")


@dataclass
class RemapCallback:
    """Manages remap callback with dynamic scoping."""

    callback: RemapCallbackFunction

    @staticmethod
    def _get() -> "RemapCallback":
        """
        Get the current RemapCallback.

        Returns:
            The current callback holder.

        Raises:
            LookupError: If no callback has been set.
        """
        return _remap_callback.get()

    @staticmethod
    def get() -> RemapCallbackFunction:
        """
        Get the current remap callback function.

        Returns:
            The current callback function, or None if not set.
        """
        try:
            callback_obj = RemapCallback._get().callback
            return callback_obj
        except LookupError:
            return None

    @staticmethod
    @contextmanager
    def using(callback: RemapCallbackFunction):
        """
        Context manager to set the callback for the current context.

        Args:
            callback: The callback function to use, or None.

        Example:
            def my_callback(msg, progress, max_progress):
                print(f"{msg}: {progress}/{max_progress}")

            with RemapCallback.using(my_callback):
                # Code here can access the callback via RemapCallback.get()
                remap(...)
        """
        token = _remap_callback.set(RemapCallback(callback))
        try:
            yield
        finally:
            _remap_callback.reset(token)
