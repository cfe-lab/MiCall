"""Tests for the RemapCallback dynamic scoping module."""

from micall.utils.remap_callback import RemapCallback


def test_remap_callback_get_without_context():
    """Test that RemapCallback.get() returns None when no context is set."""
    result = RemapCallback.get()
    assert result is None


def test_remap_callback_get_with_context():
    """Test that RemapCallback.get() returns the callback function when context is set."""
    def my_callback(message, progress, max_progress):
        pass

    with RemapCallback.using(my_callback):
        result = RemapCallback.get()
        assert result is my_callback


def test_remap_callback_nested_contexts():
    """Test that nested contexts work correctly."""
    def callback1(message, progress, max_progress):
        return 1

    def callback2(message, progress, max_progress):
        return 2

    with RemapCallback.using(callback1):
        assert RemapCallback.get() is callback1

        with RemapCallback.using(callback2):
            assert RemapCallback.get() is callback2

        # After exiting inner context, should return to outer callback
        assert RemapCallback.get() is callback1

    # After exiting all contexts, should return to None
    assert RemapCallback.get() is None


def test_remap_callback_none_value():
    """Test that None can be explicitly set as the callback."""
    with RemapCallback.using(None):
        result = RemapCallback.get()
        assert result is None


def test_remap_callback_callable():
    """Test that the callback can be called if it's not None."""
    results = []

    def my_callback(message, progress, max_progress):
        results.append((message, progress, max_progress))

    with RemapCallback.using(my_callback):
        callback = RemapCallback.get()
        if callback:
            callback("test", 5, 10)

    assert results == [("test", 5, 10)]
