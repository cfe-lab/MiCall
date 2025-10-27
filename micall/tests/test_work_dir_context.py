"""
This test shows how WorkDir provides dynamic scoping, allowing functions
deep in the call chain to access the work_dir without it being explicitly
passed through all intermediate functions.
"""

from pathlib import Path
from micall.utils.work_dir import WorkDir


def test_work_dir_basic():
    """Test basic get/set functionality of WorkDir."""
    test_path = Path("/tmp/test_work")

    with WorkDir.using(test_path):
        # Inside the scope, we can retrieve the work_dir
        assert WorkDir.get() == test_path

    # Outside the scope, get_path raises LookupError
    try:
        WorkDir.get()
        assert False, "Should have raised LookupError"
    except LookupError:
        pass  # Expected


def test_work_dir_nested():
    """Test nested scopes - inner scope should override outer."""
    outer_path = Path("/tmp/outer")
    inner_path = Path("/tmp/inner")

    with WorkDir.using(outer_path):
        assert WorkDir.get() == outer_path

        with WorkDir.using(inner_path):
            # Inner scope overrides
            assert WorkDir.get() == inner_path

        # Back to outer scope after exiting inner
        assert WorkDir.get() == outer_path


def test_work_dir_dynamic_scoping():
    """Test dynamic scoping - deeply nested function can access work_dir."""

    def level3_function():
        """Deepest function that needs work_dir but doesn't take it as parameter."""
        return WorkDir.get()

    def level2_function():
        """Middle function that doesn't care about work_dir."""
        return level3_function()

    def level1_function():
        """Top function that sets the scope."""
        test_path = Path("/tmp/dynamic_scope_test")
        with WorkDir.using(test_path):
            result = level2_function()
        return result

    # level3_function gets the path without level2_function passing it
    result = level1_function()
    assert result == Path("/tmp/dynamic_scope_test")
