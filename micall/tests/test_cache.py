"""
Tests for the cache module that handles disk-based caching of function results.
"""

import json
from pathlib import Path
from unittest.mock import patch

import pytest

from micall.utils import cache


@pytest.fixture
def temp_cache_dir(tmp_path):
    """Create a temporary cache directory and set it as CACHE_FOLDER."""
    cache_dir = tmp_path / "cache"
    cache_dir.mkdir()

    # Patch the CACHE_FOLDER for the duration of the test
    with patch.object(cache, "CACHE_FOLDER", str(cache_dir)):
        yield cache_dir


@pytest.fixture
def sample_files(tmp_path):
    """Create sample input and output files for testing."""
    input_dir = tmp_path / "inputs"
    input_dir.mkdir()
    output_dir = tmp_path / "outputs"
    output_dir.mkdir()

    # Create sample input files
    input1 = input_dir / "input1.txt"
    input1.write_text("This is input file 1")

    input2 = input_dir / "input2.txt"
    input2.write_text("This is input file 2")

    # Create sample output files
    output1 = output_dir / "output1.txt"
    output1.write_text("This is output file 1")

    output2 = output_dir / "output2.txt"
    output2.write_text("This is output file 2")

    return {
        "input1": input1,
        "input2": input2,
        "output1": output1,
        "output2": output2,
    }


@pytest.fixture(autouse=True)
def clear_hashes():
    """Clear the HASHES cache before each test."""
    cache.HASHES.clear()
    yield
    cache.HASHES.clear()


class TestHashFile:
    """Tests for the hash_file function."""

    def test_hash_file_basic(self, tmp_path):
        """Test that hash_file computes correct MD5 hash."""
        test_file = tmp_path / "test.txt"
        test_file.write_text("Hello, World!")

        # Clear the hash cache first
        cache.HASHES.clear()

        file_hash = cache.hash_file(test_file)

        # Verify it's a valid MD5 hash (32 hex characters)
        assert len(file_hash) == 32
        assert all(c in "0123456789abcdef" for c in file_hash)

        # Verify it matches expected hash for this content
        import hashlib

        expected_hash = hashlib.md5(b"Hello, World!").hexdigest()
        assert file_hash == expected_hash

    def test_hash_file_caching(self, tmp_path):
        """Test that hash_file caches results."""
        test_file = tmp_path / "test.txt"
        test_file.write_text("Content")

        cache.HASHES.clear()

        # First call computes hash
        hash1 = cache.hash_file(test_file)
        assert test_file in cache.HASHES

        # Second call should return cached value
        hash2 = cache.hash_file(test_file)
        assert hash1 == hash2
        assert cache.HASHES[test_file] == hash1

    def test_hash_file_different_content(self, tmp_path):
        """Test that different content produces different hashes."""
        file1 = tmp_path / "file1.txt"
        file1.write_text("Content A")

        file2 = tmp_path / "file2.txt"
        file2.write_text("Content B")

        cache.HASHES.clear()

        hash1 = cache.hash_file(file1)
        hash2 = cache.hash_file(file2)

        assert hash1 != hash2

    def test_hash_file_large_file(self, tmp_path):
        """Test hashing a file larger than the chunk size."""
        test_file = tmp_path / "large.txt"
        # Write more than 8192 bytes (the chunk size)
        test_file.write_text("X" * 20000)

        cache.HASHES.clear()

        file_hash = cache.hash_file(test_file)
        assert len(file_hash) == 32


class TestMakeCacheKey:
    """Tests for the _make_cache_key function."""

    def test_make_cache_key_basic(self, sample_files):
        """Test basic cache key generation."""
        inputs = {"input1": sample_files["input1"]}
        key = cache._make_cache_key("test_procedure", inputs)

        # Should be a valid MD5 hash
        assert len(key) == 32
        assert all(c in "0123456789abcdef" for c in key)

    def test_make_cache_key_includes_procedure(self, sample_files):
        """Test that cache key changes with different procedures."""
        inputs = {"input1": sample_files["input1"]}

        key1 = cache._make_cache_key("procedure_a", inputs)
        key2 = cache._make_cache_key("procedure_b", inputs)

        assert key1 != key2

    def test_make_cache_key_includes_input_files(self, sample_files):
        """Test that cache key changes with different input files."""
        inputs1 = {"input1": sample_files["input1"]}
        inputs2 = {"input1": sample_files["input2"]}

        key1 = cache._make_cache_key("test_procedure", inputs1)
        key2 = cache._make_cache_key("test_procedure", inputs2)

        assert key1 != key2

    def test_make_cache_key_sorted_inputs(self, sample_files):
        """Test that input order doesn't affect cache key."""
        inputs1 = {"a": sample_files["input1"], "b": sample_files["input2"]}
        inputs2 = {"b": sample_files["input2"], "a": sample_files["input1"]}

        key1 = cache._make_cache_key("test_procedure", inputs1)
        key2 = cache._make_cache_key("test_procedure", inputs2)

        assert key1 == key2

    def test_make_cache_key_with_none(self, sample_files):
        """Test cache key generation with None values."""
        inputs = {"input1": sample_files["input1"], "input2": None}

        key = cache._make_cache_key("test_procedure", inputs)
        assert len(key) == 32

    def test_make_cache_key_deterministic(self, sample_files):
        """Test that cache key is deterministic."""
        inputs = {"input1": sample_files["input1"]}

        key1 = cache._make_cache_key("test_procedure", inputs)
        key2 = cache._make_cache_key("test_procedure", inputs)

        assert key1 == key2


class TestLoadSaveCacheData:
    """Tests for _load_cache_data and _save_cache_data functions."""

    def test_load_cache_data_no_folder(self):
        """Test loading when CACHE_FOLDER is not set."""
        with patch.object(cache, "CACHE_FOLDER", None):
            data = cache._load_cache_data()
            assert data == {}

    def test_load_cache_data_nonexistent_file(self, temp_cache_dir):
        """Test loading when data.json doesn't exist."""
        data = cache._load_cache_data()
        assert data == {}

    def test_save_and_load_cache_data(self, temp_cache_dir):
        """Test saving and loading cache data."""
        test_data = {"key1": "value1", "key2": {"nested": "data"}}

        cache._save_cache_data(test_data)

        # Verify file was created
        data_file = temp_cache_dir / "data.json"
        assert data_file.exists()

        # Load and verify
        loaded_data = cache._load_cache_data()
        assert loaded_data == test_data

    def test_save_cache_data_creates_folder(self, tmp_path):
        """Test that save_cache_data creates the cache folder if needed."""
        cache_dir = tmp_path / "new_cache"

        with patch.object(cache, "CACHE_FOLDER", str(cache_dir)):
            test_data = {"key": "value"}
            cache._save_cache_data(test_data)

            assert cache_dir.exists()
            assert (cache_dir / "data.json").exists()

    def test_save_cache_data_no_folder(self):
        """Test saving when CACHE_FOLDER is not set."""
        with patch.object(cache, "CACHE_FOLDER", None):
            # Should not raise an error
            cache._save_cache_data({"key": "value"})

    def test_load_cache_data_invalid_json(self, temp_cache_dir):
        """Test loading when data.json contains invalid JSON."""
        data_file = temp_cache_dir / "data.json"
        data_file.write_text("invalid json {")

        # Should raise JSONDecodeError since there's no error handling
        with pytest.raises(json.JSONDecodeError):
            cache._load_cache_data()


class TestCacheGet:
    """Tests for the get function."""

    def test_get_no_cache_folder(self, sample_files):
        """Test get when CACHE_FOLDER is not set."""
        with patch.object(cache, "CACHE_FOLDER", None):
            inputs = {"input1": sample_files["input1"]}
            result = cache.get("test_procedure", inputs)
            assert result is None

    def test_get_cache_miss(self, temp_cache_dir, sample_files):
        """Test get when cache entry doesn't exist."""
        inputs = {"input1": sample_files["input1"]}
        result = cache.get("test_procedure", inputs)
        assert result is None

    def test_get_single_file(self, temp_cache_dir, sample_files):
        """Test retrieving a single cached file."""
        # Set up cache
        inputs = {"input1": sample_files["input1"]}
        output = sample_files["output1"]

        cache.set("test_procedure", inputs, output)

        # Retrieve from cache
        result = cache.get("test_procedure", inputs)

        assert result is not None
        assert isinstance(result, Path)
        assert result.exists()
        assert result.read_text() == output.read_text()

    def test_get_multiple_files(self, temp_cache_dir, sample_files):
        """Test retrieving multiple cached files."""
        # Set up cache
        inputs = {"input1": sample_files["input1"]}
        outputs = {
            "output1": sample_files["output1"],
            "output2": sample_files["output2"],
        }

        cache.set("test_procedure", inputs, outputs)

        # Retrieve from cache
        result = cache.get("test_procedure", inputs)

        assert result is not None
        assert isinstance(result, dict)
        assert "output1" in result
        assert "output2" in result
        assert result["output1"].read_text() == sample_files["output1"].read_text()
        assert result["output2"].read_text() == sample_files["output2"].read_text()

    def test_get_with_none_value(self, temp_cache_dir, sample_files):
        """Test retrieving cached results with None values."""
        inputs = {"input1": sample_files["input1"]}
        outputs = {"output1": sample_files["output1"], "output2": None}

        cache.set("test_procedure", inputs, outputs)

        result = cache.get("test_procedure", inputs)

        assert result is not None
        assert isinstance(result, dict)
        assert result["output1"] is not None
        assert result["output2"] is None

    def test_get_missing_cached_file(self, temp_cache_dir, sample_files):
        """Test get when cached file has been deleted."""
        # Set up cache
        inputs = {"input1": sample_files["input1"]}
        output = sample_files["output1"]

        cache.set("test_procedure", inputs, output)

        # Delete the cached file
        cache_key = cache._make_cache_key("test_procedure", inputs)
        data = cache._load_cache_data()
        cached_file_hash = data[cache_key]
        cached_file = temp_cache_dir / cached_file_hash
        cached_file.unlink()

        # Should return None since file is missing
        result = cache.get("test_procedure", inputs)
        assert result is None


class TestCacheSet:
    """Tests for the set function."""

    def test_set_no_cache_folder(self, sample_files):
        """Test set when CACHE_FOLDER is not set."""
        with patch.object(cache, "CACHE_FOLDER", None):
            inputs = {"input1": sample_files["input1"]}
            output = sample_files["output1"]

            # Should not raise an error
            cache.set("test_procedure", inputs, output)

    def test_set_single_file(self, temp_cache_dir, sample_files):
        """Test caching a single file."""
        inputs = {"input1": sample_files["input1"]}
        output = sample_files["output1"]

        cache.set("test_procedure", inputs, output)

        # Verify data.json was created
        data_file = temp_cache_dir / "data.json"
        assert data_file.exists()

        # Verify cache entry
        data = json.loads(data_file.read_text())
        cache_key = cache._make_cache_key("test_procedure", inputs)
        assert cache_key in data

        # Verify cached file exists
        file_hash = data[cache_key]
        cached_file = temp_cache_dir / file_hash
        assert cached_file.exists()
        assert cached_file.read_text() == output.read_text()

    def test_set_multiple_files(self, temp_cache_dir, sample_files):
        """Test caching multiple files."""
        inputs = {"input1": sample_files["input1"]}
        outputs = {
            "output1": sample_files["output1"],
            "output2": sample_files["output2"],
        }

        cache.set("test_procedure", inputs, outputs)

        # Verify cache entry
        data = cache._load_cache_data()
        cache_key = cache._make_cache_key("test_procedure", inputs)
        assert cache_key in data

        cache_entry = data[cache_key]
        assert isinstance(cache_entry, dict)
        assert "output1" in cache_entry
        assert "output2" in cache_entry

        # Verify cached files exist
        for key, file_hash in cache_entry.items():
            cached_file = temp_cache_dir / file_hash
            assert cached_file.exists()

    def test_set_with_none_value(self, temp_cache_dir, sample_files):
        """Test caching with None values."""
        inputs = {"input1": sample_files["input1"]}
        outputs = {"output1": sample_files["output1"], "output2": None}

        cache.set("test_procedure", inputs, outputs)

        data = cache._load_cache_data()
        cache_key = cache._make_cache_key("test_procedure", inputs)
        cache_entry = data[cache_key]

        assert cache_entry["output1"] is not None
        assert cache_entry["output2"] is None

    def test_set_deduplicates_files(self, temp_cache_dir, sample_files):
        """Test that files with same content are deduplicated."""
        # Create two files with identical content
        file1 = sample_files["output1"]
        file2 = sample_files["output1"].parent / "duplicate.txt"
        file2.write_text(file1.read_text())

        inputs1 = {"input": sample_files["input1"]}
        inputs2 = {"input": sample_files["input2"]}

        cache.set("test_procedure", inputs1, file1)
        cache.set("test_procedure", inputs2, file2)

        # Both should point to the same cached file
        hash1 = cache.hash_file(file1)
        hash2 = cache.hash_file(file2)
        assert hash1 == hash2

        # Only one cached file should exist
        cached_file = temp_cache_dir / hash1
        assert cached_file.exists()

    def test_set_invalid_output_type(self, temp_cache_dir, sample_files):
        """Test set with invalid output type raises ValueError."""
        inputs = {"input1": sample_files["input1"]}

        with pytest.raises(ValueError, match="Unsupported output type"):
            cache.set("test_procedure", inputs, "invalid_output")

    def test_set_creates_cache_folder(self, tmp_path, sample_files):
        """Test that set creates cache folder if it doesn't exist."""
        cache_dir = tmp_path / "new_cache"

        with patch.object(cache, "CACHE_FOLDER", str(cache_dir)):
            inputs = {"input1": sample_files["input1"]}
            output = sample_files["output1"]

            cache.set("test_procedure", inputs, output)

            assert cache_dir.exists()


class TestCachedDecorator:
    """Tests for the @cached decorator."""

    def test_cached_basic(self, temp_cache_dir, sample_files):
        """Test basic decorator functionality."""
        call_count = 0

        @cache.cached("test_proc")
        def test_func(input1: Path) -> Path:
            nonlocal call_count
            call_count += 1
            return sample_files["output1"]

        # First call executes function
        result1 = test_func(input1=sample_files["input1"])
        assert call_count == 1
        assert result1.read_text() == sample_files["output1"].read_text()

        # Second call with same input should use cache
        result2 = test_func(input1=sample_files["input1"])
        assert call_count == 1  # Not called again
        assert result2.read_text() == sample_files["output1"].read_text()

    def test_cached_different_inputs(self, temp_cache_dir, sample_files):
        """Test that different inputs bypass cache."""
        call_count = 0

        @cache.cached("test_proc")
        def test_func(input1: Path) -> Path:
            nonlocal call_count
            call_count += 1
            return sample_files["output1"]

        # First call
        test_func(input1=sample_files["input1"])
        assert call_count == 1

        # Second call with different input
        test_func(input1=sample_files["input2"])
        assert call_count == 2  # Called again

    def test_cached_multiple_inputs(self, temp_cache_dir, sample_files):
        """Test decorator with multiple input parameters."""
        call_count = 0

        @cache.cached("test_proc")
        def test_func(input1: Path, input2: Path) -> Path:
            nonlocal call_count
            call_count += 1
            return sample_files["output1"]

        test_func(input1=sample_files["input1"], input2=sample_files["input2"])
        assert call_count == 1

        test_func(input1=sample_files["input1"], input2=sample_files["input2"])
        assert call_count == 1  # Cached

    def test_cached_with_none(self, temp_cache_dir, sample_files):
        """Test decorator with None input values."""
        call_count = 0

        @cache.cached("test_proc")
        def test_func(input1: Path, input2: Path | None) -> Path:
            nonlocal call_count
            call_count += 1
            return sample_files["output1"]

        test_func(input1=sample_files["input1"], input2=None)
        assert call_count == 1

        test_func(input1=sample_files["input1"], input2=None)
        assert call_count == 1  # Cached

    def test_cached_invalid_input_type(self, temp_cache_dir, sample_files):
        """Test decorator raises error for invalid input types."""

        @cache.cached("test_proc")
        def test_func(input1: Path) -> Path:
            return sample_files["output1"]

        with pytest.raises(ValueError, match="must be of type Path or None"):
            test_func(input1="not_a_path")

    def test_cached_no_cache_folder(self, sample_files):
        """Test decorator when CACHE_FOLDER is not set."""
        with patch.object(cache, "CACHE_FOLDER", None):
            call_count = 0

            @cache.cached("test_proc")
            def test_func(input1: Path) -> Path:
                nonlocal call_count
                call_count += 1
                return sample_files["output1"]

            # Should execute every time without caching
            test_func(input1=sample_files["input1"])
            assert call_count == 1

            test_func(input1=sample_files["input1"])
            assert call_count == 2  # Called again, no caching

    def test_cached_multiple_outputs(self, temp_cache_dir, sample_files):
        """Test decorator with dict return value."""
        call_count = 0

        @cache.cached("test_proc")
        def test_func(input1: Path) -> dict[str, Path]:
            nonlocal call_count
            call_count += 1
            return {
                "output1": sample_files["output1"],
                "output2": sample_files["output2"],
            }

        result1 = test_func(input1=sample_files["input1"])
        assert call_count == 1
        assert "output1" in result1
        assert "output2" in result1

        result2 = test_func(input1=sample_files["input1"])
        assert call_count == 1  # Cached
        assert result2["output1"].read_text() == sample_files["output1"].read_text()


class TestIntegration:
    """Integration tests for the cache module."""

    def test_full_workflow(self, temp_cache_dir, sample_files):
        """Test complete cache workflow: set, get, use."""
        procedure = "denovo[iva]"
        inputs = {"fastq1": sample_files["input1"], "fastq2": sample_files["input2"]}
        output = sample_files["output1"]

        # Initially no cache
        result = cache.get(procedure, inputs)
        assert result is None

        # Cache the result
        cache.set(procedure, inputs, output)

        # Now it should be cached
        result = cache.get(procedure, inputs)
        assert result is not None
        assert result.read_text() == output.read_text()

        # Verify cache structure
        data_file = temp_cache_dir / "data.json"
        assert data_file.exists()

        data = json.loads(data_file.read_text())
        assert len(data) == 1

        # Verify cached file exists with correct name
        cache_key = cache._make_cache_key(procedure, inputs)
        file_hash = data[cache_key]
        cached_file = temp_cache_dir / file_hash
        assert cached_file.exists()

    def test_multiple_procedures(self, temp_cache_dir, sample_files):
        """Test caching results from different procedures."""
        inputs = {"input1": sample_files["input1"]}

        cache.set("procedure_a", inputs, sample_files["output1"])
        cache.set("procedure_b", inputs, sample_files["output2"])

        # Both should be cached separately
        result_a = cache.get("procedure_a", inputs)
        result_b = cache.get("procedure_b", inputs)

        assert result_a is not None
        assert result_b is not None
        assert result_a.read_text() == sample_files["output1"].read_text()
        assert result_b.read_text() == sample_files["output2"].read_text()

    def test_cache_persistence(self, temp_cache_dir, sample_files):
        """Test that cache persists across module reloads."""
        inputs = {"input1": sample_files["input1"]}
        output = sample_files["output1"]

        # Cache something
        cache.set("test_procedure", inputs, output)

        # Simulate module reload by clearing in-memory cache
        cache.HASHES.clear()

        # Should still be able to retrieve from disk
        result = cache.get("test_procedure", inputs)
        assert result is not None
        assert result.read_text() == output.read_text()
