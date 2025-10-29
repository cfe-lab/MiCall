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


class TestMakeInputKey:
    """Tests for the _make_input_key function."""

    def test_make_input_key_basic(self, sample_files):
        """Test basic input key generation."""
        inputs = {"input1": sample_files["input1"]}
        key = cache._make_cache_key(inputs)

        # Should be a dict with input names and hashes
        assert isinstance(key, dict)
        assert "input1" in key
        assert isinstance(key["input1"], str)
        assert len(key["input1"]) == 32  # MD5 hash

    def test_make_input_key_multiple_inputs(self, sample_files):
        """Test input key with multiple inputs."""
        inputs = {"input1": sample_files["input1"], "input2": sample_files["input2"]}
        key = cache._make_cache_key(inputs)

        assert len(key) == 2
        assert "input1" in key
        assert "input2" in key
        assert key["input1"] != key["input2"]

    def test_make_input_key_different_files(self, sample_files):
        """Test that different files produce different hashes."""
        inputs1 = {"input1": sample_files["input1"]}
        inputs2 = {"input1": sample_files["input2"]}

        key1 = cache._make_cache_key(inputs1)
        key2 = cache._make_cache_key(inputs2)

        assert key1["input1"] != key2["input1"]

    def test_make_input_key_sorted(self, sample_files):
        """Test that input keys are sorted consistently."""
        inputs1 = {"a": sample_files["input1"], "b": sample_files["input2"]}
        inputs2 = {"b": sample_files["input2"], "a": sample_files["input1"]}

        key1 = cache._make_cache_key(inputs1)
        key2 = cache._make_cache_key(inputs2)

        assert key1 == key2

    def test_make_input_key_with_none(self, sample_files):
        """Test input key generation with None values."""
        inputs = {"input1": sample_files["input1"], "input2": None}

        key = cache._make_cache_key(inputs)
        assert key["input1"] is not None
        assert key["input2"] is None

    def test_make_input_key_deterministic(self, sample_files):
        """Test that input key is deterministic."""
        inputs = {"input1": sample_files["input1"]}

        key1 = cache._make_cache_key(inputs)
        key2 = cache._make_cache_key(inputs)

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

        # Get the stored data structure
        data = cache._load_cache_data()
        # Find the entry for this procedure
        entries = data.get("test_procedure", [])
        assert len(entries) > 0

        # Delete the cached file
        file_hash = entries[0]["outputs"]
        cached_file = temp_cache_dir / file_hash
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

        # Verify cache entry structure
        data = json.loads(data_file.read_text())
        assert "test_procedure" in data
        assert isinstance(data["test_procedure"], list)
        assert len(data["test_procedure"]) == 1

        entry = data["test_procedure"][0]
        assert "inputs" in entry
        assert "outputs" in entry

        # Verify cached file exists
        file_hash = entry["outputs"]
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
        assert "test_procedure" in data
        entries = data["test_procedure"]
        assert len(entries) == 1

        entry = entries[0]
        assert "inputs" in entry
        assert "outputs" in entry

        cache_entry = entry["outputs"]
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
        entries = data["test_procedure"]
        cache_entry = entries[0]["outputs"]

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

        with pytest.raises(ValueError, match="must be either a Path, None, listed in parameters, or listed in outputs"):
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

    def test_cached_positional_arguments(self, temp_cache_dir, sample_files):
        """Test decorator with positional arguments (not just kwargs)."""
        call_count = 0

        @cache.cached("test_proc")
        def test_func(input1: Path, input2: Path, input3: Path | None = None) -> Path:
            nonlocal call_count
            call_count += 1
            return sample_files["output1"]

        # Call with positional arguments (like denovo does)
        test_func(sample_files["input1"], sample_files["input2"], None)
        assert call_count == 1

        # Second call with same positional arguments should use cache
        test_func(sample_files["input1"], sample_files["input2"], None)
        assert call_count == 1  # Cached

        # Verify the inputs were captured correctly in the cache
        data = cache._load_cache_data()
        assert "test_proc" in data
        entries = data["test_proc"]
        assert len(entries) == 1

        # Check that all three parameters are in the inputs
        inputs = entries[0]["inputs"]
        assert "input1" in inputs
        assert "input2" in inputs
        assert "input3" in inputs
        assert inputs["input1"] is not None
        assert inputs["input2"] is not None
        assert inputs["input3"] is None


class TestCachedDecoratorWithParameters:
    """Tests for the @cached decorator with parameters argument."""

    def test_parameters_basic(self, temp_cache_dir, sample_files, tmp_path):
        """Test caching with simple parameters."""
        call_count = 0
        output_file = tmp_path / "result.txt"

        @cache.cached("test_proc", parameters=["threshold"], outputs=["output"])
        def process(input_file: Path, output: Path, threshold: int) -> None:
            nonlocal call_count
            call_count += 1
            output.write_text(f"Processed {input_file.name} with threshold {threshold}")

        # First call
        process(sample_files["input1"], output_file, 10)
        assert call_count == 1
        assert output_file.exists()
        content1 = output_file.read_text()

        # Delete output to verify it gets restored from cache
        output_file.unlink()

        # Second call with same parameters - should use cache
        process(sample_files["input1"], output_file, 10)
        assert call_count == 1  # Not called again
        assert output_file.exists()
        assert output_file.read_text() == content1

        # Call with different threshold - should execute again
        process(sample_files["input1"], output_file, 20)
        assert call_count == 2
        assert "threshold 20" in output_file.read_text()

    def test_parameters_multiple(self, temp_cache_dir, sample_files, tmp_path):
        """Test caching with multiple parameters."""
        call_count = 0
        output_file = tmp_path / "result.txt"

        @cache.cached("test_proc",
                     parameters=["method", "threshold", "enabled"],
                     outputs=["output"])
        def process(input_file: Path, output: Path, method: str,
                   threshold: int, enabled: bool) -> None:
            nonlocal call_count
            call_count += 1
            output.write_text(f"{method}-{threshold}-{enabled}")

        # First call
        process(sample_files["input1"], output_file, "fast", 10, True)
        assert call_count == 1

        # Same parameters - cached
        process(sample_files["input1"], output_file, "fast", 10, True)
        assert call_count == 1

        # Different method - execute
        process(sample_files["input1"], output_file, "slow", 10, True)
        assert call_count == 2

        # Different threshold - execute
        process(sample_files["input1"], output_file, "fast", 20, True)
        assert call_count == 3

        # Different enabled - execute
        process(sample_files["input1"], output_file, "fast", 10, False)
        assert call_count == 4

    def test_parameters_with_set(self, temp_cache_dir, sample_files, tmp_path):
        """Test caching with set parameter (like excluded_seeds)."""
        call_count = 0
        output_file = tmp_path / "result.txt"

        @cache.cached("test_proc", parameters=["excluded"], outputs=["output"])
        def process(input_file: Path, output: Path, excluded: set[str]) -> None:
            nonlocal call_count
            call_count += 1
            output.write_text(f"Excluded: {sorted(excluded)}")

        # First call
        process(sample_files["input1"], output_file, {"a", "b", "c"})
        assert call_count == 1

        # Same set (different order) - should be cached
        process(sample_files["input1"], output_file, {"c", "a", "b"})
        assert call_count == 1

        # Different set - execute
        process(sample_files["input1"], output_file, {"a", "b"})
        assert call_count == 2

        # Empty set - execute
        process(sample_files["input1"], output_file, set())
        assert call_count == 3

    def test_parameters_with_none(self, temp_cache_dir, sample_files, tmp_path):
        """Test caching with None parameter values."""
        call_count = 0
        output_file = tmp_path / "result.txt"

        @cache.cached("test_proc", parameters=["optional"], outputs=["output"])
        def process(input_file: Path, output: Path, optional: str | None) -> None:
            nonlocal call_count
            call_count += 1
            output.write_text(f"Optional: {optional}")

        # First call with None
        process(sample_files["input1"], output_file, None)
        assert call_count == 1

        # Same (None) - cached
        process(sample_files["input1"], output_file, None)
        assert call_count == 1

        # Different value - execute
        process(sample_files["input1"], output_file, "value")
        assert call_count == 2

    def test_parameters_serialization(self, temp_cache_dir, sample_files, tmp_path):
        """Test that parameters are correctly serialized in cache."""
        output_file = tmp_path / "result.txt"

        @cache.cached("test_proc", parameters=["count", "name"], outputs=["output"])
        def process(input_file: Path, output: Path, count: int, name: str) -> None:
            output.write_text(f"{count}-{name}")

        process(sample_files["input1"], output_file, 42, "test")

        # Check cache data
        data = cache._load_cache_data()
        assert "test_proc" in data
        entry = data["test_proc"][0]
        inputs = entry["inputs"]

        # Parameters should be in inputs
        assert inputs["count"] == 42
        assert inputs["name"] == "test"

    def test_parameters_without_outputs(self, temp_cache_dir, sample_files):
        """Test parameters with function that returns a Path (no outputs list)."""
        call_count = 0

        @cache.cached("test_proc", parameters=["multiplier"])
        def process(input_file: Path, multiplier: int) -> Path:
            nonlocal call_count
            call_count += 1
            return sample_files["output1"]

        # First call
        process(sample_files["input1"], 2)
        assert call_count == 1

        # Same parameters - cached
        process(sample_files["input1"], 2)
        assert call_count == 1

        # Different multiplier - execute
        process(sample_files["input1"], 3)
        assert call_count == 2


class TestCachedDecoratorWithOutputs:
    """Tests for the @cached decorator with outputs argument."""

    def test_outputs_single(self, temp_cache_dir, sample_files, tmp_path):
        """Test caching with single output file."""
        call_count = 0
        output_file = tmp_path / "result.txt"

        @cache.cached("test_proc", outputs=["output"])
        def process(input_file: Path, output: Path) -> None:
            nonlocal call_count
            call_count += 1
            content = input_file.read_text()
            output.write_text(f"Processed: {content}")

        # First call
        process(sample_files["input1"], output_file)
        assert call_count == 1
        assert output_file.exists()
        original_content = output_file.read_text()

        # Delete output file
        output_file.unlink()
        assert not output_file.exists()

        # Second call - should restore from cache
        process(sample_files["input1"], output_file)
        assert call_count == 1  # Not executed again
        assert output_file.exists()
        assert output_file.read_text() == original_content

    def test_outputs_multiple(self, temp_cache_dir, sample_files, tmp_path):
        """Test caching with multiple output files."""
        call_count = 0
        output1 = tmp_path / "out1.txt"
        output2 = tmp_path / "out2.txt"

        @cache.cached("test_proc", outputs=["out1", "out2"])
        def process(input_file: Path, out1: Path, out2: Path) -> None:
            nonlocal call_count
            call_count += 1
            content = input_file.read_text()
            out1.write_text(f"Output1: {content}")
            out2.write_text(f"Output2: {content}")

        # First call
        process(sample_files["input1"], output1, output2)
        assert call_count == 1
        assert output1.exists()
        assert output2.exists()
        content1 = output1.read_text()
        content2 = output2.read_text()

        # Delete outputs
        output1.unlink()
        output2.unlink()

        # Second call - should restore both from cache
        process(sample_files["input1"], output1, output2)
        assert call_count == 1
        assert output1.read_text() == content1
        assert output2.read_text() == content2

    def test_outputs_different_locations(self, temp_cache_dir, sample_files, tmp_path):
        """Test that cached output can be written to different locations."""
        call_count = 0

        @cache.cached("test_proc", outputs=["output"])
        def process(input_file: Path, output: Path) -> None:
            nonlocal call_count
            call_count += 1
            output.write_text(f"Result: {input_file.read_text()}")

        # First call to location A
        output_a = tmp_path / "a" / "result.txt"
        output_a.parent.mkdir()
        process(sample_files["input1"], output_a)
        assert call_count == 1

        # Second call to location B - should use cache
        output_b = tmp_path / "b" / "result.txt"
        output_b.parent.mkdir()
        process(sample_files["input1"], output_b)
        assert call_count == 1  # Cached
        assert output_b.read_text() == output_a.read_text()

    def test_outputs_excludes_from_cache_key(self, temp_cache_dir, sample_files, tmp_path):
        """Test that output paths don't affect cache key."""
        call_count = 0

        @cache.cached("test_proc", outputs=["output"])
        def process(input_file: Path, output: Path) -> None:
            nonlocal call_count
            call_count += 1
            output.write_text("result")

        # Call with different output paths but same input
        output1 = tmp_path / "out1.txt"
        output2 = tmp_path / "out2.txt"

        process(sample_files["input1"], output1)
        assert call_count == 1

        process(sample_files["input1"], output2)
        assert call_count == 1  # Cached - output path doesn't matter

        # Verify both files exist with same content
        assert output1.read_text() == output2.read_text()


class TestCachedDecoratorCombined:
    """Tests for @cached decorator with both parameters and outputs."""

    def test_combined_basic(self, temp_cache_dir, sample_files, tmp_path):
        """Test caching with both parameters and outputs."""
        call_count = 0
        output_file = tmp_path / "result.txt"

        @cache.cached("test_proc",
                     parameters=["threshold"],
                     outputs=["output"])
        def process(input_file: Path, output: Path, threshold: int) -> None:
            nonlocal call_count
            call_count += 1
            output.write_text(f"Processed with {threshold}")

        # First call
        process(sample_files["input1"], output_file, 10)
        assert call_count == 1

        # Same input and parameters - cached
        process(sample_files["input1"], output_file, 10)
        assert call_count == 1

        # Different threshold - execute
        process(sample_files["input1"], output_file, 20)
        assert call_count == 2

        # Back to first threshold - should use original cache
        output_file.unlink()
        process(sample_files["input1"], output_file, 10)
        assert call_count == 2  # Cached
        assert "with 10" in output_file.read_text()

    def test_combined_complex(self, temp_cache_dir, sample_files, tmp_path):
        """Test complex caching scenario like prelim_map."""
        call_count = 0

        @cache.cached("prelim_map",
                     parameters=["nthreads", "rdgopen", "rfgopen", "gzip", "excluded_seeds"],
                     outputs=["prelim_csv"])
        def prelim_map(fastq1: Path, fastq2: Path, prelim_csv: Path,
                      nthreads: int, rdgopen: int, rfgopen: int,
                      gzip: bool, excluded_seeds: set[str] | None) -> None:
            nonlocal call_count
            call_count += 1
            result = f"Mapped {fastq1.name} and {fastq2.name}\n"
            result += f"Threads: {nthreads}, RDG: {rdgopen}, RFG: {rfgopen}\n"
            result += f"Gzip: {gzip}, Excluded: {excluded_seeds}"
            prelim_csv.write_text(result)

        output = tmp_path / "prelim.csv"

        # First call
        prelim_map(sample_files["input1"], sample_files["input2"], output,
                  nthreads=4, rdgopen=10, rfgopen=10, gzip=False, excluded_seeds={"seed1"})
        assert call_count == 1
        content1 = output.read_text()

        # Delete output
        output.unlink()

        # Same everything - cached
        prelim_map(sample_files["input1"], sample_files["input2"], output,
                  nthreads=4, rdgopen=10, rfgopen=10, gzip=False, excluded_seeds={"seed1"})
        assert call_count == 1
        assert output.read_text() == content1

        # Different nthreads - execute
        prelim_map(sample_files["input1"], sample_files["input2"], output,
                  nthreads=8, rdgopen=10, rfgopen=10, gzip=False, excluded_seeds={"seed1"})
        assert call_count == 2

        # Different excluded_seeds - execute
        prelim_map(sample_files["input1"], sample_files["input2"], output,
                  nthreads=4, rdgopen=10, rfgopen=10, gzip=False, excluded_seeds={"seed2"})
        assert call_count == 3

        # Different input file - execute
        prelim_map(sample_files["input2"], sample_files["input1"], output,
                  nthreads=4, rdgopen=10, rfgopen=10, gzip=False, excluded_seeds={"seed1"})
        assert call_count == 4

    def test_error_handling_wrong_type(self, temp_cache_dir, sample_files, tmp_path):
        """Test error handling when parameter type is wrong."""
        @cache.cached("test_proc", parameters=["threshold"], outputs=["output"])
        def process(input_file: Path, output: Path, threshold: int) -> None:
            output.write_text("result")

        output = tmp_path / "result.txt"

        # This should work
        process(sample_files["input1"], output, 10)

        # This should fail - output must be a Path when listed in outputs
        with pytest.raises(ValueError, match="Output argument .* must be a Path"):
            @cache.cached("bad", outputs=["output"])
            def bad_func(input_file: Path, output: str) -> None:
                pass
            bad_func(sample_files["input1"], "not_a_path")

    def test_error_handling_unlisted_parameter(self, temp_cache_dir, sample_files, tmp_path):
        """Test error when non-Path argument not in parameters list."""
        with pytest.raises(ValueError, match="must be either a Path, None, listed in parameters, or listed in outputs"):
            @cache.cached("test_proc", outputs=["output"])
            def process(input_file: Path, output: Path, threshold: int) -> None:
                output.write_text("result")

            output = tmp_path / "result.txt"
            # threshold is not a Path and not in parameters - should error
            process(sample_files["input1"], output, 10)


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
        # Should have one procedure
        assert procedure in data
        # With one entry
        assert len(data[procedure]) == 1

        entry = data[procedure][0]
        # Should have structured inputs
        assert "inputs" in entry
        assert "fastq1" in entry["inputs"]
        assert "fastq2" in entry["inputs"]
        # Should have outputs
        assert "outputs" in entry
        file_hash = entry["outputs"]
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


class TestClearCache:
    """Tests for the clear_cache function."""

    def test_clear_cache_no_cache_folder(self):
        """Test that clear_cache handles missing CACHE_FOLDER gracefully."""
        with patch.object(cache, "CACHE_FOLDER", None):
            count = cache.clear_cache()
            assert count == 0

    def test_clear_cache_nonexistent_folder(self, tmp_path):
        """Test that clear_cache handles nonexistent cache folder."""
        cache_dir = tmp_path / "nonexistent_cache"
        with patch.object(cache, "CACHE_FOLDER", str(cache_dir)):
            count = cache.clear_cache()
            assert count == 0

    def test_clear_entire_cache(self, temp_cache_dir, sample_files):
        """Test clearing entire cache without pattern."""
        # Set up cache with multiple procedures
        inputs1 = {"input": sample_files["input1"]}
        output1 = sample_files["output1"]
        cache.set("procedure1", inputs1, output1)

        inputs2 = {"input": sample_files["input2"]}
        output2 = sample_files["output2"]
        cache.set("procedure2", inputs2, output2)

        # Verify cache has data
        data_file = temp_cache_dir / "data.json"
        assert data_file.exists()
        with data_file.open("r") as f:
            data = json.load(f)
        assert len(data) == 2

        # Clear entire cache
        count = cache.clear_cache()
        assert count == 2

        # Verify cache is empty
        with data_file.open("r") as f:
            data = json.load(f)
        assert len(data) == 0

    def test_clear_cache_exact_match(self, temp_cache_dir, sample_files):
        """Test clearing cache with exact procedure name match."""
        # Set up cache with multiple procedures
        inputs1 = {"input": sample_files["input1"]}
        output1 = sample_files["output1"]
        cache.set("denovo[iva]", inputs1, output1)

        inputs2 = {"input": sample_files["input2"]}
        output2 = sample_files["output2"]
        cache.set("other_procedure", inputs2, output2)

        # Clear only denovo[iva]
        count = cache.clear_cache("denovo[iva]")
        assert count == 1

        # Verify only denovo[iva] was removed
        data_file = temp_cache_dir / "data.json"
        with data_file.open("r") as f:
            data = json.load(f)
        assert "denovo[iva]" not in data
        assert "other_procedure" in data

    def test_clear_cache_regex_match(self, temp_cache_dir, sample_files):
        """Test clearing cache with regex pattern."""
        # Set up cache with multiple procedures
        inputs1 = {"input": sample_files["input1"]}
        output1 = sample_files["output1"]
        cache.set("denovo[iva]", inputs1, output1)

        inputs2 = {"input": sample_files["input2"]}
        output2 = sample_files["output2"]
        cache.set("denovo[hiv]", inputs2, output2)

        # Add another procedure that doesn't match
        cache.set("other_procedure", inputs1, output1)

        # Clear all procedures starting with "denovo"
        count = cache.clear_cache("denovo.*")
        assert count == 2

        # Verify only denovo procedures were removed
        data_file = temp_cache_dir / "data.json"
        with data_file.open("r") as f:
            data = json.load(f)
        assert "denovo[iva]" not in data
        assert "denovo[hiv]" not in data
        assert "other_procedure" in data

    def test_clear_cache_no_match(self, temp_cache_dir, sample_files):
        """Test clearing cache with pattern that doesn't match."""
        # Set up cache
        inputs = {"input": sample_files["input1"]}
        output = sample_files["output1"]
        cache.set("procedure1", inputs, output)

        # Try to clear with non-matching pattern
        count = cache.clear_cache("nonexistent")
        assert count == 0

        # Verify cache still has data
        data_file = temp_cache_dir / "data.json"
        with data_file.open("r") as f:
            data = json.load(f)
        assert "procedure1" in data

    def test_clear_cache_cleans_orphaned_files(self, temp_cache_dir, sample_files):
        """Test that clearing cache also removes orphaned files."""
        # Set up cache with two procedures
        inputs1 = {"input": sample_files["input1"]}
        output1 = sample_files["output1"]
        cache.set("procedure1", inputs1, output1)

        inputs2 = {"input": sample_files["input2"]}
        output2 = sample_files["output2"]
        cache.set("procedure2", inputs2, output2)

        # Get the hash files
        data_file = temp_cache_dir / "data.json"
        with data_file.open("r") as f:
            data = json.load(f)

        hash1 = data["procedure1"][0]["outputs"]
        hash2 = data["procedure2"][0]["outputs"]

        # Verify both files exist
        assert (temp_cache_dir / hash1).exists()
        assert (temp_cache_dir / hash2).exists()

        # Clear only procedure1
        cache.clear_cache("procedure1")

        # Verify hash1 file was removed but hash2 still exists
        assert not (temp_cache_dir / hash1).exists()
        assert (temp_cache_dir / hash2).exists()


class TestCLI:
    """Tests for the CLI main function."""

    def test_main_clear_no_pattern(self, temp_cache_dir, sample_files):
        """Test CLI clear command without pattern."""
        # Set up cache
        inputs = {"input": sample_files["input1"]}
        output = sample_files["output1"]
        cache.set("procedure1", inputs, output)

        # Run CLI
        result = cache.main(["clear"])
        assert result == 0

        # Verify cache is empty
        data_file = temp_cache_dir / "data.json"
        with data_file.open("r") as f:
            data = json.load(f)
        assert len(data) == 0

    def test_main_clear_with_pattern(self, temp_cache_dir, sample_files):
        """Test CLI clear command with pattern."""
        # Set up cache
        inputs = {"input": sample_files["input1"]}
        output = sample_files["output1"]
        cache.set("denovo[iva]", inputs, output)
        cache.set("other", inputs, output)

        # Run CLI
        result = cache.main(["clear", "denovo[iva]"])
        assert result == 0

        # Verify only denovo[iva] was removed
        data_file = temp_cache_dir / "data.json"
        with data_file.open("r") as f:
            data = json.load(f)
        assert "denovo[iva]" not in data
        assert "other" in data

    def test_main_no_command(self):
        """Test CLI with no command shows help."""
        result = cache.main([])
        assert result == 1

    def test_main_invalid_command(self):
        """Test CLI with invalid command shows help."""
        # argparse will handle invalid commands by printing error and exiting
        # We just verify the function can be called
        with pytest.raises(SystemExit):
            cache.main(["invalid"])
