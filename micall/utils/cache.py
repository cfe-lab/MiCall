#! /usr/bin/env python3

"""
Disk-based caching system module.

This module provides a persistent caching mechanism that stores function results
on disk, enabling expensive computations to be skipped when the same inputs are
encountered again. The cache is content-addressable, meaning files are identified
by their content hash rather than their path.

The cache requires the `MICALL_CACHE_FOLDER` environment variable to be set.
If not configured, all cache operations will be no-ops.

Command-line Interface:
    The module provides CLI commands for cache management:

    - `micall cache clear`: Remove all cached entries
    - `micall cache clear <procedure>`: Remove entries for a specific procedure
    - `micall cache clear <pattern>`: Remove entries matching a regex pattern

Example:
    >>> from micall.utils.cache import cached
    >>> from pathlib import Path
    >>>
    >>> @cached("my_analysis")
    >>> def analyze_sequence(input_file: Path) -> Path:
    ...     # Expensive computation here
    ...     result = process(input_file)
    ...     return result
    >>>
    >>> # First call computes and caches
    >>> result1 = analyze_sequence(Path("data.fastq"))
    >>>
    >>> # Second call with same input returns cached result instantly
    >>> result2 = analyze_sequence(Path("data.fastq"))
"""

import argparse
import sys
from typing import Mapping, Optional, Sequence, Union, Callable, TypeVar
from pathlib import Path
import os
import hashlib
import json
import shutil
from functools import wraps
import inspect
import re
import builtins


CACHE_FOLDER = os.environ.get("MICALL_CACHE_FOLDER")
HASHES: dict[Path, str] = {}

T = TypeVar("T")


def hash_file(file_path: Path) -> str:
    """Compute an MD5 hash of the given file.

    Args:
        file_path: Path to the file to hash.

    Returns:
        The MD5 hash as a hexadecimal string.
    """

    existing = HASHES.get(file_path)
    if existing:
        return existing

    hasher = hashlib.md5()
    with file_path.open("rb") as f:
        chunk: bytes = b"111"
        while chunk:
            chunk = f.read(8192)
            hasher.update(chunk)

    ret = hasher.hexdigest()
    HASHES[file_path] = ret
    return ret


def _load_cache_data() -> dict:
    """Load the cache data from data.json."""

    if not CACHE_FOLDER:
        return {}

    cache_folder = Path(CACHE_FOLDER)
    data_file = cache_folder / "data.json"

    if not data_file.exists():
        return {}

    with data_file.open("r") as f:
        return json.load(f)


def _save_cache_data(data: dict) -> None:
    """Save the cache data to data.json."""
    if not CACHE_FOLDER:
        return

    cache_folder = Path(CACHE_FOLDER)
    cache_folder.mkdir(parents=True, exist_ok=True)

    data_file = cache_folder / "data.json"

    with data_file.open("w") as f:
        json.dump(data, f, indent="\t")


def _make_cache_key(
    inputs: Mapping[str, Optional[Path]],
    parameters: Mapping[str, object] = {},
) -> dict[str, Optional[str] | object]:
    """Create a structured input key from input file hashes and serializable parameters.

    Args:
        inputs: A mapping of input identifiers to their corresponding file paths.
        parameters: Mapping of parameter names to their values (for non-Path arguments).

    Returns:
        A dict mapping input names to their file hashes (or None) and parameter names to their values.
    """
    input_key: dict[str, Optional[str] | object] = {}

    # Add file hashes
    for key in sorted(inputs.keys()):
        value = inputs[key]
        if value is None:
            input_key[key] = None
        else:
            input_key[key] = hash_file(value)

    # Add serializable parameters
    for key in sorted(parameters.keys()):
        param_value = parameters[key]
        # Convert sets to sorted lists for consistent serialization
        if isinstance(param_value, builtins.set):
            input_key[key] = sorted(list(param_value))
        else:
            input_key[key] = param_value

    return input_key


def _find_cache_entry(
    procedure: str, input_key: dict[str, Optional[str] | object]
) -> Optional[dict]:
    """Find a cache entry matching the procedure and input key.

    Args:
        procedure: The procedure name.
        input_key: Dict mapping input names to their hashes or parameter values.

    Returns:
        The cache entry dict if found, None otherwise.
    """
    cache_data = _load_cache_data()

    if procedure not in cache_data:
        return None

    # Search through all entries for this procedure
    for entry in cache_data[procedure]:
        if entry.get("inputs") == input_key:
            return entry

    return None


def _add_cache_entry(
    procedure: str,
    input_key: dict[str, Optional[str] | object],
    outputs: Union[str, dict[str, Optional[str]]],
) -> None:
    """Add or update a cache entry.

    Args:
        procedure: The procedure name.
        input_key: Dict mapping input names to their hashes.
        outputs: Either a single hash or dict of output hashes.
    """
    cache_data = _load_cache_data()

    if procedure not in cache_data:
        cache_data[procedure] = []

    # Remove any existing entry with the same inputs
    cache_data[procedure] = [
        entry for entry in cache_data[procedure] if entry.get("inputs") != input_key
    ]

    # Add the new entry
    cache_data[procedure].append({"inputs": input_key, "outputs": outputs})

    _save_cache_data(cache_data)


def get(
    procedure: str,
    inputs: Mapping[str, Optional[Path]],
    parameters: Mapping[str, object] = {},
) -> Optional[Union[Path, Mapping[str, Optional[Path]]]]:
    """Retrieve cached data for the given inputs and parameters, if available.

    Args:
        procedure: The name of the procedure for which to retrieve cached data.
        inputs: A mapping of input identifiers to their corresponding file paths.
        parameters: Mapping of parameter names to their values.

    Returns:
        A mapping of output identifiers to their cached file paths, or None if
        the inputs are not cached.
    """

    if not CACHE_FOLDER:
        return None

    input_key = _make_cache_key(inputs, parameters)
    entry = _find_cache_entry(procedure, input_key)

    if entry is None:
        return None

    outputs = entry.get("outputs")
    cache_folder = Path(CACHE_FOLDER)

    # Check if this is a single Path or a mapping
    if isinstance(outputs, str):
        # Single output file
        cached_file = cache_folder / outputs
        if not cached_file.exists():
            return None
        return cached_file
    elif isinstance(outputs, dict):
        # Multiple output files
        result: dict[str, Optional[Path]] = {}
        for key, file_hash in outputs.items():
            if file_hash is None:
                result[key] = None
            else:
                cached_file = cache_folder / file_hash
                if not cached_file.exists():
                    return None
                result[key] = cached_file
        return result

    return None


def set(
    procedure: str,
    inputs: Mapping[str, Optional[Path]],
    output: Union[Path, Mapping[str, Optional[Path]]],
    parameters: Mapping[str, object] = {},
) -> None:
    """Cache the outputs for the given inputs and parameters.

    Args:
        procedure: The name of the procedure for which to cache data.
        inputs: A mapping of input identifiers to their corresponding file paths.
        output: Either a single output file path or a mapping of output identifiers.
        parameters: Mapping of parameter names to their values.
    """

    if not CACHE_FOLDER:
        return

    cache_folder = Path(CACHE_FOLDER)
    cache_folder.mkdir(parents=True, exist_ok=True)

    input_key = _make_cache_key(inputs, parameters)

    # Process the output and copy files to cache
    if isinstance(output, Path):
        # Single output file
        file_hash = hash_file(output)
        cached_file = cache_folder / file_hash

        # Copy file to cache if not already there
        if not cached_file.exists():
            shutil.copy2(output, cached_file)

        _add_cache_entry(procedure, input_key, file_hash)
    elif isinstance(output, dict):
        # Multiple output files
        cache_entry: dict[str, Optional[str]] = {}
        for key, file_path in output.items():
            if file_path is None:
                cache_entry[key] = None
            else:
                file_hash = hash_file(file_path)
                cached_file = cache_folder / file_hash

                # Copy file to cache if not already there
                if not cached_file.exists():
                    shutil.copy2(file_path, cached_file)

                cache_entry[key] = file_hash

        _add_cache_entry(procedure, input_key, cache_entry)
    else:
        raise ValueError(f"Unsupported output type: {type(output)}")


def cached(
    procedure: str,
    parameters: Sequence[str] = (),
    outputs: Optional[Sequence[str]] = None,
) -> Callable[[Callable[..., T]], Callable[..., T]]:
    """Decorator for caching function results based on file path inputs and parameters.

    Args:
        procedure: Name of the cached procedure.
        parameters: Optional list of parameter names that should be serialized into the cache key
                   instead of being treated as Paths. These parameters will have their values
                   directly added to the cache key.
        outputs: Optional list of Path parameter names that are outputs (not inputs).
                These are excluded from the cache key and are where results are written.

    Example:
        @cached("my_proc", parameters=['threshold'], outputs=['output_file'])
        def process(input_file: Path, output_file: Path, threshold: int) -> None:
            # input_file is hashed for cache key
            # output_file is where results are written (excluded from cache key)
            # threshold is serialized as a parameter
            pass
    """

    def decorator(func: Callable[..., T]) -> Callable[..., T]:
        @wraps(func)
        def wrapper(*args: Sequence[object], **kwargs: dict[str, object]) -> T:
            # Build cache key from both args and kwargs
            # We need to map positional args to parameter names
            sig = inspect.signature(func)
            bound_args = sig.bind(*args, **kwargs)
            bound_args.apply_defaults()

            # Separate Path inputs from serializable parameters and outputs
            path_inputs: dict[str, Optional[Path]] = {}
            param_values: dict[str, object] = {}
            output_paths: dict[str, Path] = {}
            parameter_set = builtins.set(parameters)
            output_set = builtins.set(outputs) if outputs else builtins.set()

            for key, value in bound_args.arguments.items():
                if key in output_set:
                    # This is an output Path
                    if not isinstance(value, Path):
                        raise ValueError(
                            f"Output argument '{key}' must be a Path, got {type(value).__name__}"
                        )
                    output_paths[key] = value
                elif key in parameter_set:
                    # This is a serializable parameter
                    param_values[key] = value
                elif isinstance(value, Path):
                    # This is an input Path
                    path_inputs[key] = value
                elif value is None:
                    # None is allowed for optional Path arguments
                    path_inputs[key] = None
                else:
                    # Not a Path and not in parameters/outputs list - error
                    raise ValueError(
                        f"Argument '{key}' must be either a Path, None, listed in parameters, or listed in outputs. "
                        f"Got type {type(value).__name__}"
                    )

            # Try to get from cache
            cached_result = get(procedure, path_inputs, param_values)
            if cached_result is not None:
                # If no outputs specified, return cached result directly
                if not output_paths:
                    return cached_result  # type: ignore

                # Copy cached result to output location(s)
                if len(output_paths) == 1:
                    # Single output
                    output_path = list(output_paths.values())[0]
                    if isinstance(cached_result, Path):
                        shutil.copy2(cached_result, output_path)
                    else:
                        raise ValueError(f"Expected cached result to be a Path, got {type(cached_result)}")
                elif len(output_paths) > 1:
                    # Multiple outputs
                    if not isinstance(cached_result, dict):
                        raise ValueError(f"Expected cached result to be a dict, got {type(cached_result)}")
                    for key, output_path in output_paths.items():
                        cached_file = cached_result.get(key)
                        if cached_file:
                            shutil.copy2(cached_file, output_path)

                # Return the original return value of the function
                return None  # type: ignore

            # Execute function
            result = func(*args, **kwargs)

            # Store output file(s) in cache
            output_to_cache: Union[Path, Mapping[str, Path]]
            if len(output_paths) == 1:
                output_to_cache = list(output_paths.values())[0]
            elif len(output_paths) > 1:
                output_to_cache = output_paths
            else:
                # No outputs specified, use the function's return value
                output_to_cache = result  # type: ignore

            set(procedure, path_inputs, output_to_cache, param_values)

            return result

        return wrapper

    return decorator


def clear_cache(pattern: Optional[str] = None) -> int:
    """Clear cache entries matching the given pattern.

    Args:
        pattern: Optional pattern to match procedure names.
                 If None, clears all cache.
                 Can be an exact match or a regex pattern.

    Returns:
        Number of entries cleared.
    """
    qpattern = json.dumps(pattern)

    if not CACHE_FOLDER:
        print("No cache folder configured (MICALL_CACHE_FOLDER not set).")
        return 0

    cache_folder = Path(CACHE_FOLDER)

    if not cache_folder.exists():
        print("Cache folder does not exist.")
        return 0

    cache_data = _load_cache_data()

    if not pattern:
        # Clear entire cache
        count = sum(len(entries) for entries in cache_data.values())

        # Remove all files in cache folder
        for item in cache_folder.iterdir():
            if item.is_file():
                item.unlink()
            elif item.is_dir():
                shutil.rmtree(item)

        # Clear data.json
        _save_cache_data({})

        print(f"Cleared entire cache ({count} entries).")
        return count

    # Clear by pattern (exact match or regex)
    cleared_count = 0
    procedures_to_remove = []

    for procedure in cache_data.keys():
        # Try exact match first
        if procedure == pattern:
            cleared_count += len(cache_data[procedure])
            procedures_to_remove.append(procedure)
        else:
            # Try regex match
            try:
                if re.match(pattern, procedure):
                    cleared_count += len(cache_data[procedure])
                    procedures_to_remove.append(procedure)
            except re.error:
                # Not a valid regex, skip
                pass

    # Remove matched procedures
    for procedure in procedures_to_remove:
        del cache_data[procedure]

    # Save updated cache data
    _save_cache_data(cache_data)

    # Clean up orphaned files (files no longer referenced in data.json)
    if cache_data:
        referenced_hashes: builtins.set[str] = builtins.set()
        for entries in cache_data.values():
            for entry in entries:
                outputs = entry.get("outputs")
                if isinstance(outputs, str):
                    referenced_hashes.add(outputs)
                elif isinstance(outputs, dict):
                    for file_hash in outputs.values():
                        if file_hash:
                            referenced_hashes.add(file_hash)

        # Remove unreferenced files
        for item in cache_folder.iterdir():
            if item.is_file() and item.name != "data.json":
                if item.name not in referenced_hashes:
                    item.unlink()

    if cleared_count > 0:
        print(f"Cleared {cleared_count} cache entries matching {qpattern}.")
    else:
        print(f"No cache entries found matching {qpattern}.")

    return cleared_count


def main(argv: Sequence[str]) -> int:
    parser = argparse.ArgumentParser(
        description="Utility to manage MiCall cache.", prog="micall cache"
    )

    subparsers = parser.add_subparsers(dest="command", help="Cache commands")

    # Clear command
    clear_parser = subparsers.add_parser("clear", help="Clear cache entries")
    clear_parser.add_argument(
        "pattern",
        nargs="?",
        help="Optional pattern to match procedure names (exact or regex). If omitted, clears entire cache.",
    )

    args = parser.parse_args(argv)

    if args.command == "clear":
        clear_cache(args.pattern)
        return 0
    else:
        parser.print_help()
        return 1

    return 0


def entry():
    sys.exit(main(sys.argv[1:]))


if __name__ == "__main__": entry()  # noqa
