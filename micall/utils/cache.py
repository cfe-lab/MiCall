#! /usr/bin/env python3

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


def _make_cache_key(inputs: Mapping[str, Optional[Path]]) -> dict[str, Optional[str]]:
    """Create a structured input key from input file hashes.

    Args:
        inputs: A mapping of input identifiers to their corresponding file paths.

    Returns:
        A dict mapping input names to their file hashes (or None).
    """
    input_key: dict[str, Optional[str]] = {}
    for key in sorted(inputs.keys()):
        value = inputs[key]
        if value is None:
            input_key[key] = None
        else:
            input_key[key] = hash_file(value)

    return input_key


def _find_cache_entry(
    procedure: str, input_key: dict[str, Optional[str]]
) -> Optional[dict]:
    """Find a cache entry matching the procedure and input key.

    Args:
        procedure: The procedure name.
        input_key: Dict mapping input names to their hashes.

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
    input_key: dict[str, Optional[str]],
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
    procedure: str, inputs: Mapping[str, Optional[Path]]
) -> Optional[Union[Path, Mapping[str, Optional[Path]]]]:
    """Retrieve cached data for the given inputs, if available.

    Args:
        procedure: The name of the procedure for which to retrieve cached data.
        inputs: A mapping of input identifiers to their corresponding file paths.

    Returns:
        A mapping of output identifiers to their cached file paths, or None if
        the inputs are not cached.
    """

    if not CACHE_FOLDER:
        return None

    input_key = _make_cache_key(inputs)
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
) -> None:
    """Cache the outputs for the given inputs.

    Args:
        procedure: The name of the procedure for which to cache data.
        inputs: A mapping of input identifiers to their corresponding file paths.
        output: Either a single output file path or a mapping of output identifiers.
    """

    if not CACHE_FOLDER:
        return

    cache_folder = Path(CACHE_FOLDER)
    cache_folder.mkdir(parents=True, exist_ok=True)

    input_key = _make_cache_key(inputs)

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


def cached(procedure: str) -> Callable[[Callable[..., T]], Callable[..., T]]:
    """Decorator for caching function results based on file path inputs."""

    def decorator(func: Callable[..., T]) -> Callable[..., T]:
        @wraps(func)
        def wrapper(*args: Sequence[object], **kwargs: dict[str, object]) -> T:
            # Build cache key from both args and kwargs
            # We need to map positional args to parameter names
            sig = inspect.signature(func)
            bound_args = sig.bind(*args, **kwargs)
            bound_args.apply_defaults()

            cache_key: dict[str, Optional[Path]] = {}
            for key, value in bound_args.arguments.items():
                if not isinstance(value, (Path, type(None))):
                    raise ValueError(
                        f"Input key '{key}' must be of type Path or None, got {type(value)}"
                    )
                cache_key[key] = value

            # Try to get from cache
            cached_result = get(procedure, cache_key)
            if cached_result is not None:
                return cached_result  # type: ignore

            # Execute function
            result = func(*args, **kwargs)

            # Store in cache
            set(procedure, cache_key, result)  # type: ignore

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
