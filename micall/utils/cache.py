from typing import Mapping, Optional, Sequence, Union, Callable, TypeVar
from pathlib import Path
import os
import hashlib
import json
import shutil
from functools import wraps


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


def _make_cache_key(procedure: str, inputs: Mapping[str, Optional[Path]]) -> str:
    """Create a cache key from procedure and input file hashes.

    Args:
        procedure: The name of the procedure.
        inputs: A mapping of input identifiers to their corresponding file paths.

    Returns:
        A unique cache key string.
    """
    # Build a deterministic key from procedure and sorted input hashes
    key_parts = [procedure]

    for key in sorted(inputs.keys()):
        value = inputs[key]
        if value is None:
            key_parts.append(f"{key}:None")
        else:
            file_hash = hash_file(value)
            key_parts.append(f"{key}:{file_hash}")

    key_string = "|".join(key_parts)
    # Hash the key string to get a fixed-length key
    return hashlib.md5(key_string.encode()).hexdigest()


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

    cache_data = _load_cache_data()
    cache_key = _make_cache_key(procedure, inputs)

    if cache_key not in cache_data:
        return None

    cache_entry = cache_data[cache_key]
    cache_folder = Path(CACHE_FOLDER)

    # Check if this is a single Path or a mapping
    if isinstance(cache_entry, str):
        # Single output file
        cached_file = cache_folder / cache_entry
        if not cached_file.exists():
            return None
        return cached_file
    elif isinstance(cache_entry, dict):
        # Multiple output files
        result: dict[str, Optional[Path]] = {}
        for key, file_hash in cache_entry.items():
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

    cache_data = _load_cache_data()
    cache_key = _make_cache_key(procedure, inputs)

    # Process the output and copy files to cache
    if isinstance(output, Path):
        # Single output file
        file_hash = hash_file(output)
        cached_file = cache_folder / file_hash

        # Copy file to cache if not already there
        if not cached_file.exists():
            shutil.copy2(output, cached_file)

        cache_data[cache_key] = file_hash
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

        cache_data[cache_key] = cache_entry
    else:
        raise ValueError(f"Unsupported output type: {type(output)}")

    _save_cache_data(cache_data)


def cached(procedure: str) -> Callable[[Callable[..., T]], Callable[..., T]]:
    """Decorator for caching function results based on file path inputs."""

    def decorator(func: Callable[..., T]) -> Callable[..., T]:
        @wraps(func)
        def wrapper(*args: Sequence[object], **kwargs: dict[str, object]) -> T:
            # Build cache key from kwargs
            cache_key: dict[str, Optional[Path]] = {}

            for key, value in kwargs.items():
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
