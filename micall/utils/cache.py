from typing import Mapping, Optional, Sequence, Union, Callable, TypeVar
from pathlib import Path
import os
import hashlib
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

    raise NotImplementedError("Cache retrieval is not implemented yet.")


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

    raise NotImplementedError("Cache setting is not implemented yet.")


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
