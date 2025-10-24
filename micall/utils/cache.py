from typing import Mapping, Optional
from pathlib import Path
import os
import hashlib


CACHE_FOLDER = os.environ.get("MICALL_CACHE_FOLDER")
HASHES: dict[Path, str] = {}


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
) -> Optional[Mapping[str, Optional[Path]]]:
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
    procedure: str, inputs: Mapping[str, Path], outputs: Mapping[str, Path]
) -> None:
    """Cache the outputs for the given inputs.

    Args:
        procedure: The name of the procedure for which to cache data.
        inputs: A mapping of input identifiers to their corresponding file paths.
        outputs: A mapping of output identifiers to their corresponding file paths.
    """

    if not CACHE_FOLDER:
        return

    raise NotImplementedError("Cache setting is not implemented yet.")
