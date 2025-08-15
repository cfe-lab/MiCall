
from importlib.metadata import version
from typing import Sequence
from functools import cache


@cache
def get_version() -> str:
    full_package_name = __package__
    if full_package_name is None:
        return "development"
    else:
        root_package_name = full_package_name.split('.')[0]
        return str(version(root_package_name))


def main(argv: Sequence[str]) -> int:
    print(get_version())
    return 0


if __name__ == "__main__":
    import sys
    sys.exit(main(sys.argv[1:]))
