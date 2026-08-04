"""

This test is supposed to verify that installation of MiCall is not broken.

This tests assumes Debian-compatible operating system, such as Ubuntu.
It also assumes that python3 and python3-venv are installed.

It then:
  1. Creates a temporary virtual environment.
  2. Activates the environment.
  3. Installs MiCall via pip.
  4. Runs various shell commands to check the installation.

"""


import subprocess
import venv
from pathlib import Path
from typing import Sequence
import pytest
import shlex
import os
from itertools import groupby
from packaging.version import Version, InvalidVersion

from micall.utils.get_list_of_executables import iterate_executables
from micall.__main__ import EXECUTABLES


# Function to quote shell arguments.
def quote(s: object) -> str:
    return shlex.quote(str(s))


@pytest.fixture(scope="session")
def temp_venv(tmpdir_factory):
    """
    Fixture for creating and cleaning up a virtual environment.

    This fixture creates a virtual environment in a temporary directory,
    provides context to run commands in this environment, and cleans up after the test.
    """

    # Create the virtual environment
    venv_dir = tmpdir_factory.mktemp("venv")
    venv.create(venv_dir, with_pip=True)

    # Yield the environment setup to the test function
    yield venv_dir / "bin" / "activate"


@pytest.fixture(scope="session")
def micall_installation(temp_venv: Path):
    """
    Ensures an installed micall executable.
    """

    q = quote

    # Check that MiCall is not installed.
    stdout, stderr, returncode = run_command(f"export PATH= ; . {q(temp_venv)} && command -v micall")
    assert returncode != 0, "Unexpected MiCall installation."

    # Path to MiCall directory (3 levels up from the current script file)
    script_path = Path(__file__).resolve()
    micall_path = script_path.parent.parent.parent

    # Install MiCall using pip from the local path
    stdout, stderr, returncode = run_command(f". {q(temp_venv)} && pip install -- {q(micall_path)}")
    assert returncode == 0, f"Failed to install MiCall:\n{stderr}"

    yield "micall"


def run_command(command: Sequence[str]):
    """Executes a shell command within a provided environment and returns output, error, and return code."""

    result = subprocess.run(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    return result.stdout.decode('utf-8').strip(), result.stderr.decode('utf-8').strip(), result.returncode


def test_micall_installation(temp_venv, micall_installation):
    """
    Test to verify installation of MiCall.

    This test installs MiCall in an isolated virtual environment and verifies the installation
    by executing the command `command -v micall`.
    """

    # Check MiCall executable path to verify installation
    q = quote
    stdout, stderr, returncode = run_command(f"export PATH= ; . {q(temp_venv)} && command -v micall")
    assert returncode == 0, f"Cound not find MiCall installation:\n{stderr}"
    assert stdout.endswith('micall'), "Unexpected output for micall path check."


def test_micall_version(temp_venv, micall_installation):
    """
    Test to verify installation of MiCall.

    This test installs MiCall in an isolated virtual environment and verifies the installation
    by executing the command `micall --version`.
    """

    # Check MiCall version to verify installation
    q = quote
    stdout, stderr, returncode = run_command(f"export PATH= ; . {q(temp_venv)} && micall --version")
    assert returncode == 0, f"MiCall version command failed:\n{stderr}"
    lines = [line.strip() for line in stdout.split('\n')]
    first_line = lines[0].strip()
    try:
        Version(first_line)
    except InvalidVersion:
        raise AssertionError(f"Unexpected output for micall --version:\n{stdout}")


def test_micall_help(temp_venv, micall_installation):
    """
    Test to verify installation of MiCall.

    This test installs MiCall in an isolated virtual environment and verifies the installation
    by executing the command `micall --help`.
    """

    # These are supposed to be listed in output of --help.
    executables = [os.path.splitext(path.name)[0] for path in iterate_executables()]

    # Check MiCall help to verify installation
    q = quote
    stdout, stderr, returncode = run_command(f"export PATH= ; . {q(temp_venv)} && micall --help")
    assert returncode == 0, f"MiCall help command failed:\n{stderr}"

    for executable in executables:
        assert executable in stdout, f"Executable {executable!r} not listed in micall --help."


def test_executables_names():
    """
    Verify that all and only those executables found by `iterate_executables()` are used in micall/__main__.py.
    """

    assert set(EXECUTABLES) == set(map(str, iterate_executables()))


def test_executables_duplicates():
    """
    Verify that there is no duplication in names of executables.
    """

    def get_name(path: Path) -> str:
        return os.path.splitext(path.name)[0]

    executables = list(iterate_executables())

    for key, group in groupby(executables, key=get_name):
        paths = list(map(str, group))
        assert len(paths) == 1, f"Scripts {group!r} share the same executable name."
