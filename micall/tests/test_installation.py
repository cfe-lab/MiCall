#! /usr/bin/env python

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
import re
import os
from micall.utils.get_list_of_executables import iterate_executables


@pytest.fixture(scope="function")
def temp_venv(tmpdir: Path):
    """
    Fixture for creating and cleaning up a virtual environment.

    This fixture creates a virtual environment in a temporary directory,
    provides context to run commands in this environment, and cleans up after the test.
    """

    # Create the virtual environment
    venv_dir = tmpdir / "temp_test_venv"
    venv.create(venv_dir, with_pip=True)

    # Yield the environment setup to the test function
    yield venv_dir / "bin" / "activate"


def run_command(command: Sequence[str]):
    """Executes a shell command within a provided environment and returns output, error, and return code."""

    result = subprocess.run(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    return result.stdout.decode('utf-8').strip(), result.stderr.decode('utf-8').strip(), result.returncode


def test_micall_installation(temp_venv):
    """
    Test to verify installation of MiCall.

    This test installs MiCall in an isolated virtual environment and verifies the installation
    by executing the command `command -v micall`.
    """

    # Path to MiCall directory (3 levels up from the current script file)
    script_path = Path(__file__).resolve()
    micall_path = script_path.parent.parent.parent

    # Function to quote shell arguments.
    def q(s: object) -> str:
        return shlex.quote(str(s))

    # Check that MiCall is not installed.
    stdout, stderr, returncode = run_command(f"export PATH= ; . {q(temp_venv)} && command -v micall")
    assert returncode != 0, "Unexpected MiCall installation."

    # Install MiCall using pip from the local path
    stdout, stderr, returncode = run_command(f". {q(temp_venv)} && pip install -- {q(micall_path)}")
    assert returncode == 0, f"Failed to install MiCall:\n{stderr}"

    # Check MiCall executable path to verify installation
    stdout, stderr, returncode = run_command(f"export PATH= ; . {q(temp_venv)} && command -v micall")
    assert returncode == 0, f"Cound not find MiCall installation:\n{stderr}"
    assert stdout.endswith('micall'), "Unexpected output for micall path check."


def test_micall_version(temp_venv):
    """
    Test to verify installation of MiCall.

    This test installs MiCall in an isolated virtual environment and verifies the installation
    by executing the command `micall --version`.
    """

    # Path to MiCall directory (3 levels up from the current script file)
    script_path = Path(__file__).resolve()
    micall_path = script_path.parent.parent.parent

    # Function to quote shell arguments.
    def q(s: object) -> str:
        return shlex.quote(str(s))

    # Install MiCall using pip from the local path
    stdout, stderr, returncode = run_command(f". {q(temp_venv)} && pip install -- {q(micall_path)}")
    assert returncode == 0, f"Failed to install MiCall:\n{stderr}"

    # Check MiCall version to verify installation
    stdout, stderr, returncode = run_command(f"export PATH= ; . {q(temp_venv)} && micall --version")
    assert returncode == 0, f"MiCall version command failed:\n{stderr}"
    lines = [line.strip() for line in stdout.split('\n')]
    first_line = lines[0].strip()
    assert re.match(r'(\d+[.]\d+[.]\d+)|development', first_line), "Unexpected output for micall version check."


def test_micall_help(temp_venv):
    """
    Test to verify installation of MiCall.

    This test installs MiCall in an isolated virtual environment and verifies the installation
    by executing the command `micall --help`.
    """

    # These are supposed to be listed in output of --help.
    executables = [os.path.splitext(path.name)[0] for path in iterate_executables()]

    # Path to MiCall directory (3 levels up from the current script file)
    script_path = Path(__file__).resolve()
    micall_path = script_path.parent.parent.parent

    # Function to quote shell arguments.
    def q(s: object) -> str:
        return shlex.quote(str(s))

    # Install MiCall using pip from the local path
    stdout, stderr, returncode = run_command(f". {q(temp_venv)} && pip install -- {q(micall_path)}")
    assert returncode == 0, f"Failed to install MiCall:\n{stderr}"

    # Check MiCall help to verify installation
    stdout, stderr, returncode = run_command(f"export PATH= ; . {q(temp_venv)} && micall --help")
    assert returncode == 0, f"MiCall help command failed:\n{stderr}"

    for executable in executables:
        assert executable in stdout, f"Executable {executable!r} not listed in micall --help."
