#! /usr/bin/env python

"""

This test is supposed to verify that installation of MiCall is not broken.

This tests assumes Debian-compatible operating system, such as Ubuntu.
It also assumes that python3 and python3-venv are installed.

It then:
  1. Creates a temporary virtual environment.
  2. Activates the environment.
  3. Installs MiCall via pip.
  4. Runs `command -v micall` to check the installation.

"""


import subprocess
import venv
import shutil
from pathlib import Path
from typing import Sequence
import pytest
import shlex


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

    # Cleanup the virtual environment after the test
    shutil.rmtree(venv_dir)


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
    stdout, stderr, returncode = run_command(f"export PATH= ; . {q(temp_venv)} && pip install -- {q(micall_path)}")
    assert returncode == 0, f"Failed to install MiCall:\n{stderr}"

    # Check MiCall version to verify installation
    stdout, stderr, returncode = run_command(f"export PATH= ; . {q(temp_venv)} && command -v micall")
    assert returncode == 0, f"MiCall version command failed:\n{stderr}"
    assert stdout.endswith('micall'), "Unexpected output for micall path check."