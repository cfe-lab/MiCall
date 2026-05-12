#!/usr/bin/env python3
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
from typing import Sequence
import os
import subprocess

import errno
import sys

from micall.utils.externals import root_directory


def get_parser() -> ArgumentParser:
    parser = ArgumentParser(
        description="Build docker image from local source code or a tag.",
        formatter_class=ArgumentDefaultsHelpFormatter)

    parser.add_argument('-a',
                        '--agent_id',
                        default=os.environ.get('BASESPACE_AGENT_ID'),
                        help='local agent id from BaseSpace Form Builder')
    parser.add_argument('--nopush',
                        action='store_true')
    return parser


def get_repository_name() -> str:
    """
    Get the repository name for the docker image, which is based on the version from git describe.
    """

    with root_directory() as source_path:
        version = subprocess.check_output(['git',
                                            '-c', 'core.fileMode=false',  # Ignore file mode
                                            '-c', 'safe.directory=*',  # Allow git commands in any directory
                                            'describe',
                                            '--tags'],
                                            cwd=source_path,
                                            text=True).strip()

    if version.startswith('v'):
        version = version[1:]

    return f'docker.illumina.com/cfe_lab/micall:{version}'


def build() -> str:
    """
    Build the docker image from the local source code.
    The image will be tagged with the version from git describe, and pushed to docker.illumina.com/cfe_lab/micall.
    If the version is not tagged, it will be tagged with the latest tag plus a suffix indicating the number of commits since the latest tag and the short sha of the current commit.
    If there are uncommitted changes, the version will have a -dirty suffix.
    Returns the repository name of the built image.
    """

    repository_name = get_repository_name()

    try:
        with root_directory() as source_path:
            try:
                subprocess.check_call(['docker',
                                    'build',
                                    '--tag', repository_name,
                                    '--', str(source_path)])
            except subprocess.CalledProcessError as ex:
                if ex.returncode == errno.EPERM:
                    raise PermissionError(
                        'Docker build failed. Do you have root permission?') from ex
                raise
    except subprocess.CalledProcessError as ex:
        if ex.returncode == errno.EPERM:
            raise PermissionError(
                'Docker build failed. Do you have root permission?') from ex
        raise

    return repository_name


def main(argv: Sequence[str]) -> None:
    parser = get_parser()
    args = parser.parse_args(argv)
    build()
    if args.nopush:
        print("Rerun without --nopush to attempt to push the docker image to illumina and launch spacedock")
        exit()
    repository_name = get_repository_name()
    subprocess.check_call(['docker', 'push', '--', repository_name])
    if args.agent_id is None:
        print('Docker image pushed. Include an agent id to launch spacedock.')
    else:
        try:
            subprocess.check_call(['spacedock', '-a', args.agent_id, '-m', 'https://mission.basespace.illumina.com'])
        except KeyboardInterrupt:
            print()  # Clear the line after Ctrl-C.
            print('Shutting down.')


def entry() -> None:
    main(sys.argv[1:])


if __name__ == '__main__':
    entry()
