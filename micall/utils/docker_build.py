#!/usr/bin/env python3
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
from typing import Sequence
import os
import subprocess

import errno
import sys


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


def main(argv: Sequence[str]) -> None:
    parser = get_parser()
    args = parser.parse_args(argv)
    source_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
    version = subprocess.check_output(['git',
                                        '-c', 'core.fileMode=false',  # Ignore file mode
                                        '-c', 'safe.directory=*',  # Allow git commands in any directory
                                        'describe',
                                        '--tags',
                                        '--dirty'],
                                        cwd=source_path,
                                        text=True).strip()

    if version.startswith('v'):
        version = version[1:]

    repository_name = 'docker.illumina.com/cfe_lab/micall'
    if args.tag is not None:
        repository_name += ':' + version

    try:
        subprocess.check_call(['docker',
                               'build',
                               '--tag', repository_name,
                               '--', source_path])
    except subprocess.CalledProcessError as ex:
        if ex.returncode == errno.EPERM:
            raise PermissionError(
                'Docker build failed. Do you have root permission?') from ex
        raise
    if args.nopush:
        print("Rerun without --nopush to attempt to push the docker image to illumina and launch spacedock")
        exit()
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
