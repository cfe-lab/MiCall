#!/usr/bin/env python3
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
import os
from subprocess import check_call, CalledProcessError

import errno


def parse_args():
    # noinspection PyTypeChecker
    parser = ArgumentParser(
        description="Build docker image from local source code or a tag.",
        formatter_class=ArgumentDefaultsHelpFormatter)

    parser.add_argument('-a',
                        '--agent_id',
                        default=os.environ.get('BASESPACE_AGENT_ID'),
                        help='local agent id from BaseSpace Form Builder')
    parser.add_argument('-t',
                        '--tag',
                        help='Docker tag to apply (vX.Y.Z)')
    return parser.parse_args()


def main():
    args = parse_args()
    repository_name = 'docker.illumina.com/cfe_lab/micall'
    if args.tag is not None:
        repository_name += ':' + args.tag

    source_path = os.path.abspath(os.path.dirname(__file__))
    version_filename = os.path.join(source_path, 'version.txt')
    with open(version_filename, 'w') as version_file:
        check_call(['git', 'describe', '--tags', '--dirty'],
                   cwd=source_path,
                   stdout=version_file)
    try:
        try:
            check_call(['docker',
                        'build',
                        '-t',
                        repository_name,
                        source_path])
        finally:
            os.remove(version_filename)
    except CalledProcessError as ex:
        if ex.returncode == errno.EPERM:
            raise PermissionError(
                'Docker build failed. Do you have root permission?') from ex
        raise
    check_call(['docker', 'push', repository_name])
    if args.agent_id is None:
        print('Docker image pushed. Include an agent id to launch spacedock.')
    else:
        command = ('spacedock -a ' + args.agent_id +
                   ' -m https://mission.basespace.illumina.com')
        try:
            check_call(command, shell=True)
        except KeyboardInterrupt:
            print()  # Clear the line after Ctrl-C.
            print('Shutting down.')


main()
