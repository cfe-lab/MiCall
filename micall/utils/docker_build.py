#!/usr/bin/env python3
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
from typing import Sequence
import logging
import os
import subprocess

import errno
import sys

from micall.utils.externals import root_directory


logger = logging.getLogger(__name__)


def get_parser() -> ArgumentParser:
    parser = ArgumentParser(
        description="Build docker image from local source code or a tag.",
        formatter_class=ArgumentDefaultsHelpFormatter)

    verbosity_group = parser.add_mutually_exclusive_group()
    verbosity_group.add_argument('--verbose', action='store_true', help='Increase output verbosity.')
    verbosity_group.add_argument('--no-verbose', action='store_true', help='Normal output verbosity.', default=True)
    verbosity_group.add_argument('--debug', action='store_true', help='Maximum output verbosity.')
    verbosity_group.add_argument('--quiet', action='store_true', help='Minimize output verbosity.')

    parser.add_argument('-a',
                        '--agent_id',
                        default=os.environ.get('BASESPACE_AGENT_ID'),
                        help='local agent id from BaseSpace Form Builder')
    parser.add_argument('--nopush',
                        action='store_true')
    return parser


def configure_logging(args) -> None:
    if args.quiet:
        logger.setLevel(logging.ERROR)
    elif args.verbose:
        logger.setLevel(logging.INFO)
    elif args.debug:
        logger.setLevel(logging.DEBUG)
    else:
        logger.setLevel(logging.WARN)

    logging.basicConfig(
        level=logger.level,
        format='%(asctime)s[%(levelname)s]%(name)s: %(message)s',
    )


def get_latest_git_tag() -> str:
    """
    Get the version string from git describe, which is used for tagging the docker image.
    If the version is not tagged, it will be tagged with the latest tag plus a suffix indicating the number of commits since the latest tag and the short sha of the current commit.
    """

    logger.debug('Resolving latest git tag with git describe.')
    with root_directory() as source_path:
        git_root = source_path.parent
        version = subprocess.check_output(['git',
                                            '-c', 'core.fileMode=false',  # Ignore file mode
                                            '-c', 'safe.directory=*',  # Allow git commands in any directory
                                            'describe',
                                            '--tags'],
                                            cwd=git_root,
                                            text=True).strip()
    logger.debug('Resolved git tag: %s', version)

    return version


def get_repository_name() -> str:
    """
    Get the repository name for the docker image, which is based on the version from git describe.
    """

    version = get_latest_git_tag()

    if version.startswith('v'):
        version = version[1:]

    repository_name = f'docker.illumina.com/cfe_lab/micall:{version}'
    logger.debug('Resolved docker repository name: %s', repository_name)
    return repository_name


def build(verbose: bool) -> str:
    """
    Build the docker image from the local source code.
    The image will be tagged with the version from git describe, and pushed to docker.illumina.com/cfe_lab/micall.
    If the version is not tagged, it will be tagged with the latest tag plus a suffix indicating the number of commits since the latest tag and the short sha of the current commit.
    If there are uncommitted changes, the version will have a -dirty suffix.
    Returns the repository name of the built image.
    """

    repository_name = get_repository_name()
    logger.info('Building Docker image %s.', repository_name)

    try:
        with root_directory() as source_path:
            git_root = source_path.parent
            try:
                subprocess.check_call(['docker',
                                    'build',
                                    '--tag', repository_name,
                                    '--', str(git_root)],
                                    stdout=sys.stderr if not verbose else None,
                                    stderr=None if verbose else subprocess.DEVNULL,
                                    )
                logger.info('Docker build completed: %s', repository_name)
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


def main(argv: Sequence[str]) -> int:
    parser = get_parser()
    args = parser.parse_args(argv)
    configure_logging(args)
    verbose = args.verbose or args.debug

    logger.info('Starting Docker build workflow.')
    repository_name = build(verbose=verbose)

    if args.nopush:
        logger.info('Skipping docker push and spacedock launch because --nopush was provided.')
        return 0

    logger.info('Pushing Docker image: %s', repository_name)
    stderr_stream = None if verbose else subprocess.DEVNULL
    stdout_stream = stderr_stream if stderr_stream is not None else sys.stderr
    subprocess.check_call(['docker', 'push', '--', repository_name],
                          stdout=stdout_stream,
                          stderr=stderr_stream)

    if args.agent_id is None:
        logger.info('Docker image pushed. Include an agent id to launch spacedock.')
    else:
        try:
            logger.info('Launching spacedock for agent %s.', args.agent_id)
            subprocess.check_call(['spacedock', '-a', args.agent_id, '-m', 'https://mission.basespace.illumina.com'],
                                  stdout=stdout_stream,
                                  stderr=stderr_stream)
        except KeyboardInterrupt:
            logger.warning('Shutting down.')

    logger.info('Docker workflow completed successfully.')
    return 0


def entry() -> None:
    sys.exit(main(sys.argv[1:]))


if __name__ == '__main__':
    entry()
