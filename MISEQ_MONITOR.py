#! /usr/bin/env python

"""
MISEQ_MONITOR.py
1) For runs flagged 'needsprocessing' that have not yet been processed, copy fastqs to local disk
2) Process the run by calling MISEQ_MONITOR.py
3) Upload results to the network drive
"""

import logging
import sys
from time import sleep

import micall.settings as settings
from micall.monitor.kive_loader import KiveLoader
from argparse import ArgumentParser


if sys.version_info[:2] != (2, 7):
    raise Exception("Python 2.7 not detected")


def parse_args():
    parser = ArgumentParser(description='Process MiSeq data on Kive.')
    parser.add_argument('--max',
                        '-x',
                        type=int,
                        default=settings.kive_max_runs,
                        help='maximum number of active samples on Kive')
    parser.add_argument('--status_delay',
                        '-s',
                        type=int,
                        default=settings.kive_status_delay,
                        help='seconds between checking run status')
    parser.add_argument('--folder_delay',
                        '-f',
                        type=int,
                        default=settings.kive_status_delay,
                        help='seconds between checking for new folders')
    return parser.parse_args()


def mark_run_as_disabled(root, message, exc_info=None):
    """ Mark a run that failed, so it won't be processed again.

    @param root: path to the run folder that had an error
    @param message: a description of the error
    @param exc_info: details about the error's exception in the standard tuple,
        True to look up the current exception, or None if there is no exception
        to report
    """
    failure_message = message + " - skipping run " + root
    logger.error(failure_message, exc_info=exc_info)
    if settings.production:
        with open(root + settings.ERROR_PROCESSING, 'w') as f:
            f.write(message)
    else:
        # in development mode - exit the monitor if a run fails
        sys.exit()
    return failure_message

logger = logging.getLogger('MISEQ_MONITOR')


def main():
    args = parse_args()
    logger.info('Starting up.')
    try:
        loader = KiveLoader(launch_limit=args.max,
                            status_delay=args.status_delay,
                            folder_delay=args.folder_delay)
        while True:
            delay = loader.poll()
            sleep(delay)
    except:
        logger.error('Fatal error.', exc_info=True)
main()
