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

logger = logging.getLogger('MISEQ_MONITOR')
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
                        default=settings.kive_folder_delay,
                        help='seconds between checking for new folders')
    parser.add_argument('--retry_delay',
                        '-r',
                        type=int,
                        default=settings.kive_retry_delay,
                        help='seconds to continue retrying after error')
    return parser.parse_args()


def main():
    args = parse_args()
    logger.info('Starting up.')
    try:
        loader = KiveLoader(launch_limit=args.max,
                            status_delay=args.status_delay,
                            folder_delay=args.folder_delay,
                            retry_delay=args.retry_delay)
        while True:
            delay = loader.poll()
            sleep(delay)
    except:
        logger.error('Fatal error.', exc_info=True)
main()
