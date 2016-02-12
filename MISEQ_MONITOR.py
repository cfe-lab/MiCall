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


if sys.version_info[:2] != (2, 7):
    raise Exception("Python 2.7 not detected")


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
    logger.info('Starting up.')
    try:
        loader = KiveLoader()
        while True:
            loader.poll()
            sleep(5)
    except:
        logger.error('Fatal error.', exc_info=True)
main()
