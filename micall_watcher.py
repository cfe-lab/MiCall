from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter, SUPPRESS

import os

from pathlib import Path
from queue import Queue, Empty
from threading import Thread
from time import sleep

from micall.monitor.kive_watcher import find_samples, KiveWatcher

POLLING_DELAY = 10  # seconds between scans for new samples or finished runs


def parse_args(argv=None):
    parser = ArgumentParser(description='Watch the raw data folder for new runs.',
                            formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument(
        '--micall_filter_quality_pipeline_id',
        type=int,
        default=os.environ.get('MICALL_FILTER_QUALITY_PIPELINE_ID', None),
        help="id of filter quality pipeline")
    parser.add_argument(
        '--micall_main_pipeline_id',
        type=int,
        default=os.environ.get('MICALL_MAIN_PIPELINE_ID', None),
        help="id of main pipeline")
    parser.add_argument(
        '--micall_resistance_pipeline_id',
        type=int,
        default=os.environ.get('MICALL_RESISTANCE_PIPELINE_ID', None),
        help="id of resistance pipeline")
    parser.add_argument(
        '--mixed_hcv_pipeline_id',
        type=int,
        default=os.environ.get('MICALL_MIXED_HCV_PIPELINE_ID', None),
        help='id of Mixed HCV pipeline')
    parser.add_argument(
        '--raw_data',
        type=Path,
        default=os.environ.get('MICALL_RAW_DATA',
                               Path.home() / "data/RAW_DATA"),
        help='folder to scan for raw data files, containing MiSeq/runs')
    parser.add_argument(
        '--pipeline_version',
        default='0-dev',
        help='version suffix for batch names and folder names')
    parser.add_argument(
        '--kive_server',
        default=os.environ.get('MICALL_KIVE_SERVER', 'http://localhost:8000'),
        help='server to send runs to')
    parser.add_argument(
        '--kive_user',
        default=os.environ.get('MICALL_KIVE_USER', 'kive'),
        help='user name for Kive server')
    parser.add_argument(
        '--kive_password',
        default=SUPPRESS,
        help='password for Kive server (default not shown)')
    parser.add_argument(
        '--max_active',
        type=int,
        default=os.environ.get('MICALL_MAX_ACTIVE', 5),
        help='maximum number of active jobs in Kive')

    args = parser.parse_args(argv)
    if not hasattr(args, 'kive_password'):
        args.kive_password = os.environ.get('MICALL_KIVE_PASSWORD', 'kive')
    if args.micall_filter_quality_pipeline_id is None:
        parser.error('Missing micall_filter_quality_pipeline_id. Set the '
                     'argument or environment variable.')
    return args


def main():
    args = parse_args()
    kive_watcher = KiveWatcher(args)
    sample_queue = Queue(maxsize=2)
    finder_thread = Thread(target=find_samples,
                           args=(args.raw_data, sample_queue, False),
                           daemon=True)
    finder_thread.start()
    while True:
        kive_watcher.poll_runs()
        if kive_watcher.is_full():
            sleep(POLLING_DELAY)
        else:
            try:
                base_calls, sample_group = sample_queue.get(
                    timeout=POLLING_DELAY)
                kive_watcher.add_sample_group(base_calls, sample_group)
            except Empty:
                pass


if __name__ == '__main__':
    main()
