from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter, SUPPRESS
from functools import partial
import logging.config
import os
from pathlib import Path
from queue import Queue, Empty
from threading import Thread
from time import sleep

from micall.monitor.kive_watcher import find_samples, KiveWatcher, FolderEventType
from micall.monitor import update_qai
try:
    from micall_logging_override import LOGGING
except ImportError:
    from micall_logging_config import LOGGING

POLLING_DELAY = 10  # seconds between scans for new samples or finished runs
logger = logging.getLogger(__name__)


def parse_args(argv=None):
    parser = ArgumentParser(description='Watch the raw data folder for new runs.',
                            formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument(
        '--micall_filter_quality_pipeline_id',
        type=int,
        default=os.environ.get('MICALL_FILTER_QUALITY_PIPELINE_ID', None),
        help="id of filter quality pipeline's container app")
    parser.add_argument(
        '--micall_main_pipeline_id',
        type=int,
        default=os.environ.get('MICALL_MAIN_PIPELINE_ID', None),
        help="id of main pipeline's container app")
    parser.add_argument(
        '--micall_resistance_pipeline_id',
        type=int,
        default=os.environ.get('MICALL_RESISTANCE_PIPELINE_ID', None),
        help="id of resistance pipeline's container app")
    parser.add_argument(
        '--mixed_hcv_pipeline_id',
        type=int,
        default=os.environ.get('MICALL_MIXED_HCV_PIPELINE_ID', None),
        help="id of Mixed HCV pipeline's container app")
    parser.add_argument(
        '--denovo_pipeline_id',
        type=int,
        default=os.environ.get('MICALL_DENOVO_PIPELINE_ID', None),
        help="id of Denovo pipeline's container app")
    parser.add_argument(
        '--denovo_combined_pipeline_id',
        type=int,
        default=os.environ.get('MICALL_DENOVO_COMBINED_PIPELINE_ID', None),
        help="id of Denovo pipeline's container app for combining two samples")
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
        '--qai_server',
        default=os.environ.get('MICALL_QAI_SERVER', 'http://localhost:3000'),
        help='server to post reviews on')
    parser.add_argument(
        '--qai_user',
        default=os.environ.get('MICALL_QAI_USER', 'bob'),
        help='user name for QAI server')
    parser.add_argument(
        '--qai_password',
        default=SUPPRESS,
        help='password for QAI server (default not shown)')
    parser.add_argument(
        '--max_active',
        type=int,
        default=os.environ.get('MICALL_MAX_ACTIVE', 5),
        help='maximum number of active jobs in Kive')
    parser.add_argument(
        '--quit',
        action='store_true',
        help='quit when all runs are processed')

    args = parser.parse_args(argv)
    if not hasattr(args, 'kive_password'):
        args.kive_password = os.environ.get('MICALL_KIVE_PASSWORD', 'kive')
    if not hasattr(args, 'qai_password'):
        args.qai_password = os.environ.get('MICALL_QAI_PASSWORD', 'testing')
    if args.micall_main_pipeline_id is None:
        parser.error("Argument --micall_main_pipeline_id not set and "
                     "$MICALL_MAIN_PIPELINE_ID environment variable not set.")
    return args


def main_loop(args, sample_queue):
    result_handler = partial(update_qai.process_folder,
                             qai_server=args.qai_server,
                             qai_user=args.qai_user,
                             qai_password=args.qai_password,
                             pipeline_version=args.pipeline_version)
    kive_watcher = KiveWatcher(args, result_handler, retry=True)
    while True:
        kive_watcher.poll_runs()
        if kive_watcher.is_full():
            sleep(POLLING_DELAY)
        else:
            try:
                folder_event = sample_queue.get(timeout=POLLING_DELAY)
                if folder_event.type == FolderEventType.ADD_SAMPLE:
                    kive_watcher.add_sample_group(folder_event.base_calls,
                                                  folder_event.sample_group)
                else:
                    kive_watcher.finish_folder(folder_event.base_calls)
            except Empty:
                if args.quit and kive_watcher.is_idle():
                    logger.info('Shutting down.')
                    break


def main():
    args = parse_args()
    logger.info('Starting up with server %s', args.kive_server)

    sample_queue = Queue(maxsize=2)
    wait = True
    finder_thread = Thread(target=find_samples,
                           args=(args.raw_data,
                                 args.pipeline_version,
                                 sample_queue,
                                 wait),
                           daemon=True)
    finder_thread.start()

    try:
        main_loop(args, sample_queue)
    except KeyboardInterrupt:
        logger.info('Shutting down.')


if __name__ == '__main__':
    logging.config.dictConfig(LOGGING)
    main()
