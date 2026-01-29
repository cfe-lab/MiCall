from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter, SUPPRESS
import logging.config
import os
from pathlib import Path
from queue import Queue, Empty
from threading import Thread
from time import sleep
from typing import Tuple

from micall.monitor.sample_watcher import PipelineType
from micall.utils.version import get_version
from micall.monitor.kive_watcher import find_samples, KiveWatcher, FolderEventType, FolderEvent
from micall.monitor import update_qai
from micall.monitor.run_completion_watcher import monitor_run_completion
try:
    from micall.utils.micall_logging_override import LOGGING  # type: ignore[import]
    is_logging_override = True
except ImportError:
    from micall.utils.micall_logging_config import LOGGING
    is_logging_override = False

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
        '--denovo_main_pipeline_id',
        type=int,
        default=os.environ.get('MICALL_DENOVO_MAIN_PIPELINE_ID', None),
        help="id of general purpose denovo pipeline's container app")
    parser.add_argument(
        '--proviral_pipeline_id',
        type=int,
        default=os.environ.get('MICALL_PROVIRAL_MAIN_PIPELINE_ID', None),
        help="id of general purpose proviral pipeline's container app")
    # noinspection PyTypeChecker
    parser.add_argument(
        '--raw_data',
        type=Path,
        default=os.environ.get('MICALL_RAW_DATA',
                               Path.home() / "data/RAW_DATA"),
        help='folder to scan for raw data files, containing MiSeq/runs')
    parser.add_argument(
        '--pipeline_version',
        default=os.environ.get('MICALL_PIPELINE_VERSION', '0-dev'),
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
        default=os.environ.get('MICALL_QAI_SERVER', 'http://localhost:4567'),
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
    if args.micall_filter_quality_pipeline_id is None:
        parser.error("Argument --micall_filter_quality_pipeline_id not set and "
                     "$MICALL_FILTER_QUALITY_PIPELINE_ID environment variable "
                     "not set.")
    main_pipeline_ids = ('micall_main_pipeline_id',
                         'denovo_main_pipeline_id',
                         'mixed_hcv_pipeline_id')
    for arg_name in main_pipeline_ids:
        pipeline_id = getattr(args, arg_name)
        if pipeline_id is not None:
            break
    else:
        parser.error(f"No arguments or environment variables set for main "
                     f"pipeline ids ({', '.join(main_pipeline_ids)}).")

    return args


def main_loop(args, sample_queue, qai_upload_queue):
    kive_watcher = KiveWatcher(args, qai_upload_queue=qai_upload_queue, retry=True)
    while True:
        kive_watcher.poll_runs()
        if kive_watcher.is_full():
            sleep(POLLING_DELAY)
        else:
            try:
                while True:
                    folder_event = sample_queue.get(timeout=POLLING_DELAY)
                    if folder_event.type == FolderEventType.ADD_SAMPLE:
                        sample_watcher = kive_watcher.add_sample_group(
                            folder_event.base_calls,
                            folder_event.sample_group)
                        if sample_watcher is None or not sample_watcher.runs:
                            # It's a new sample, so go through the polling loop.
                            break
                    else:
                        kive_watcher.finish_folder(folder_event.base_calls)
            except Empty:
                if args.quit and kive_watcher.is_idle():
                    logger.info('Shutting down.')
                    break


def main():
    args = parse_args()
    logger.info('Starting up with server %s and version %s', args.kive_server, get_version())
    logger.info('Logging override: %s', is_logging_override)

    sample_queue: Queue[FolderEvent] = Queue(maxsize=2)
    qai_upload_queue: Queue[None | Tuple[Path, PipelineType]] = Queue()  # [Path] for results folders or [None] to quit.
    wait = True

    # Start the run completion watcher thread
    # This monitors MiSeq runs and creates 'needsprocessing' markers when complete
    runs_dir = args.raw_data / 'MiSeq' / 'runs'
    completion_thread = Thread(target=monitor_run_completion,
                              args=(runs_dir,),
                              daemon=True,
                              name='RunCompletionWatcher')
    completion_thread.start()

    finder_thread = Thread(target=find_samples,
                           args=(args.raw_data,
                                 args.pipeline_version,
                                 sample_queue,
                                 qai_upload_queue,
                                 wait,
                                 True),  # retry
                           daemon=True)
    finder_thread.start()

    qai_upload_thread = Thread(target=update_qai.upload_loop,
                               args=(args.qai_server,
                                     args.qai_user,
                                     args.qai_password,
                                     args.pipeline_version,
                                     qai_upload_queue),
                               daemon=True)
    qai_upload_thread.start()

    try:
        main_loop(args, sample_queue, qai_upload_queue)

        qai_upload_queue.put(None)
        qai_upload_thread.join()
    except KeyboardInterrupt:
        logger.info('Shut down requested.')


if __name__ == '__main__':
    logging.config.dictConfig(LOGGING)
    main()
