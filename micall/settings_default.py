"""
To make pipeline portable, allow user to specify local paths and thread counts.
"""
import logging.config

pipeline_version = '7.4'        # Change for each release

instrument_number = 'M01841'  # for Illumina MiSeq, second item in run folder name
production = False  # set this to True to push results to NAS
are_temp_folders_deleted = True  # Should FIFO worker clean up working folders?

# Local path for writing data
home = '/data/miseq/'
nruns_to_store = 20  # protect X most recent runs from cleaning up intermediate files

# Scheduling processes: these should be a multiple of the total number of slots
# in your hostfile.
mapping_processes = 36
counting_processes = 144

rawdata_mount = '/media/RAW_DATA/'  # NAS
delay = 3600  # Delay (seconds) for polling NAS for unprocessed runs
NEEDS_PROCESSING = 'needsprocessing'  # File flags
ERROR_PROCESSING = 'errorprocessing'
DONE_PROCESSING = 'doneprocessing'
QC_UPLOADED = 'qc_uploaded'

# Kive settings
kive_server_url = 'http://127.0.0.1:8000/'
kive_user = 'FILLINUSERNAME'
kive_password = '*******'
kive_groups_allowed = ['Everyone']
kive_max_runs = 50  # Number of sample runs to have active at one time
kive_status_delay = 30  # seconds between checking run status
kive_folder_delay = 60*60  # seconds between scanning for new folders
kive_retry_delay = 60*60  # seconds to continue retrying after error
kive_pipelines = {000: dict(inputs=['quality', 'fastq1', 'fastq2'],
                            format='MiSeq - {sample} ({folder})')}  # Change for each release
quality_cdt_kive_id = 25        # Kive ID for CompoundDatatype (tile:integer, cycle:integer, errorrate:float?)
"""
This can be retrieved by entering the Django shell with './manage.py shell'
and using the following script:
from metadata.models import CompoundDatatype
map(lambda x: x.id, filter(lambda x: 'tile' in x.short_name, CompoundDatatype.objects.all()))
"""

# Connection to QAI RESTful API for uploading results
qai_user = "FILLINUSERNAME"
qai_password = "****"
qai_path = "http://192.168.X.Y:port"

# Connection to QAI RESTful API for dumping project configuration (read-only)
qai_project_user = "FILLINUSERNAME"
qai_project_password = "****"
qai_project_path = "http://192.168.X.Y:port"

logging.config.dictConfig({
    'version': 1,
    'formatters': {'basic': {
        'format': '%(asctime)s[%(levelname)s]%(name)s.%(funcName)s(): %(message)s',
        'datefmt': '%Y-%m-%d %H:%M:%S'}},
    'handlers': {'console': {'class': 'logging.StreamHandler',
                             'level': 'DEBUG',
                             'formatter': 'basic'},
                 'file': {'class': 'logging.handlers.RotatingFileHandler',
                          'level': 'DEBUG',
                          'formatter': 'basic',
                          'filename': '/data/miseq/micall.log',
                          'maxBytes': 1024*1024*15,  # 15MB
                          'backupCount': 10}},
    # This is the default logger.
    'root': {'handlers': ['console', 'file'],
             'level': 'WARN'},
    'loggers': {"kive_loader": {"level": "INFO"},
                "MISEQ_MONITOR": {"level": "INFO"}}})
