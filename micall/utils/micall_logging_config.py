""" Logging configuration for the MiCall pipeline.

You can override the settings in micall_logging_config.py by copying the whole
file to micall_logging_override.py. The original settings in
micall_logging_config.py should be useful defaults for developers, with extra
settings available for production use, but disabled.

For a detailed description of the settings, see the Python documentation:
https://docs.python.org/3/library/logging.config.html#logging-config-dictschema

Do not commit micall_logging_override.py to source control.
"""

# Production server probably needs /var/log/micall/micall.log
# Don't forget to create the folder and change owner to micall.
LOG_FILE = '/tmp/micall.log'

LOGGING = {
    # This is the default logger. Probably want to switch console to mail.
    'root': {'handlers': ['console', 'file'],
             'level': 'INFO'},
    'loggers': {
        "__main__": {"level": "INFO"},
        
        # Suppress routine Kive API credential refresh noise while preserving fatal errors
        "kiveapi": {"level": "ERROR"},

        # Suppress urllib3 connection retry warnings (temporary DNS/network failures)
        "urllib3.connectionpool": {"level": "ERROR"},
    },

    # This lets you call logging.getLogger() before the configuration is done.
    'disable_existing_loggers': False,

    'version': 1,
    'formatters': {'basic': {
        'format': '%(asctime)s[%(levelname)s]%(name)s.%(funcName)s(): %(message)s',
        'datefmt': '%Y-%m-%d %H:%M:%S'}},
    'filters': {
        'rate_limit': {'()': 'micall.utils.ratelimitingfilter.RateLimitingFilter',
                       'rate': 1,
                       'per': 300,
                       'burst': 5}
    },
    'handlers': {'console': {'class': 'logging.StreamHandler',
                             'level': 'DEBUG',
                             'formatter': 'basic'},
                 'file': {'class': 'logging.handlers.RotatingFileHandler',
                          'level': 'DEBUG',
                          'formatter': 'basic',
                          'filename': LOG_FILE,
                          'maxBytes': 1024*1024*15,  # 15MB
                          'backupCount': 10},
                 'mail': {'class': 'logging.handlers.SMTPHandler',
                          'filters': ['rate_limit'],
                          'level': 'WARN',
                          'formatter': 'basic',
                          'mailhost': 'localhost',  # Needs postfix to forward.
                          'fromaddr': 'no.reply.micall.server@FILLINDOMAIN.com',
                          'toaddrs': ['admin.team@FILLINDOMAIN.com'],
                          'subject': 'Error logged in MiCall Watcher'}},
}
