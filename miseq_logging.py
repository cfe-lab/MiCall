import logging

class Timestamp(logging.Formatter):
	"""
	Extended logging.Formatter to use DateTime instead of struct_time which doesn't support milliseconds.
	http://stackoverflow.com/questions/6290739/python-logging-use-milliseconds-in-time-format
	"""
	import datetime

	converter=datetime.datetime.fromtimestamp

	def formatTime(self, record, datefmt=None):
		ct = self.converter(record.created)
		if datefmt:
			s = ct.strftime(datefmt)
		else:
			t = ct.strftime("%Y-%m-%d %H:%M:%S")
			s = "%s,%03d" % (t, record.msecs)
		return s


def init_logging(logging_path, file_log_level=logging.DEBUG, console_log_level=logging.DEBUG):
	"""
	Creates a logger object which will log to the console and to the path specified.
	Formatting of logging is performed by an extended version of logging.Formatter

	The logging_level is assumed to be at the DEBUG level, but can be set to:
	logging.DEBUG, logging.INFO, logging.WARN, logging.ERROR, logging.CRITICAL
	"""
	import sys

	logger = logging.getLogger()
	logger.setLevel(logging.DEBUG)

	# Logging will go to 2 different places (The console and a log file)
	file_logger = logging.FileHandler(logging_path)
	file_logger.setLevel(file_log_level)
	console_logger = logging.StreamHandler(sys.stdout)
	console_logger.setLevel(console_log_level)

	# Format the handlers
	formatter = Timestamp('%(asctime)s - [%(levelname)s] %(message)s', datefmt="%Y-%m-%d %H:%M:%S.%f")
	console_logger.setFormatter(formatter)
	file_logger.setFormatter(formatter)
	logger.addHandler(console_logger)
	logger.addHandler(file_logger)

	return logger
