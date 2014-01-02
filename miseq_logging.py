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

def init_logging_console_only(log_level=logging.DEBUG):
	import logging,sys
        logger = logging.getLogger()
	logger.handlers = []				# Clear previous handlers
        logger.setLevel(logging.DEBUG)
	console_logger = logging.StreamHandler(sys.stdout)
	console_logger.setLevel(log_level)
	formatter = Timestamp('%(asctime)s - [%(levelname)s] %(message)s', datefmt="%Y-%m-%d %H:%M:%S.%f")
	console_logger.setFormatter(formatter)
	logger.addHandler(console_logger)
	return logger	

def init_logging(logging_path, file_log_level=logging.DEBUG, console_log_level=logging.DEBUG):
	"""
	Creates a logger object which will log to the console and to the path specified.
	Formatting of logging is performed by an extended version of logging.Formatter

	The logging_level is assumed to be at the DEBUG level, but can be set to:
	logging.DEBUG, logging.INFO, logging.WARN, logging.ERROR, logging.CRITICAL
	"""
	import logging,sys

	logger = logging.getLogger()
	logger.handlers = []				# Clear previous handlers
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

def collate_logs(path, extension, final_log):
	"""
	Combine logs of a given extension into final_log in a sorted manner.
	Delete original source logs.
	"""

	import datetime, os
	from glob import glob

	logs = []
	for log_path in glob("{}/*.{}".format(path, extension)):
		with open(log_path, "r") as f:
			for line in f:

				# FIXME: We are LOSING improperly formed logs - we need to retain them
				try:
					date_string = line.split(" - ")[0]
					date_time = datetime.datetime.strptime(date_string, "%Y-%m-%d %H:%M:%S.%f")
					logs.append((date_time, line))
				except:
					pass
		os.remove(log_path)

	with open("{}/{}".format(path, final_log), "w") as collated_log:
		for date_time, message in logs:
			collated_log.write(message)
