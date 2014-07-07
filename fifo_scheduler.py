from Queue import Queue, Empty
import subprocess, sys
from multiprocessing.pool import ThreadPool
import threading
import traceback

class Job:
	"""
	A job is a shell command + optional file paths for stdout/stderr.
	If file paths are undefined, outputs default to the console.
	"""

	def __init__(self, shell_command, standard_output=None, standard_error=None):
		self.shell_command = shell_command
		self.standard_output = standard_output
		self.standard_error = standard_error
		
	def __repr__(self):
		return "Job(%r)" % self.shell_command
		
class Worker:
	def __init__(self, resource="", launch_callback=None):
		""" Create a worker object that can launch one process at a time.
		
		@param resource: A special command that will prefix each process
		command. For example, ssh or bpsh can be used to launch processes on
		another host.
		@param launch_callback: a method that will be called with the full
		command just before it is launched.
		"""
		
		self.resource_allocated = resource
		self.launch_callback = launch_callback

	def run_job(self, job):
		"""Run a job's command, and open files to receive stdout/error.
		
		Prefix command with resource.
		"""
		
		command = "{} {}".format(self.resource_allocated, job.shell_command)
		
		# If job standard in/out is specified to be sent to a file, do so - else default to sys stdout/stderr
		stdout = open(job.standard_output, "a", 0) if job.standard_output != None else None
		stderr = open(job.standard_error, "a", 0) if job.standard_error != None else None
		
		try:
			if self.launch_callback is not None:
				try:
					self.launch_callback(command)
				except:
					etype, value, _ = sys.exc_info()
					original_message = traceback.format_exception_only(
						etype, 
						value)[-1].strip()
					message = "Failed callback for {!r}. {}".format(
						command,
						original_message)
					raise RuntimeError(message)
			
			# We use subprocess but the pipeline still expects shell=True
			# access for bash operators (>, >>, |, etc)
			# FIXME: Change to shell=False but only after the pipeline has been migrated
			return subprocess.check_call(
				command, 
				shell=True,
				stdout=stdout,
				stderr=stderr)
		finally:
			if stdout is not None:
				stdout.close()
			if stderr is not None:
				stderr.close()
		
class Factory:
	"""Factories have a queue of jobs and workers to work on them"""	

	def __init__(self, assigned_resources=[], launch_callback=None):
		""" Create a factory with the requested resources.
		
		@param assigned_resources: a list of (resource, worker_count) tuples.
		resource is passed to each worker, and may be the empty string if there
		are no special requirements for launching jobs.
		@param launch_callback: A callback method that is called with a
		(process, command) pair every time a job is started.
		"""
		self.workers = []
		self.workerQueue = Queue()
		self.results = Queue()
		self.local_data = threading.local()

		for resource, number in assigned_resources:
			for _ in range(number):
				worker = Worker(
				    resource, 
				    launch_callback=launch_callback)
				self.workers.append(worker)
				self.workerQueue.put(worker)
		self.pool = ThreadPool(len(self.workers))

	def queue_work(self, shell_command, standard_out=None, standard_error=None):
		"""Submit a job to the pool of workers.
		
		Returns a multiprocessing.pool.AsyncResult object. Calling get() on
		that object will return 0 or raise a subprocess.CalledProcessError, 
		and successful() will return True if the return code was zero.
		"""
		job = Job(shell_command, standard_out, standard_error)
		result = self.pool.apply_async(self._run_job, (job, ))
		self.results.put(result)
		return result
	
	def _run_job(self, job):
		if not hasattr(self.local_data, 'worker'):
			self.local_data.worker = self.workerQueue.get()
		return self.local_data.worker.run_job(job)

	def wait(self):
		"""Wait for all queued jobs to complete and check for failures.
		
		Raise the last exception encountered, with any other exceptions logged
		at the error level."""
		
		toRaise = None
		while True:
			is_blocking = False
			try:
				result = self.results.get(is_blocking)
			except Empty:
				if toRaise is None:
					return
				raise toRaise
			try:
				result.get()
			except Exception as e:
				toRaise = e
