import datetime
from Queue import Queue
import subprocess, sys
from multiprocessing.pool import ThreadPool
import threading

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
	def __init__(self, resource="", stderr=None, launch_callback=None):
		""" Create a worker object that can launch one process at a time.
		
		@param resource: A special command that will prefix each process
		command. For example, ssh or bpsh can be used to launch processes on
		another host.
		@param stderr: An open file if you want to redirect error reports. If
		this is None, then sys.stderr will be used.
		@param launch_callback: a method that will be called with the full
		command just before it is launched.
		"""
		
		self.resource_allocated = resource
		self.stderr = stderr
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
				self.launch_callback(command)
			
			# We use subprocess but the pipeline still expects shell=True
			# access for bash operators (>, >>, |, etc)
			# FIXME: Change to shell=False but only after the pipeline has been migrated
			try:
				returncode = subprocess.check_call(
					command, 
					shell=True,
					stdout=stdout,
					stderr=stderr)
				return returncode
			except subprocess.CalledProcessError as e:
				self._report_error(e)
				raise
		finally:
			if stdout is not None:
				stdout.close()
			if stderr is not None:
				stderr.close()
				
	def _report_error(self, exception):
		stderr = self.stderr if self.stderr is not None else sys.stderr
		curr_datetime = datetime.datetime.now()
		formatted_datetime = curr_datetime.strftime("%Y-%m-%d %H:%M:%S.%f")
		stderr.write(
            "{} - {}\n".format(
                formatted_datetime,
                exception))
		
class Factory:
	"""Factories have a queue of jobs and workers to work on them"""	

	def __init__(self, assigned_resources=[], stderr=None, launch_callback=None):
		""" Create a factory with the requested resources.
		
		@param assigned_resources: a list of (resource, worker_count) tuples.
		resource is passed to each worker, and may be the empty string if there
		are no special requirements for launching jobs.
		@param stderr: A file object that lets you redirect any error messages 
		*from the workers*, not from the jobs.
		@param launch_callback: A callback method that is called with a
		(process, command) pair every time a job is started.
		"""
		self.workers = []
		self.local_data = threading.local()

		for resource, number in assigned_resources:
			for _ in range(number):
				self.workers.append(Worker(
					resource,
					stderr=stderr,
					launch_callback=launch_callback))
		self._create_pool()
		
	def _create_pool(self):
		queue = Queue()
		for worker in self.workers:
			queue.put(worker)
		self.workerQueue = queue
		self.pool = ThreadPool(len(self.workers))

	def queue_work(self, shell_command, standard_out=None, standard_error=None):
		"""Submit a job to the pool of workers.
		
		Returns a multiprocessing.pool.AsyncResult object. Calling get() on
		that object will return 0 or raise a subprocess.CalledProcessError, 
		and successful() will return True if the return code was zero.
		"""
		job = Job(shell_command, standard_out, standard_error)
		return self.pool.apply_async(self._run_job, (job,))
	
	def _run_job(self, job):
		if not hasattr(self.local_data, 'worker'):
			self.local_data.worker = self.workerQueue.get()
		return self.local_data.worker.run_job(job)

	def wait(self, sleep_seconds=1):
		"""Wait for all queued jobs to complete.
		
		No new jobs can be queued while this method is waiting, and it's not
		safe to have two threads waiting on the same factory object."""
		
		self.pool.close()
		self.pool.join()
		self._create_pool()
