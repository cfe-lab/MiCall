import subprocess, sys
from Queue import Queue

class Job:
	"""
	A job is composed of a shell command, and the paths to
	which standard_output and standard_error are saved.
	If either path is undefined, the default is console.
	"""

	def __init__(self, shell_command, standard_output=None, standard_error=None):
		self.shell_command = shell_command
		self.standard_output = None if standard_output == None else open(standard_output, "wb")
		self.standard_error = None if standard_error == None else open(standard_error,"wb")


class Job_Queue:
	"""FIFO queue of jobs to be completed"""

	def __init__(self):
		self.jobs = Queue()

	def add_job(self,job):
		self.jobs.put(job)

	def num_jobs(self):
		return int(self.jobs.qsize())

	def jobs_outstanding(self):
		return self.num_jobs() != 0

 	def get_next_job(self):
		"""Return next Job + remove it from the queue"""
 		if self.jobs_outstanding(): return self.jobs.get()
		raise Exception("No jobs in queue")

class Worker:

	def __init__(self, resource=""):
		self.resource_allocated = resource
		self.process = None

 	def available_for_work(self):
		if self.process == None or self.process.poll() == 0: return True
		return False

	def start_job(self,job):
		if not self.available_for_work(): raise Exception('Worker already allocated')

		command = "{} {}".format(self.resource_allocated, job.shell_command)
		standard_out_f = job.standard_output if job.standard_output != None else sys.stdout
		standard_error_f = job.standard_error if job.standard_error != None else sys.stderr
		self.process = subprocess.Popen(job.shell_command, shell=True, stdout=standard_out_f, stderr=standard_error_f)

		return self.process

class Factory:

	def __init__(self):
		self.workers = []
		self.jobs = Job_Queue()

	def hire_worker(self,resource=""):
		"""Give factory additional resources"""
		self.workers.append(Worker(resource))

	def fire_worker(index):
		"""Remove a resource"""
		del self.workers[index]

	def queue_work(self, shell_command, standard_out=None, standard_error=None):
		new_job = Job(shell_command, standard_out, standard_error)
		for worker in self.workers:
			if worker.available_for_work():
				return worker.start_job(new_job)
		self.jobs.add_job(new_job)

	def supervise(self):
		processes = []
		for worker in self.workers:
			if worker.available_for_work() and self.jobs.jobs_outstanding():
				processes.append(worker.start_job(self.jobs.get_next_job()))
		return processes

	def completely_idle(self):
		if self.jobs.jobs_outstanding(): return False
		for worker in self.workers:
			if not worker.available_for_work():
				return False
		return True
