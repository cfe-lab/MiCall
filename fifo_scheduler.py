from Queue import Queue, Empty
import subprocess, sys
from multiprocessing.pool import ThreadPool
import threading
import traceback
import logging
from sys import exc_info
import os
import tempfile
import shutil

logger = logging.getLogger(__name__)

class Job:
    """
    A job is a shell command + optional file paths for stdout/stderr.
    If file paths are undefined, outputs default to the console.
    """

    def __init__(self,
                 shell_command=None,
                 stdout=None,
                 stderr=None,
                 script=None,
                 helpers=[],
                 args=[]):
        """ Create a new job definition.
        
        @param shell_command: a simple shell command that does not need to be
        executed in a new folder. Either pass this or script, not both.
        @param stdout: a file path to redirect standard output to
        @param stderr: a file path to redirect standard error to
        @param script: an executable script that needs to be executed in a new
        folder. A temporary folder will be created and used as the current
        directory.
        @param helpers: a sequence of helper files, like libraries or data
        files. These files will be copied into the temporary folder with the
        script.
        @param args: a sequence of arguments to pass to the script
        """
        self.shell_command = shell_command
        self.script = script
        self.stdout = stdout
        self.stderr = stderr
        self.helpers = helpers[:]
        self.args = args[:]
        
    def __repr__(self):
        return "Job(%r)" % self.shell_command
        
class Worker:
    WORKING_PATH = 'working'
    
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

        if not os.path.isdir(self.WORKING_PATH):
            os.mkdir(self.WORKING_PATH)

    def run_job(self, job):
        """Run a job's command, and open files to receive stdout/error.
        
        Prefix command with resource.
        """
        
        # If job standard in/out is specified to be sent to a file, do so - else default to sys stdout/stderr
        stdout = open(job.stdout, "a", 0) if job.stdout != None else None
        stderr = open(job.stderr, "a", 0) if job.stderr != None else None
        tempdir = None
        
        try:
            if job.shell_command is not None:
                shell_command = job.shell_command
            else:
                tempdir = tempfile.mkdtemp(dir=self.WORKING_PATH)
                shutil.copy2(job.script, tempdir)
                shell_command = os.path.join('.', os.path.basename(job.script))
                for helper in job.helpers:
                    shutil.copy2(helper, tempdir)
                for arg in job.args:
                    shell_command += ' ' + arg
    
            command = "{} {}".format(self.resource_allocated, shell_command)
            
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
                cwd=tempdir,
                stdout=stdout,
                stderr=stderr)
        finally:
            if stdout is not None:
                stdout.close()
            if stderr is not None:
                stderr.close()
            if tempdir is not None:
                shutil.rmtree(tempdir)
        
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
        """Submit a command to the pool of workers.
        
        Returns a multiprocessing.pool.AsyncResult object. Calling get() on
        that object will return 0 or raise a subprocess.CalledProcessError, 
        and successful() will return True if the return code was zero.
        """
        job = Job(shell_command, standard_out, standard_error)
        result = self.queue_job(job)
        return result
    
    def queue_job(self, job):
        """Submit a job to the pool of workers.
        
        Returns a multiprocessing.pool.AsyncResult object. Calling get() on
        that object will return 0 or raise a subprocess.CalledProcessError, 
        and successful() will return True if the return code was zero.
        """
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
        
        to_raise = None
        while True:
            is_blocking = False
            try:
                result = self.results.get(is_blocking)
            except Empty:
                if to_raise is None:
                    return
                etype, value, tb = to_raise
                raise etype, value, tb
            try:
                result.get()
            except:
                if to_raise is not None:
                    logger.error(str(to_raise[1]), exc_info=to_raise)
                to_raise = exc_info()
