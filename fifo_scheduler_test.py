import unittest

from fifo_scheduler import Job, Job_Queue, Worker, Factory
from cStringIO import StringIO
import os
import time

def purge_files(*filenames):
    for filename in filenames:
        try:
            os.remove(filename)
        except OSError:
            pass
    
class JobTest(unittest.TestCase):
    OUT_FILE = "working/test_out.log"
    ERROR_FILE = "working/test_error.log"
    
    def setUp(self):
        purge_files(self.OUT_FILE, self.ERROR_FILE)
        self.open_files = []
    
    def tearDown(self):
        for f in self.open_files:
            f.close()
        purge_files(self.OUT_FILE, self.ERROR_FILE)
    
    def test_init(self):
        job = Job("python -c pass", self.OUT_FILE, self.ERROR_FILE)
        
        self.assertEqual(job.standard_output, self.OUT_FILE)
        
    def test_repr(self):
        job = Job("python -c pass", self.OUT_FILE, self.ERROR_FILE)

        self.assertEquals(repr(job), "Job('python -c pass')")
        
    def test_teardown(self):
        job = Job("python -c pass", self.OUT_FILE, self.ERROR_FILE)
        out_file = open(job.standard_output, "a", 0)
        #TODO: make the Job class responsible for opening the file.
        self.open_files.append(out_file)
        
        job.standard_output_f = out_file
        job.teardown()
        
        self.assertTrue(out_file.closed, "Output file is closed")
        
class Job_QueueTest(unittest.TestCase):
    def test_add_job(self):
        job = Job("python -c pass")
        queue = Job_Queue()
        
        start_count = queue.num_jobs()
        queue.add_job(job)
        end_count = queue.num_jobs()
        
        self.assertEquals(start_count, 0)
        self.assertEquals(end_count, 1)
        
    def test_get_next_job(self):
        job = Job("python -c pass")
        queue = Job_Queue()
        
        queue.add_job(job)
        next_job = queue.get_next_job()
        
        self.assertIs(next_job, job)
        
    def test_get_next_job_from_empty_queue(self):
        queue = Job_Queue()
        
        # TODO: switch this to RuntimeError or Queue.Empty exception
        with self.assertRaises(Exception) as context:
            queue.get_next_job()
        
        self.assertEquals(context.exception.message, "No jobs in queue")

class WorkerTest(unittest.TestCase):
    OUT_FILE = 'working/test_out.log'
    ERROR_FILE = 'working/test_error.log'
    
    def setUp(self):
        purge_files(self.OUT_FILE, self.ERROR_FILE)
        
    def tearDown(self):        
        purge_files(self.OUT_FILE, self.ERROR_FILE)        
        
    def test_start_job(self):
        command = """python -c "print('Hello, World!')" """
        job = Job(command, 
                  standard_output=self.OUT_FILE, 
                  standard_error=self.ERROR_FILE)
        worker = Worker()
         
        process, full_command = worker.start_job(job)        
        process.wait()
        
        with open(self.OUT_FILE) as f:
            self.assertEquals(f.read(), "Hello, World!\n")
        self.assertEquals(full_command, " " + command)
        
    def test_start_job_no_redirect(self):
        job = Job('python -c pass')
        worker = Worker()
         
        process, _ = worker.start_job(job)        
        process.wait()

    def test_available_for_work(self):
        job = Job('python -c "import time; time.sleep(0.5)"')
        worker = Worker()
        
        worker.start_job(job)
        is_available_immediately = worker.available_for_work()
        while not worker.available_for_work():
            time.sleep(0.05)
            
        self.assertFalse(is_available_immediately)
        
    def test_zero_return(self):
        job = Job('python -c "exit(0)"')
        expected_error = ""
        stderr = StringIO()
        worker = Worker(stderr=stderr)
        
        worker.start_job(job)
        while not worker.available_for_work():
            time.sleep(0.05)
        
        error = stderr.getvalue()
        error_tail = error[-len(expected_error):]
        self.assertEquals(error_tail, expected_error)
        
    def test_nonzero_return(self):
        job = Job('python -c "exit(1)"')
        expected_error = """ - non-zero exit code 1 from command 'python -c "exit(1)"'
"""
        stderr = StringIO()
        worker = Worker(stderr=stderr)
        
        worker.start_job(job)
        while not worker.available_for_work():
            time.sleep(0.05)
        
        error = stderr.getvalue()
        error_tail = error[-len(expected_error):]
        self.assertEquals(error_tail, expected_error)

    def test_clean_terminate(self):
        timeout = time.time() + 1 # 1 second in the future
        job = Job('python -c "import time; time.sleep(30)"')
        stderr = StringIO() # Just to suppress the error message
        worker = Worker(stderr=stderr)
        
        worker.start_job(job)
        worker.clean_terminate()
        
        while not worker.available_for_work():
            if time.time() > timeout:
                self.fail("Timed out waiting for process to end.")
            time.sleep(0.1)
        
    def test_clean_terminate_after_process_ends(self):
        job = Job('python -c pass')
        worker = Worker()
         
        process, _ = worker.start_job(job)
        process.wait()

        worker.clean_terminate()
        
class FactoryTest(unittest.TestCase):
    def test_init(self):
        expected_resources = ["res1", "res1", "res2", "res2", "res2"]
        factory = Factory([("res1", 2), ("res2", 3)])
        
        resources = [worker.resource_allocated for worker in factory.workers]
        self.assertEquals(resources, expected_resources)
    
    def test_hire_worker(self):
        expected_resources = ["res1", "res1", "res2"]
        factory = Factory([("res1", 2)])
        factory.hire_worker("res2")
        
        resources = [worker.resource_allocated for worker in factory.workers]
        self.assertEquals(resources, expected_resources)

    def test_queue_work(self):
        factory = Factory([("", 1)])
        
        process, _ = factory.queue_work('python -c pass')
        
        process.wait()
        self.assertTrue(factory.completely_idle())
        
    def test_supervise(self):
        timeout = time.time() + 1 # 1 second in the future
        stderr = StringIO() # suppress error output
        factory = Factory([("", 1)], stderr=stderr)
        
        slow_process, _ = factory.queue_work(
            'python -c "import time; time.sleep(30)"')
        result2 = factory.queue_work('python -c pass')
        
        slow_process.terminate()
        
        while not factory.completely_idle():
            if time.time() > timeout:
                self.fail("Timed out waiting for factory to finish.")
                
            factory.supervise()
            
            time.sleep(0.1)
        
        self.assertIs(result2, None) # result when queuing to busy factory
    
    def test_wait(self):
        notified_commands = StringIO()
        expected_commands = """\
 python -c "import time; time.sleep(0.5)"
 python -c "x = 1"
 python -c "x = 2"
"""
        def callback(process, command):
            notified_commands.write(command + "\n")
            
        factory = Factory([("", 1)], launch_callback=callback)
        
        factory.queue_work('python -c "import time; time.sleep(0.5)"')
        factory.queue_work('python -c "x = 1"')
        factory.queue_work('python -c "x = 2"')
        
        sleep_seconds = 0.1
        factory.wait(sleep_seconds)
        
        self.assertEquals(notified_commands.getvalue(), expected_commands)
