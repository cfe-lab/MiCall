import os
from StringIO import StringIO
from subprocess import CalledProcessError
import time
import unittest

from fifo_scheduler import Factory, Job, Worker
from testfixtures.logcapture import LogCapture

def purge_files(*filenames):
    for filename in filenames:
        try:
            os.remove(filename)
        except OSError:
            pass
    
class JobTest(unittest.TestCase):
    OUT_FILE = "working/test_out.log"
    ERROR_FILE = "working/test_error.log"
    
    def test_init(self):
        job = Job("python -c pass", self.OUT_FILE, self.ERROR_FILE)
        
        self.assertEqual(job.standard_output, self.OUT_FILE)
        
    def test_repr(self):
        job = Job("python -c pass", self.OUT_FILE, self.ERROR_FILE)

        self.assertEquals(repr(job), "Job('python -c pass')")
        
class WorkerTest(unittest.TestCase):
    OUT_FILE = 'working/test_out.log'
    ERROR_FILE = 'working/test_error.log'
    
    def setUp(self):
        purge_files(self.OUT_FILE, self.ERROR_FILE)
        
    def tearDown(self):        
        purge_files(self.OUT_FILE, self.ERROR_FILE)        
        
    def test_run_job(self):
        command = """python -c "print('Hello, World!')" """
        job = Job(command, 
                  standard_output=self.OUT_FILE, 
                  standard_error=self.ERROR_FILE)
        worker = Worker()
        
        worker.run_job(job)        
        
        with open(self.OUT_FILE) as f:
            self.assertEquals(f.read(), "Hello, World!\n")
        
    def test_run_job_no_redirect(self):
        job = Job('python -c pass')
        worker = Worker()
          
        worker.run_job(job)        
         
    def test_callback(self):
        self.actual_command = None
        def callback(command):
            self.actual_command = command
        requested_command = 'python -c pass'
        job = Job(requested_command)
        expected_command = ' ' + requested_command
        worker = Worker(launch_callback=callback)
           
        worker.run_job(job)
        
        self.assertEqual(self.actual_command, expected_command)
         
    def test_callback_with_resource(self):
        self.actual_command = None
        def callback(command):
            self.actual_command = command
        requested_command = 'python -c pass'
        job = Job(requested_command)
        resource = 'MYVAR=23'
        expected_command = resource + ' ' + requested_command
        worker = Worker(resource=resource, launch_callback=callback)
           
        worker.run_job(job)
        
        self.assertEqual(self.actual_command, expected_command)
         
    def test_zero_return(self):
        job = Job('python -c "exit(0)"')
        worker = Worker()
        
        returncode = worker.run_job(job)
        
        self.assertEqual(returncode, 0)
         
    def test_nonzero_return(self):
        job = Job('python -c "exit(1)"')
        expected_error = """Command ' python -c "exit(1)"' returned non-zero exit status 1"""
        worker = Worker()
        
        exception_message = None
        try:
            worker.run_job(job)
            
            self.fail("Should have thrown")
        except CalledProcessError as e:
            exception_message = str(e)
         
        self.assertEquals(exception_message, expected_error)

class FactoryTest(unittest.TestCase):
    def test_init(self):
        expected_resources = ["res1", "res1", "res2", "res2", "res2"]
        factory = Factory([("res1", 2), ("res2", 3)])
        
        resources = [worker.resource_allocated for worker in factory.workers]
        self.assertEquals(resources, expected_resources)
    
    def test_wait(self):
        notified_commands = StringIO()
        expected_commands = """\
 python -c "import time; time.sleep(0.5)"
 python -c "x = 1"
 python -c "x = 2"
"""
        def callback(command):
            notified_commands.write(command + "\n")
            
        factory = Factory([("", 1)], launch_callback=callback)
        
        factory.queue_work('python -c "import time; time.sleep(0.5)"')
        factory.queue_work('python -c "x = 1"')
        factory.queue_work('python -c "x = 2"')
        
        factory.wait()
        
        self.assertEquals(notified_commands.getvalue(), expected_commands)
    
    def test_wait_nonzero(self):
        expected_error = """Command ' python -c "exit(1)"' returned non-zero exit status 1"""
        factory = Factory([("", 1)])
        
        factory.queue_work('python -c "exit(1)"')

        exception_message = None
        try:
            factory.wait()
            
            self.fail("Should have thrown")
        except CalledProcessError as e:
            exception_message = str(e)
         
        self.assertEquals(exception_message, expected_error)
    
    def test_wait_multiple_nonzero(self):
        expected_error = """Command ' python -c "exit(2)"' returned non-zero exit status 2"""
        factory = Factory([("", 1)])
        
        factory.queue_work('python -c "exit(1)"')
        factory.queue_work('python -c "exit(2)"')

        exception_message = None
        try:
            with LogCapture() as log:
                factory.wait()
            
            self.fail("Should have thrown")
        except CalledProcessError as e:
            exception_message = str(e)
         
        self.assertEquals(exception_message, expected_error)
        log.check(('fifo_scheduler', 
                   'ERROR', 
                   """Command ' python -c "exit(1)"' returned non-zero exit status 1"""))
    
    def test_parallel(self):
        notified_commands = []
        expected_commands = """\
 python -c "x=0; import time; time.sleep(0.5)"
 python -c "x=1; import time; time.sleep(0.5)"
 python -c "x=2; import time; time.sleep(0.5)"
 python -c "x=3; import time; time.sleep(0.5)"\
"""

        def callback(command):
            notified_commands.append(command)
            
        factory = Factory([("", 2)], launch_callback=callback)
        start = time.time()
        
        for i in range(4):
            script = 'python -c "x=%d; import time; time.sleep(0.5)"' % i
            factory.queue_work(script)
        
        factory.wait()

        duration = time.time() - start
        notified_commands.sort()
        self.assertEquals('\n'.join(notified_commands), expected_commands)
        self.assertGreater(duration, 0.5)
        self.assertLess(duration, 1.5)

    def test_result(self):
        factory = Factory([("", 1)])
        
        async_result = factory.queue_work('python -c "exit(0)"')
        
        returncode = async_result.get()
        
        self.assertEquals(returncode, 0)

    def test_result_nonzero(self):
        expected_error = """Command ' python -c "exit(1)"' returned non-zero exit status 1"""
        factory = Factory([("", 1)])
        
        async_result = factory.queue_work('python -c "exit(1)"')

        exception_message = None
        try:
            async_result.get()
            
            self.fail("Should have thrown")
        except CalledProcessError as e:
            exception_message = str(e)
         
        self.assertEquals(exception_message, expected_error)

    def test_failing_callback(self):
        command = 'python -c pass'
        expected_message = (
            "Failed callback for ' python -c pass'. KeyError: 'missing_key'")

        def failing_callback(command):
            a = {}
            a['x'] = a['missing_key']
            
        factory = Factory([("", 1)], launch_callback=failing_callback)
        
        factory.queue_work(command)
        
        exception_message = None
        try:
            factory.wait()
            
            self.fail("Should have thrown")
        except RuntimeError as e:
            exception_message = str(e)
        
        self.assertEquals(exception_message, expected_message)
