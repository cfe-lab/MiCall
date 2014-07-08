import os
from StringIO import StringIO
from subprocess import CalledProcessError
import time
import unittest

from fifo_scheduler import Factory, Job, Worker
from testfixtures.logcapture import LogCapture
import stat

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
        
        self.assertEqual(job.stdout, self.OUT_FILE)
        
    def test_repr(self):
        job = Job("python -c pass", self.OUT_FILE, self.ERROR_FILE)

        self.assertEquals(repr(job), "Job('python -c pass')")
        
class WorkerTest(unittest.TestCase):
    IN_FILE = 'working/test_in.log'
    OUT_FILE = 'working/test_out.log'
    OUT_FILE2 = 'working/test_out2.log'
    ERROR_FILE = 'working/test_error.log'
    HELPER_FILE = 'working/helper.py'
    SCRIPT_FILE = 'working/script.py'
    ALL_FILES = [IN_FILE,
                 OUT_FILE,
                 OUT_FILE2,
                 ERROR_FILE,
                 HELPER_FILE,
                 SCRIPT_FILE,
                 'working.txt']
     
    def setUp(self):
        purge_files(*self.ALL_FILES)
         
    def tearDown(self):        
        purge_files(*self.ALL_FILES)
    
    def assertFileContents(self, filename, expectedContents):
        with open(filename, "r") as f:
            self.assertEquals(f.read(), expectedContents)
             
    def prepareFileContents(self, filename, contents, mode=None):
        with open(filename, "w") as f:
            f.write(contents)
        if mode is not None:
            os.chmod(filename, mode)

    def test_run_job(self):
        command = """python -c "print('Hello, World!')" """
        job = Job(command, 
                  stdout=self.OUT_FILE, 
                  stderr=self.ERROR_FILE)
        worker = Worker()
        
        worker.run_job(job)
        
        self.assertFileContents(self.OUT_FILE, "Hello, World!\n")
        
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
        
    def test_resource(self):
        job = Job("""python -c "import os; print os.environ['TEST_VALUE']" """,
                  stdout=self.OUT_FILE)
        expected_text = 'ExpectedText\n'
        worker = Worker("TEST_VALUE=ExpectedText")
        
        worker.run_job(job)
        
        self.assertFileContents(self.OUT_FILE, expected_text)
         
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
        
    def test_script_file(self):
        """ Define a job as a script file to execute. """
        
        self.prepareFileContents(
            self.SCRIPT_FILE,
            """\
#!/usr/bin/env python

print 'expected text'
""",
            mode=stat.S_IRWXU)
        expected_output = 'expected text\n'
        worker = Worker()
             
        worker.run_job(Job(script=self.SCRIPT_FILE, 
                           stdout=self.OUT_FILE))
        
        self.assertFileContents(self.OUT_FILE, expected_output)
        
    def test_script_file_with_resource(self):
        self.prepareFileContents(
            self.SCRIPT_FILE,
            """\
#!/usr/bin/env python

import os

print os.environ['TEST_VALUE']
""",
            mode=stat.S_IRWXU)
        expected_output = 'ExpectedText\n'
        
        worker = Worker("TEST_VALUE=ExpectedText")
        
        worker.run_job(Job(script=self.SCRIPT_FILE, 
                           stdout=self.OUT_FILE))
        
        self.assertFileContents(self.OUT_FILE, expected_output)

    def test_local_data_file(self):
        """ Scripts should be allowed to create any temporary files they want
        in their current directory. Parallel executions should each get a
        separate directory."""
         
        self.prepareFileContents(
            self.SCRIPT_FILE,
            """\
#!/usr/bin/env python

with open('working.txt', 'a') as f:
    f.write('expected text')
  
with open('working.txt', 'r') as f:
    print f.read()
""",
            mode=stat.S_IRWXU)
        expected_output = 'expected text\n'
              
        worker = Worker()
        
        worker.run_job(Job(script=self.SCRIPT_FILE, 
                           stdout=self.OUT_FILE))
        worker.run_job(Job(script=self.SCRIPT_FILE, 
                           stdout=self.OUT_FILE2))
        
        self.assertFileContents(self.OUT_FILE, expected_output)
        self.assertFileContents(self.OUT_FILE2, expected_output)        

    def test_helper_files(self):
        """ Scripts may need helper files that contain static data or library
        code."""
         
        self.prepareFileContents(
            self.HELPER_FILE,
            """\
def greet(name):
    return "Hello, %s" % name
""")
        self.prepareFileContents(
            self.SCRIPT_FILE,
            """\
#!/usr/bin/env python

import helper

print helper.greet('Bob')
""",
            mode=stat.S_IRWXU)
        expected_output = 'Hello, Bob\n'
              
        worker = Worker()
        
        worker.run_job(Job(script=self.SCRIPT_FILE,
                           helpers=[self.HELPER_FILE],
                           stdout=self.OUT_FILE))
        
        self.assertFileContents(self.OUT_FILE, expected_output)

    def test_script_args(self):
        self.prepareFileContents(
            self.SCRIPT_FILE,
            """\
#!/usr/bin/env python

import sys

print sys.argv[1]
""",
            mode=stat.S_IRWXU)
        expected_output = 'expectedText\n'
        
        worker = Worker()
        
        worker.run_job(Job(script=self.SCRIPT_FILE,
                           args=('expectedText', ),
                           stdout=self.OUT_FILE))
        
        self.assertFileContents(self.OUT_FILE, expected_output)

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
    
    def test_queue_job(self):
        notified_commands = StringIO()
        expected_commands = " python -c pass\n"
        def callback(command):
            notified_commands.write(command + "\n")
            
        factory = Factory([("", 1)], launch_callback=callback)
        
        factory.queue_job(Job('python -c pass'))
        
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
