from unittest import TestCase
from unittest.mock import patch

from micall.utils.basespace_upload import upload


class BasespaceUploadTest(TestCase):
    def setUp(self):
        self.expected_file_groups = []
        patcher = patch('micall.utils.basespace_upload.check_call',
                        self.mock_call)
        patcher.start()
        self.addCleanup(patcher.stop)

    def tearDown(self):
        self.assertEqual([], self.expected_file_groups)

    def mock_call(self, shell_args):
        file_group = shell_args[5:]
        expected_file_group = self.expected_file_groups.pop(0)
        self.assertEqual(expected_file_group, file_group)

    def test_singles(self):
        filenames = ['apple.fastq.gz', 'banana.fastq.gz']
        self.expected_file_groups = [['apple.fastq.gz'], ['banana.fastq.gz']]

        upload('project name', filenames)

    def test_groups(self):
        filenames = ['apple_R1.fastq.gz',
                     'apple_R2.fastq.gz',
                     'banana_R1.fastq.gz',
                     'banana_R2.fastq.gz']
        self.expected_file_groups = [['apple_R1.fastq.gz',
                                      'apple_R2.fastq.gz'],
                                     ['banana_R1.fastq.gz',
                                      'banana_R2.fastq.gz']]

        upload('project name', filenames)

    def test_sorted(self):
        filenames = ['foo/apple_R1.fastq.gz',
                     'foo/apple_R2.fastq.gz',
                     'bar/banana_R2.fastq.gz',
                     'bar/banana_R1.fastq.gz']
        self.expected_file_groups = [['bar/banana_R1.fastq.gz',
                                      'bar/banana_R2.fastq.gz'],
                                     ['foo/apple_R1.fastq.gz',
                                      'foo/apple_R2.fastq.gz']]

        upload('project name', filenames)
