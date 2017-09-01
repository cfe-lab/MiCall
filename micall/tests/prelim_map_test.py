import os
from random import randrange
from tempfile import NamedTemporaryFile
from unittest import TestCase

from micall.core.prelim_map import check_fastq


class CheckFastqTest(TestCase):
    def create_temp_file(self, suffix, cleanup=True):
        f = NamedTemporaryFile(mode='w', dir='.', suffix=suffix, delete=False)
        if cleanup:
            self.addCleanup(os.unlink, f.name)
        return f

    def test_simple(self):
        f = self.create_temp_file(suffix='.fastq')
        f.close()

        new_name = check_fastq(f.name)

        self.assertEqual(f.name, new_name)

    def test_missing(self):
        f = self.create_temp_file(suffix='.fastq', cleanup=False)
        f.close()
        os.unlink(f.name)

        with self.assertRaisesRegex(SystemExit, 'No FASTQ found at'):
            check_fastq(f.name)

    def test_gzipped(self):
        f = self.create_temp_file(suffix='.fastq.gz')
        f.close()

        new_name = check_fastq(f.name, gzip=True)

        self.assertEqual(f.name, new_name)

    def test_gzipped_link(self):
        f = self.create_temp_file(suffix='.txt')
        expected_content = str(randrange(1000000))
        f.write(expected_content)
        f.close()
        expected_new_name = f.name + '.gz'
        self.addCleanup(os.unlink, expected_new_name)

        new_name = check_fastq(f.name, gzip=True)

        with open(new_name) as f2:
            content = f2.read()

        self.assertEqual(expected_new_name, new_name)
        self.assertEqual(expected_content, content)

    def test(self):
        f = self.create_temp_file(suffix='.txt')
        f.close()
        expected_new_name = f.name + '.gz'
        with open(expected_new_name, 'w'):
            pass
        self.addCleanup(os.unlink, expected_new_name)

        new_name = check_fastq(f.name, gzip=True)

        self.assertEqual(expected_new_name, new_name)
