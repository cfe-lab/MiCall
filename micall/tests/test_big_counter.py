import os
from unittest import TestCase

from micall.utils.big_counter import BigCounter

FILE_PREFIX = os.path.join(os.path.dirname(__file__), 'big_counter_test_data')


class BigCounterTest(TestCase):
    def test_no_writing(self):
        expected_items = [('a', 2), ('b', 1)]

        with BigCounter(FILE_PREFIX, max_size=2) as counter:
            counter['a'] += 1
            counter['b'] += 1
            counter['a'] += 1
            final_size = len(counter.active_counts)
            items = sorted(counter.items())

        self.assertEqual(expected_items, items)
        self.assertEqual(2, final_size)

    def test_writing(self):
        expected_items = [('a', 2), ('b', 1), ('c', 1)]

        with BigCounter(FILE_PREFIX, max_size=2) as counter:
            counter['a'] += 1
            counter['b'] += 1
            counter['c'] += 1
            counter['a'] += 1
            final_size = len(counter.active_counts)
            items = sorted(counter.items())

        self.assertEqual(expected_items, items)
        self.assertEqual(1, final_size)

    def test_key_not_string(self):
        with BigCounter(FILE_PREFIX) as counter:
            with self.assertRaisesRegex(TypeError, 'Key was not a string: 23'):
                counter[23] = 5

    def test_empty(self):
        with BigCounter(FILE_PREFIX) as counter:
            items = list(counter.items())

        self.assertEqual(items, [])
