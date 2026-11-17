"""The stitcher plot must handle large shared-ancestry graphs in linear time."""

import unittest

from micall.core.plot_contigs import close_alive_set


class AncestryClosureTests(unittest.TestCase):
    def test_branching_diamond(self):
        parents = {4: [3, 2], 3: [1], 2: [1], 1: []}
        alive = {4}
        close_alive_set(alive, parents)
        self.assertEqual(alive, {1, 2, 3, 4})

    def test_cycles(self):
        alive = {1}
        close_alive_set(alive, {1: [2], 2: [1]})
        self.assertEqual(alive, {1, 2})

    def test_shared_ancestors_do_not_explode(self):
        parents = {1: []}
        for level in range(1, 250):
            left, right = 2 * level, 2 * level + 1
            prev = [2 * (level - 1), 2 * (level - 1) + 1] if level > 1 else [1]
            parents[left] = prev
            parents[right] = prev
        alive = {498, 499}
        close_alive_set(alive, parents)
        self.assertEqual(alive, set(parents))

    def test_long_chain_never_hits_recursion_limit(self):
        parents = {i: [i - 1] for i in range(1, 10001)}
        alive = {10000}
        close_alive_set(alive, parents)
        self.assertEqual(len(alive), 10001)
        self.assertIn(0, alive)

    def test_empty_source_is_empty(self):
        alive = set()
        close_alive_set(alive, {2: [3]})
        self.assertFalse(alive)


if __name__ == "__main__":
    unittest.main()
