"""
Unit tests for the generic SortedRing data structure.
"""
import pytest

from micall.utils.sorted_ring import SortedRing


@pytest.fixture
def example_ring() -> SortedRing[int]:
    """
    Provides a pre-populated SortedRing of integers 1 through 5.
    """
    # capacity must accommodate all inserted items
    ring = SortedRing[int](capacity=5)
    for value in [3, 1, 5, 2, 4]:
        ring.insert(value)
    # After insertion, data should be sorted: [1, 2, 3, 4, 5]
    return ring


def test_len_and_iter(example_ring):
    assert len(example_ring) == 5
    # Iteration yields sorted values
    assert list(example_ring) == [1, 2, 3, 4, 5]


def test_contains_and_to_list(example_ring):
    # __contains__ should reflect membership
    assert 3 in example_ring
    assert 0 not in example_ring
    # to_list should copy the internal data
    data = list(example_ring)
    assert data == [1, 2, 3, 4, 5]
    data.append(6)
    # Original ring remains unchanged
    assert len(example_ring) == 5


def test_getitem_and_delitem(example_ring):
    # __getitem__ for single index
    assert example_ring[0] == 1
    assert example_ring[-1] == 5
    # __getitem__ for slice returns a list
    assert example_ring[1:4] == [2, 3, 4]

    # __delitem__ for single index
    # copy to a new ring with sufficient capacity
    ring_copy = SortedRing[int](capacity=5)
    for x in list(example_ring):
        ring_copy.insert(x)
    del ring_copy[0]
    assert list(ring_copy) == [2, 3, 4, 5]

    # __delitem__ for slice
    del ring_copy[1:3]
    # Should remove items at positions 1 and 2 (original 3 and 4)
    assert list(ring_copy) == [2, 5]


def test_remove_and_clear(example_ring):
    ring = example_ring
    # Remove an existing element
    ring.remove(3)
    assert list(ring) == [1, 2, 4, 5]

    # Removing non-existent should raise ValueError
    with pytest.raises(ValueError):
        ring.remove(999)

    # Clear the ring
    ring.clear()
    assert len(ring) == 0
    assert list(ring) == []


def test_add_alias_and_basic_insert():
    # capacity must accommodate added items
    ring = SortedRing[int](capacity=3)
    # add() is alias for insert()
    ring.add(10)
    ring.add(5)
    ring.add(7)
    assert list(ring) == [5, 7, 10]


def test_pop_smallest_and_largest(example_ring):
    ring = example_ring
    smallest = ring.pop_smallest()
    assert smallest == 1
    assert list(ring) == [2, 3, 4, 5]

    largest = ring.pop_largest()
    assert largest == 5
    assert list(ring) == [2, 3, 4]


def test_resize_unlimited_to_limited(example_ring):
    ring = example_ring
    # initial capacity is set by fixture
    assert ring.capacity == 5

    # Resize down to 3: should keep three largest [3,4,5]
    ring.resize(3)
    assert ring.capacity == 3
    assert list(ring) == [3, 4, 5]

    # Further resize to 1: keep [5]
    ring.resize(1)
    assert list(ring) == [5]


def test_resize_increase_capacity(example_ring):
    # Start with capacity 2
    ring = SortedRing[int](capacity=2)
    ring.insert(4)
    ring.insert(1)
    # Data sorted [1,4]
    assert list(ring) == [1, 4]

    # Increase capacity to accommodate more items
    ring.resize(5)
    assert ring.capacity == 5
    # Can insert more without trimming
    for x in [2, 3, 5]:
        ring.insert(x)
    assert list(ring) == [1, 2, 3, 4, 5]


def test_repr_shows_capacity_and_data():
    ring = SortedRing[int](capacity=3)
    for x in [2, 1, 3, 4]:
        ring.insert(x)
    # Only three largest retained: [2,3,4]
    rep = repr(ring)
    assert 'SortedRing' in rep
    assert 'capacity=3' in rep
    # The internal data structure is a SortedList
    assert 'data=SortedList([2, 3, 4])' in rep
