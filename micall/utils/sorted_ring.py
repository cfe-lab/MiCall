
"""
Generic sorted ring data structure.

SortedRing maintains a collection of items in sorted order, allowing
efficient insertion and removal.
"""


from typing import TypeVar, Generic, Iterator, Sequence, Any
from collections.abc import MutableSet
from sortedcontainers import SortedList


T = TypeVar('T')


class SortedRing(Sequence, MutableSet, Generic[T]):
    """
    A minimal generic sorted ring data structure with fixed capacity.
    Inserting beyond capacity automatically removes the smallest items
    to maintain the maximum size.
    """

    def __init__(self, capacity: int) -> None:
        if capacity < 1:
            raise ValueError("SortedRing capacity must be greater than 0.", capacity)

        # capacity: maximum number of items allowed
        self._capacity: int = capacity
        # use SortedList for efficient sorted operations
        self._data: SortedList[T] = SortedList()
        # track size locally for O(1) len()
        self._size: int = 0

    def __len__(self) -> int:
        # O(1) length via local counter
        return self._size

    def __iter__(self) -> Iterator[T]:
        return iter(self._data)

    def insert(self, item: T) -> bool:
        """Insert an item, keep sorted via SortedList, and enforce capacity."""

        if self._size >= self._capacity:
            if item <= self._data[0]:
                return False

            # enforce capacity: remove the smallest item if over capacity.
            self._data.pop(0)
            self._size -= 1

        self._data.add(item)
        self._size += 1

        return True

    def remove(self, item: T) -> None:
        """Remove an item and update size."""
        self._data.remove(item)
        self._size -= 1

    def clear(self) -> None:
        """Clear all items and reset size."""
        self._data.clear()
        self._size = 0

    def add(self, item: T) -> None:
        """Alias to insert, for compatibility with SortedList."""
        self.insert(item)

    def __getitem__(self, index: Any) -> T:  # type: ignore
        """Get item or slice from the sorted data."""
        return self._data[index]

    def __delitem__(self, index: Any) -> None:  # type: ignore
        """Delete item or slice from the sorted data and update size."""
        del self._data[index]
        # recalculate size
        self._size = len(self._data)

    def resize(self, new_capacity: int) -> None:
        """
        Set a new capacity and remove smallest items if the data exceeds it.

        new_capacity: maximum number of items to keep.
        """

        if new_capacity < 1:
            raise ValueError("SortedRing capacity must be greater than 0.", new_capacity)

        self._capacity = new_capacity
        # remove smallest items until within new capacity
        while self._size > new_capacity:
            self._data.pop(0)
            self._size -= 1

    @property
    def capacity(self) -> int:
        """Return the current capacity."""
        return self._capacity

    def pop_smallest(self) -> T:
        """Remove and return the smallest item, updating size."""
        item = self._data.pop(0)
        self._size -= 1
        return item

    def pop_largest(self) -> T:
        """Remove and return the largest item, updating size."""
        item = self._data.pop()
        self._size -= 1
        return item

    def __contains__(self, item: object) -> bool:
        """Return True if item is in the sorted ring."""
        return item in self._data

    def __repr__(self) -> str:
        """Return string representation with data and capacity."""
        return (
            f"{self.__class__.__name__}(data={self._data!r}, "
            f"capacity={self._capacity})"
        )

    def discard(self, item: T) -> None:
        """Remove item if present; do nothing otherwise."""
        try:
            self.remove(item)
        except ValueError:
            pass
