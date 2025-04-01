
from typing import TypeVar, Generic, Iterable, Iterator, Collection, \
    Protocol, runtime_checkable
from sortedcontainers import SortedList


@runtime_checkable
class Comparable(Protocol):
    def __lt__(self, other) -> bool: ...
    def __gt__(self, other) -> bool: ...


T = TypeVar('T', bound=Comparable)


class SortedRing(Generic[T], Collection[T], Iterable[T]):
    def __init__(self, capacity: int):
        self.data: SortedList[T] = SortedList()
        self.size: int = 0
        self.capacity: int = capacity
        self._smallest_item: T = None  # type: ignore
        self._largest_item: T = None  # type: ignore

    def smallest_item(self) -> T:
        if self.size == 0:
            raise IndexError("Sorted ring is empty.")

        return self._smallest_item

    def largest_item(self) -> T:
        if self.size == 0:
            raise IndexError("Sorted ring is empty.")

        return self._largest_item

    def __contains__(self, item: object) -> bool:
        return item in self.data

    def __len__(self) -> int:
        return self.size

    def __iter__(self) -> Iterator[T]:
        return iter(self.data)

    def is_full(self) -> bool:
        return self.size >= self.capacity

    def add(self, item: T) -> bool:
        if self.capacity <= 0:
            return False

        if self.is_full():
            if item < self._smallest_item:
                return False

            del self.data[0]
            self.size -= 1

        if self.size > 0:
            if item < self._smallest_item:
                self._smallest_item = item

            if item > self._largest_item:
                self._largest_item = item

        else:
            self._smallest_item = item
            self._largest_item = item

        self.size += 1
        self.data.add(item)
        return True
