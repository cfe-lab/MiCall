from collections import Counter
from csv import DictWriter, DictReader
import os
from tempfile import NamedTemporaryFile


def read_cache_lines(cache_file_name: str):
    """ Yield one line at a time from a cache file.

    The point of this generator function is to only open the file when you
    read each line. Otherwise, a large set of counts with many cache files will
    open too many files at once.
    """
    pos = 0
    while True:
        with open(cache_file_name) as f:
            f.seek(pos)
            line = f.readline()
            if not line:
                break
            pos = f.tell()
        yield line


class BigCounter:
    def __init__(self, file_prefix, max_size=5000):
        self.file_prefix = os.path.abspath(file_prefix)
        self.max_size = max_size
        self.cache_files = []
        self.active_counts = Counter()

    def __enter__(self):
        return self
        
    def __exit__(self, exc_type, exc_value, traceback):
        self.clear()
    
    def __getitem__(self, key):
        return self.active_counts[key]

    def __setitem__(self, key, value):
        if not isinstance(key, str):
            raise TypeError('Key was not a string: {!r}'.format(key))

        self.active_counts[key] = value
        if len(self.active_counts) > self.max_size:
            self._write_cache()
            
    def _write_cache(self):
        with NamedTemporaryFile(mode='w+',
                                dir=os.path.dirname(self.file_prefix),
                                prefix=os.path.basename(self.file_prefix),
                                suffix='.csv',
                                delete=False) as cache:
            keys = sorted(self.active_counts)
            writer = DictWriter(cache, ['key', 'count'])
            for key in keys:
                writer.writerow(dict(key=key, count=self.active_counts[key]))
            self.active_counts.clear()
            self.cache_files.append(cache.name)

    def items(self):
        item_count = 0
        active_keys = sorted(self.active_counts)
        readers = [(dict(key=key, count=self.active_counts[key])
                    for key in active_keys)]
        for cache_file_name in self.cache_files:
            readers.append(DictReader(read_cache_lines(cache_file_name),
                                      ['key', 'count']))
        next_rows = [next(reader) for reader in readers]
        while any(next_rows):
            min_key = min(row['key']
                          for row in next_rows
                          if row is not None)
            count = 0
            for i, row in enumerate(next_rows):
                if row is not None and row['key'] == min_key:
                    count += int(row['count'])
                    try:
                        next_rows[i] = next(readers[i])
                    except StopIteration:
                        # noinspection PyTypeChecker
                        next_rows[i] = None
            yield min_key, count
            item_count += 1

    def clear(self):
        for cache_file in self.cache_files:
            os.remove(cache_file)
        self.cache_files.clear()
        self.active_counts.clear()
