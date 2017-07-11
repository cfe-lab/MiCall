from collections import Counter
from csv import DictWriter, DictReader
from tempfile import TemporaryFile


class BigCounter:
    def __init__(self, file_prefix, max_size=5000):
        self.file_prefix = file_prefix
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
        cache = TemporaryFile(mode='w+', prefix=self.file_prefix, suffix='.csv')
        self.cache_files.append(cache)
        keys = sorted(self.active_counts)
        writer = DictWriter(cache, ['key', 'count'])
        writer.writeheader()
        for key in keys:
            writer.writerow(dict(key=key, count=self.active_counts[key]))
        self.active_counts.clear()

    def items(self):
        item_count = 0
        active_keys = sorted(self.active_counts)
        readers = [(dict(key=key, count=self.active_counts[key])
                    for key in active_keys)]
        for cache in self.cache_files:
            cache.seek(0)
            # noinspection PyTypeChecker
            readers.append(DictReader(cache))
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
                        next_rows[i] = None
            yield min_key, count
            item_count += 1

    def clear(self):
        for cache in self.cache_files:
            cache.close()
        self.cache_files.clear()
        self.active_counts.clear()
