# MiCall Watcher Modules: Retry System Report

## Executive Summary

The MiCall watcher module implements a retry architecture that handles both
network (Kive API) and disk operations with exponential backoff. The retry system is
distributed across multiple modules (`kive_watcher.py`, `disk_operations.py`, and
`update_qai.py`), each with specific retry strategies tailored to their operational domains.
The system uses consistent exponential backoff parameters, but each module provides
domain-specific retry logic with different maximum retry limits.

---

## 1. Retry Mechanisms Overview

### 1.1 Core Retry Patterns

The system implements four primary retry patterns:

1. **Explicit Loop with Count**: Manual retry loops using `count()` iterator (`kive_watcher.py`, `update_qai.py`)
2. **Decorator Pattern**: Decorators wrapping functions with automatic retry (`disk_operations.py`)
3. **Context Manager Pattern**: File operations wrapped in a retry context manager (`disk_operations.py`)
4. **Helper Function Pattern**: Generic retry operation wrapper (`update_qai.py`)

### 1.2 Unified Retry Configuration

All modules share the same exponential backoff parameters:

- **`MINIMUM_RETRY_WAIT`**: 5 seconds
- **`MAXIMUM_RETRY_WAIT`**: 1 day (86,400 seconds)
- **Exponential Backoff Formula**: `min_seconds * (2 ** (attempt_count - 1))`
  - Attempt 1: 5 seconds
  - Attempt 2: 10 seconds
  - Attempt 3: 20 seconds
  - Attempt 4: 40 seconds
  - And so on until hitting the 1-day maximum

---

## 2. Module-Specific Retry Implementations

### 2.1 kive_watcher.py

#### Global-Level Retry Logic

**`find_samples()` function**

- **Purpose**: Scan for raw data folders and sample groups to process
- **Retry Mechanism**: Infinite loop with exception catching
- **Trigger**: Any unhandled exception during `scan_samples()` or QAI event sending
- **Max Attempts**: Unlimited (infinite retry, controlled by `retry` parameter)
- **Backoff**: `wait_for_retry(attempt_count, start_time)`
- **State Tracking**: `attempt_count` and `start_time` (datetime)
- **Reset Condition**: `attempt_count` and `start_time` reset to 0/None after a successful iteration
- **Controlled by**: `retry` parameter (bool) — if False, raises immediately

**`find_sample_groups()` function**

- **Purpose**: Find and parse sample groups from FASTQ files
- **Retry Mechanism**: Loop using `count(1)`
- **Trigger**: `IOError` from file access (FASTQ file reading/sample sheet parsing)
- **Max Attempts**: Unlimited
- **Backoff**: Exponential via `wait_for_retry()`
- **Exceptions Caught**: `IOError` (retried); all other `Exception` (logged, recorded in error file)
- **State Tracking**: `attempt_count` and `start_time`

**`send_event()` function**

- **Purpose**: Put folder events onto the sample processing queue with timeout
- **Retry Mechanism**: While loop with timeout
- **Trigger**: `Queue.Full` exception
- **Max Attempts**: Timed — continues until `next_scan` datetime passes
- **Backoff**: Implicit in timeout-based loop

**`send_qai_event()` function**

- **Purpose**: Send QAI upload events to the queue
- **Retry Mechanism**: Identical to `send_event()`
- **Timing**: Retries until `next_scan` datetime

#### KiveWatcher Class Methods

**`add_sample_group()` method**

- **Purpose**: Add a sample group to a folder watcher for processing
- **Retry Mechanism**: Loop using `count(1)`
- **Trigger**: Any exception if `self.retry` is True
- **Max Attempts**: Unlimited
- **Backoff**: Exponential via `wait_for_retry()`
- **State Tracking**: `attempt_count` and `start_time`
- **Operations Retried**: Session checking, creating FolderWatcher, uploading quality datasets, opening and reading FASTQ files

**`poll_runs()` method**

- **Purpose**: Poll all active Kive runs for status updates
- **Retry Mechanism**: Loop using `count(1)`
- **Trigger**: Any exception if `self.retry` is True
- **Max Attempts**: Unlimited
- **Backoff**: Exponential via `wait_for_retry()`
- **State Tracking**: `attempt_count` and `start_time`

**`kive_retry()` method**

- **Purpose**: Add a single retry to Kive API calls (handles session refresh)
- **Retry Mechanism**: Try-catch with single retry
- **Trigger**: `KiveClientException`
- **Max Attempts**: 2 (try once; refreshes session on failure; tries again)
- **Logic**: On first failure, logs in again then retries
- **No Backoff**: Immediate retry on second attempt

**`upload_filter_quality()` method**

- **Purpose**: Upload filter quality data for error metrics
- **Retry Mechanism**: Loop using `count(1)`
- **Trigger**: `OSError` when reading InterOp error metrics
- **Max Attempts**: Unlimited
- **Backoff**: Exponential via `wait_for_retry()`
- **Other Exceptions**: Any other exception is logged but not retried

#### HTTPAdapter Configuration

**`open_kive()` function**:

```python
session.mount('https://', HTTPAdapter(max_retries=20))
```

Built-in urllib3 retry mechanism with max 20 retries for HTTPS connections, independent of application-level retry logic.

#### Logging Strategy

The `wait_for_retry()` function implements smart logging:

- **Before 1 hour**: Logs at INFO level
- **After 1 hour**: Logs at WARNING level with `exc_info=True`

This helps distinguish transient issues from prolonged failures.

---

### 2.2 disk_operations.py

#### Decorator-Based Retry Pattern

**`@disk_retry()` decorator**

The decorator wraps disk operations with network-style retry logic:

```python
def disk_retry(operation_name="disk operation"):
    def decorator(func):
        def wrapper(*args, **kwargs):
            attempt_count = 0
            max_attempts = 15
            start_time = None

            while attempt_count < max_attempts:
                attempt_count += 1
                try:
                    return func(*args, **kwargs)
                except (OSError, IOError) as ex:
                    if attempt_count >= max_attempts:
                        logger.error(...)
                        raise

                    if start_time is None:
                        start_time = datetime.now()

                    wait_for_retry(attempt_count, operation_name, start_time)
                except:
                    # Non-disk errors should not be retried
                    raise
```

**Key Characteristics:**

- **Max Retries**: 15 attempts (hard limit)
- **Exceptions Caught**: `OSError` and `IOError` only
- **Other Exceptions**: Never retried; immediately raised
- **Backoff**: Exponential via `wait_for_retry()`
- **State Tracking**: `attempt_count` and `start_time`

#### Decorated Functions

| Function | Operation |
|----------|-----------|
| `mkdir_p()` | Create directory |
| `rmtree()` | Remove directory tree |
| `move()` | Move file/directory |
| `copy_file()` | Copy file |
| `write_text()` | Write text to file |
| `touch()` | Create/touch file |
| `unlink()` | Remove file |
| `rmdir()` | Remove empty directory |
| `rename()` | Rename file/directory |
| `remove_empty_directory()` | Safe empty directory removal |
| `copy_fileobj()` | Copy file object |

#### Logging Strategy

The `wait_for_retry()` function in disk_operations.py:

- **Before 1 hour**: Logs at INFO level
- **After 1 hour**: Logs at ERROR level
- Includes `exc_info=True` for full traceback

#### Context Manager Pattern

**`disk_file_operation` class**:

```python
class disk_file_operation:
    def __init__(self, path: Path, mode: str, operation_name=None):
        ...

    def __enter__(self):
        @disk_retry(self.operation_name)
        def open_file():
            return self.path.open(self.mode)

        self.file_handle = open_file()
        return self.file_handle
```

Used for file reading/writing in `kive_watcher.py`. Retries the file open operation with the same logic as other disk operations, and handles automatic file closing in `__exit__`.

---

### 2.3 update_qai.py

#### Generic Retry Wrapper Function

**`retry_operation()` function**

A generic retry mechanism for both network and disk operations:

```python
def retry_operation(operation, operation_name: str):
    attempt_count = 0
    start_time = None
    last_exception = None

    while attempt_count < MAX_RETRY_ATTEMPTS:
        try:
            return operation()
        except Exception as ex:
            attempt_count += 1
            last_exception = ex

            if start_time is None:
                start_time = datetime.now()

            if attempt_count >= MAX_RETRY_ATTEMPTS:
                logger.error(f'{operation_name} failed after {MAX_RETRY_ATTEMPTS}...')
                break

            logger.warning(f'{operation_name} failed (attempt {attempt_count}/{MAX_RETRY_ATTEMPTS})', ...)
            wait_for_retry(attempt_count, start_time)

    raise RuntimeError(...) from last_exception
```

**Key Characteristics:**

- **Max Retries**: `MAX_RETRY_ATTEMPTS = 100`
- **Exceptions Caught**: All exceptions (generic `Exception`)
- **Backoff**: Exponential via `wait_for_retry()` imported from `kive_watcher`
- **State Tracking**: `attempt_count`, `start_time`, `last_exception`
- **Final State**: Raises `RuntimeError` after all attempts are exhausted

#### Specialized Retry Wrappers

**`retry_qai_login()` function**: Wraps QAI session login using `retry_operation()`.

**`retry_find_run()` function**: Wraps QAI run lookup using `retry_operation()`. Also validates the result count (must be exactly 1 run) and raises `RuntimeError` if find_run returns 0 or multiple results.

#### Integration Points in Upload Processing

**`process_folder()` function**:

- Reads the sample sheet via `retry_operation(read_sample_sheet_with_retry, ...)`
- Logs into QAI via `retry_qai_login(session, qai_server, ...)`
- Finds the run via `retry_find_run(session, experiment_name)`
- Processes remapped data via `retry_operation(process_remapped, ...)`

**`process_remapped()` function**:

- Reads the sample sheet via `retry_operation(read_sample_sheet_with_retry, ...)`
- Reads the conseqs file via `retry_operation(read_conseqs_with_retry, ...)`
- Uploads to QAI via `retry_operation(upload_with_retry, ...)`

**`upload_loop()` function**:

- Checks coverage file via `retry_operation(check_coverage_file, ...)`
- Writes the `done_qai_upload` flag via `retry_operation(write_done_flag, ...)`
- Ensures the done flag is written even if `coverage_scores.csv` is missing

#### Exception Handling Strategy

- **Catches**: All exceptions (network, disk, QAI API failures)
- **Retries**: Everything up to `MAX_RETRY_ATTEMPTS` (100 times)
- **Logs**: WARNING for each failed attempt; ERROR when exhausted
- **Raises**: `RuntimeError` on final failure

---

## 3. Exception Handling Matrix

| Module | Exception Type | Caught | Retried | Max Attempts | Backoff |
|--------|----------------|--------|---------|--------------|---------|
| **kive_watcher** | `IOError` (find_sample_groups) | ✓ | ✓ | Unlimited | Exponential |
| **kive_watcher** | `KiveClientException` (kive_retry) | ✓ | ✓ | 2 | None |
| **kive_watcher** | `OSError` (upload_filter_quality) | ✓ | ✓ | Unlimited | Exponential |
| **kive_watcher** | All other exceptions | ✓ | ✗ | N/A | N/A |
| **disk_operations** | `OSError` | ✓ | ✓ | 15 | Exponential |
| **disk_operations** | `IOError` | ✓ | ✓ | 15 | Exponential |
| **disk_operations** | All other exceptions | ✓ | ✗ | N/A | N/A |
| **update_qai** | All exceptions | ✓ | ✓ | 100 | Exponential |

---

## 4. Retry Limits and Backoff Strategies

### 4.1 Maximum Retry Limits

| Module | Function/Context | Max Attempts | Limit Type |
|--------|-----------------|--------------|------------|
| kive_watcher | `find_samples` | Unlimited | Controlled by parameter |
| kive_watcher | `find_sample_groups` | Unlimited | No hard cap |
| kive_watcher | `add_sample_group` | Unlimited | No hard cap |
| kive_watcher | `poll_runs` | Unlimited | No hard cap |
| kive_watcher | `kive_retry` | 2 | Hard limit + session refresh |
| kive_watcher | `upload_filter_quality` | Unlimited | No hard cap |
| disk_operations | All `@disk_retry` decorated | 15 | Hard limit |
| update_qai | `retry_operation` | 100 | Hard limit (`MAX_RETRY_ATTEMPTS`) |
| update_qai | `retry_qai_login` | 100 | Wrapped by `retry_operation` |

### 4.2 Exponential Backoff Timeline

**Formula**: `delay = min(min_wait * 2^(attempt-1), max_wait)`

```
Attempt 1 failure → Wait   5s
Attempt 2 failure → Wait  10s
Attempt 3 failure → Wait  20s
Attempt 4 failure → Wait  40s
Attempt 5 failure → Wait  80s
Attempt 6 failure → Wait 160s  (2m 40s)
Attempt 7 failure → Wait 320s  (5m 20s)
Attempt 8 failure → Wait 640s  (10m 40s)
Attempt 9 failure → Wait 1280s (21m 20s)
Attempt 10 failure → Wait 2560s (42m 40s)
Attempt 11 failure → Wait 5120s (1h 25m 20s)
Attempt 12+ failure → Wait 86400s (1 day, max reached)
```

### 4.3 Time-Sensitive Logging

The `wait_for_retry()` function switches log level based on elapsed time:

- Within the first hour of failures: INFO
- After 1 hour of continuous failures: WARNING or ERROR (depending on the module)

This helps distinguish transient issues from prolonged outages.

---

## 5. Retry State Tracking

### 5.1 State Variables

| Module | Variable | Type | Purpose |
|--------|----------|------|---------|
| kive_watcher | `attempt_count` | `int` | Current retry attempt number |
| kive_watcher | `start_time` | `datetime` | When the first failure occurred |
| disk_operations | `attempt_count` | `int` | Current retry attempt number |
| disk_operations | `start_time` | `datetime` | When the first failure occurred |
| update_qai | `attempt_count` | `int` | Current retry attempt number |
| update_qai | `start_time` | `datetime` | When the first failure occurred |
| update_qai | `last_exception` | `Exception` | Stores final exception for error reporting |

### 5.2 State Lifecycle

**Initialization**:

```python
attempt_count = 0
start_time = None
```

**On First Failure**:

```python
attempt_count += 1
if start_time is None:
    start_time = datetime.now()
```

**Reset Condition** (kive_watcher scheduler loops only):

```python
# After a fully successful iteration
attempt_count = 0
start_time = None
```

**Use in Backoff Calculation**:

```python
delay = calculate_retry_wait(MINIMUM_RETRY_WAIT, MAXIMUM_RETRY_WAIT, attempt_count)
elapsed = datetime.now() - start_time  # Total time since first failure
```

### 5.3 Cumulative Wait Time Utility

`disk_operations.py` defines `calculate_cumulative_wait_time()`, which computes the total
elapsed wait time across all previous attempts. This function is currently defined but not
actively used, suggesting infrastructure for potential future timeout management.

---

## 6. Cross-Module Retry Interactions

### 6.1 kive_watcher.py → disk_operations.py

```python
# Inside add_sample_group()
disk_operations.rmtree(results_path, ignore_errors=True)   # @disk_retry (max 15)
disk_operations.unlink(results_zip, missing_ok=True)        # @disk_retry (max 15)
```

`add_sample_group()` itself has unlimited retries (when `self.retry=True`). The disk
operations it calls have their own independent 15-attempt retry. This creates nested retry
logic: an outer unlimited-retry loop calling inner 15-attempt-retry functions.

### 6.2 kive_watcher.py → update_qai.py (via queue)

```python
# kive_watcher enqueues work
self.qai_upload_queue.put((results_path, pipeline_group))

# update_qai dequeues and retries independently
upload_loop(qai_upload_queue, ...)
```

Each module handles its own retries independently. The queue decouples submission from
processing, so a failure in the update_qai processing does not directly retry the kive_watcher
submission.

### 6.3 update_qai.py → disk_operations.py (double-wrapped)

```python
def write_done_flag():
    disk_operations.write_text(done_qai_flag, status)
    # write_text is @disk_retry decorated (max 15)

retry_operation(write_done_flag, ...)  # outer: max 100
```

`write_text` is wrapped by both the `@disk_retry` decorator (15 attempts) and by
`retry_operation` (100 attempts). This double-wrapping is conservative: if 15 consecutive
disk-level attempts fail and raise, the outer `retry_operation` will catch that raised
exception and retry the entire `write_done_flag` call from scratch.

### 6.4 Queue-Based Time-Gated Retries

`find_samples()` feeds events into queues via `send_event()` / `send_qai_event()`. If the
queue is full, the event is not added and `is_complete = False` is set. The scan loop will
retry during the next scheduled scan cycle rather than immediately — a time-based retry
strategy rather than a tight retry loop.

### 6.5 Shared Configuration

```python
# kive_watcher.py (defines constants)
MINIMUM_RETRY_WAIT = timedelta(seconds=5)
MAXIMUM_RETRY_WAIT = timedelta(days=1)

# disk_operations.py (redefines with same values)
MINIMUM_RETRY_WAIT = timedelta(seconds=5)
MAXIMUM_RETRY_WAIT = timedelta(days=1)

# update_qai.py (imports from kive_watcher)
from micall.monitor.kive_watcher import wait_for_retry
```

---

## 7. Special Cases and Edge Scenarios

### 7.1 Sample Sheet Reading Retry

In `find_sample_groups()`, only `IOError` is retried:

```python
for attempt_count in count(1):
    try:
        file_names = list_fastq_file_names(...)
        sample_groups = list(find_groups(file_names, str(sample_sheet_path)))
        break  # Success
    except IOError:
        wait_for_retry(attempt_count, start_time)
        continue
    except Exception:
        # Log and write to errorprocessing, then break
        break
```

Non-IOError exceptions are logged to `run_path / "errorprocessing"` and break the retry
loop, treating them as unrecoverable input errors.

### 7.2 QAI Upload Coverage Flag Handling

In `upload_loop()`, the done flag is always written regardless of upload outcome:

```python
has_coverage_file = retry_operation(check_coverage_file, ...)

if has_coverage_file:
    try:
        process_folder(...)
        status = 'succeeded'
    except Exception:
        status = 'failed'
else:
    status = 'skipped'

# Always written with retry, preventing re-processing
retry_operation(write_done_flag, ...)
```

Even a failed or skipped upload advances the pipeline by writing the completion flag.

### 7.3 Session Refresh on KiveClientException

`kive_retry()` implements a single-retry-on-credential-failure pattern:

```python
def kive_retry(self, target: Callable[[], T]) -> T:
    try:
        return target()
    except KiveClientException:
        self.session.login(self.config.kive_user, self.config.kive_password)
        return target()
```

This assumes `KiveClientException` indicates stale credentials. It refreshes the session and
retries immediately with no backoff. If the second attempt also fails, the exception
propagates to the caller.

### 7.4 Selective Retry in upload_filter_quality()

```python
for attempt_count in count(1):
    try:
        reader = InterOpReader(folder_watcher.run_folder)
        write_phix_csv(quality_csv, records, read_lengths)
        break
    except OSError as ex:
        wait_for_retry(attempt_count, start_time)  # retry
    except Exception:
        logger.error("Finding error metrics in %s", ..., exc_info=True)
        disk_operations.write_text(folder_watcher.run_folder / "errorprocessing", ...)
        return  # no retry
```

Only `OSError` is retried; any other exception causes an immediate return with an error file
written.

---

## 8. Architectural Patterns and Design Principles

### 8.1 Defensive Design

- **Multiple Retry Layers**: Queue attempts + method retries + decorator retries
- **Exponential Backoff**: Prevents thundering herd and respects system saturation
- **Time-Aware Logging**: Distinguishes transient from persistent issues
- **Graceful Degradation**: Skipped uploads still mark the completion flag

### 8.2 Separation of Concerns

Each module handles retries only within its own domain:

- `kive_watcher.py`: Kive API interaction + retry coordination
- `disk_operations.py`: Disk I/O with a consistent retry policy
- `update_qai.py`: QAI database updates with a generic retry wrapper

### 8.3 Error Strategy

| Scenario | Strategy |
|----------|----------|
| Transient network failure | Retry with exponential backoff |
| Disk full / permission denied | Retry with exponential backoff |
| Malformed input data | Fail fast (no retry) |
| Stale session credentials | Refresh and retry once |
| QAI server 5xx error | Retry with exponential backoff |
| Missing file (graceful case) | Mark complete, no retry |

---

## 9. Monitoring and Debugging

### 9.1 Log Messages

**At INFO level** (within first hour of failures):

```
Waiting <time> before retrying.
<operation_name> failed (attempt X/MAX_ATTEMPTS)
Disk operation <op> failed, waiting <time> before retrying.
```

**At WARNING/ERROR level** (after 1 hour, or on final failure):

```
Waiting <time> before retrying.   [with exc_info=True]
<operation_name> failed after <MAX_ATTEMPTS> attempts.
Disk operation <op> failed after <max_attempts> attempts: <error>
```

### 9.2 Exception Chaining

In `update_qai.retry_operation()`, the original exception is preserved:

```python
raise RuntimeError(...) from last_exception
```

This allows debuggers and log handlers to inspect the root cause.

### 9.3 Error Files

On unrecoverable failures in `kive_watcher.py`, a file named `errorprocessing` is written to
the run folder so that downstream processes can detect the failure without parsing logs.

---

## 10. Summary Table

| Module | Function | Mechanism | Max Attempts | Backoff | Exceptions Retried |
|--------|----------|-----------|:------------:|---------|-------------------|
| kive_watcher | `find_samples` | Loop | ∞ | Exponential | All (parameter flag) |
| kive_watcher | `find_sample_groups` | `count(1)` | ∞ | Exponential | `IOError` |
| kive_watcher | `add_sample_group` | `count(1)` | ∞ | Exponential | All |
| kive_watcher | `poll_runs` | `count(1)` | ∞ | Exponential | All |
| kive_watcher | `kive_retry` | Try-catch | 2 | None | `KiveClientException` |
| kive_watcher | `upload_filter_quality` | `count(1)` | ∞ | Exponential | `OSError` |
| disk_operations | `@disk_retry` decorated | Decorator | 15 | Exponential | `OSError`, `IOError` |
| disk_operations | `disk_file_operation` | Context manager | 15 | Exponential | `OSError`, `IOError` |
| update_qai | `retry_operation` | Generic wrapper | 100 | Exponential | All exceptions |
| update_qai | `retry_qai_login` | `retry_operation` | 100 | Exponential | All |
| update_qai | `retry_find_run` | `retry_operation` | 100 | Exponential | All |

---

## Conclusion

The MiCall watcher retry system represents a mature, layered approach to handling failures in
a distributed environment. It combines:

- **Exponential backoff** to respect system resources and avoid thundering herd
- **Multiple independent retry layers** for defense in depth
- **Context-aware logging** to distinguish transient from persistent failures
- **Domain-specific retry policies** (15 for disk, 100 for QAI uploads, unlimited for scheduler loops)
- **Clear state tracking** with attempt counts and elapsed time
- **Graceful degradation** where partial success still advances progress

The system sustains robust operation despite transient network failures, disk access issues,
and API timeouts, while failing fast on unrecoverable input errors.
