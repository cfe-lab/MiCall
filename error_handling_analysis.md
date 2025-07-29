# MiCall Watcher Error Handling Analysis

## Executive Summary

The MiCall Watcher module (`micall_watcher.py`) is a bioinformatics pipeline monitoring system that watches for new sequencing data, processes it through various analysis pipelines on a Kive server, and uploads results to a QAI database. This analysis examines the module's architecture, network error handling mechanisms, and provides recommendations for improvement.

## System Architecture

### Core Components

1. **`micall_watcher.py`** - Main entry point and orchestration
2. **`kive_watcher.py`** - Kive API interaction and pipeline management
3. **`sample_watcher.py`** - Sample-level processing logic
4. **`update_qai.py`** - QAI database upload functionality

### Threading Model

The system uses a multi-threaded architecture:
- **Main thread**: Orchestrates the overall process
- **Sample finder thread**: Scans filesystem for new samples (`find_samples`)
- **QAI upload thread**: Handles database uploads (`upload_loop`)

### Data Flow

```
Raw Data Scanning → Sample Processing → Kive Pipeline Execution → Result Collection → QAI Upload
```

## Current Error Handling Mechanisms

### 1. Retry Logic with Exponential Backoff

The system implements sophisticated retry logic in `kive_watcher.py`:

```python
def wait_for_retry(attempt_count, is_logged=True):
    delay = calculate_retry_wait(MINIMUM_RETRY_WAIT, MAXIMUM_RETRY_WAIT, attempt_count)
    if is_logged:
        logger.error('Waiting %s before retrying.', delay, exc_info=True)
    sleep(delay.total_seconds())

def calculate_retry_wait(min_wait, max_wait, attempt_count):
    min_seconds = int(min_wait.total_seconds())
    seconds = min_seconds * (2 ** (attempt_count - 1))  # Exponential backoff
    seconds = min(seconds, max_wait.total_seconds())
    return timedelta(seconds=seconds)
```

**Configuration:**
- Minimum retry wait: 5 seconds
- Maximum retry wait: 1 day
- Uses exponential backoff (2^attempt_count)

### 2. Network-Specific Error Handling

#### Kive API Calls
- **`kive_retry()`**: Single retry with session refresh for Kive API failures
- **`check_session()`**: Automatic session management and re-authentication
- **HTTP Adapter**: Uses `HTTPAdapter(max_retries=20)` for low-level HTTP retries

#### Filesystem Operations
Special handling for network drive issues:
```python
# In find_samples()
is_logged = not isinstance(ex, BlockingIOError) or attempt_count > 1
wait_for_retry(attempt_count, is_logged)
```

### 3. Error State Management

#### Run-level Failures
- Failed runs are tracked and can trigger re-runs
- Cancellation handling (`state == 'X'`) with automatic re-launch
- Failure states (`state == 'F'`) raise `KiveRunFailedException`

#### Folder-level Failures  
- Error files created: `errorprocessing` flag
- Failed folders are skipped in subsequent scans
- Quality control failures prevent further processing

#### Sample-level Failures
- Individual samples can fail without stopping other samples
- Failed samples are tracked with `is_failed` flag
- Completed samples are tracked to avoid reprocessing

### 4. Exception Categories and Handling

#### Network Errors
1. **Kive API failures** → `kive_retry()` with session refresh
2. **HTTP timeouts** → HTTP adapter retries
3. **Filesystem I/O** → Exponential backoff with special `BlockingIOError` handling

#### Data Processing Errors
1. **Sample sheet parsing** → Creates error file, skips folder
2. **Quality metric parsing** → Creates error file, aborts processing
3. **File format issues** → Logged and skipped

#### System Errors
1. **Authentication failures** → Logged at startup and during retries
2. **Configuration errors** → Early validation with clear error messages
3. **Resource constraints** → Queue management and active job limiting

## Strengths of Current Implementation

### 1. Comprehensive Retry Strategy
- Multi-level retry mechanisms (HTTP, API, application)
- Intelligent backoff prevents overwhelming servers
- Configurable retry limits prevent infinite loops

### 2. Graceful Degradation
- Failed samples don't stop processing of other samples
- Failed folders are marked to prevent repeated attempts
- System can recover from transient network issues

### 3. State Persistence
- Error states are persisted to filesystem
- Processing progress is tracked
- System can resume after restarts

### 4. Logging and Observability
- Structured logging with context
- Exception details captured with `exc_info=True`
- Clear error messages for different failure modes

## Areas for Improvement

### 1. Network Error Classification

**Current Issues:**
- Generic exception handling doesn't distinguish between error types
- All network errors get the same retry treatment
- No circuit breaker pattern for persistent failures

**Recommendations:**
```python
from requests.exceptions import (
    ConnectionError, Timeout, HTTPError, 
    RequestException, ReadTimeout, ConnectTimeout
)

class NetworkErrorHandler:
    def __init__(self):
        self.circuit_breaker = CircuitBreaker()
        
    def handle_network_error(self, error, operation_context):
        if isinstance(error, (ConnectTimeout, ConnectionError)):
            # Infrastructure issues - longer backoff
            return self.handle_infrastructure_error(error, operation_context)
        elif isinstance(error, ReadTimeout):
            # Server overload - shorter backoff
            return self.handle_server_overload(error, operation_context)
        elif isinstance(error, HTTPError):
            # Application errors - different strategy
            return self.handle_application_error(error, operation_context)
        else:
            # Unknown errors - conservative approach
            return self.handle_unknown_error(error, operation_context)
```

### 2. Circuit Breaker Pattern

**Current Issues:**
- System will retry indefinitely for persistent failures
- No mechanism to temporarily disable failing services
- Resource waste on hopeless retry attempts

**Recommendations:**
```python
class CircuitBreaker:
    def __init__(self, failure_threshold=5, timeout=300):
        self.failure_count = 0
        self.failure_threshold = failure_threshold
        self.timeout = timeout
        self.last_failure_time = None
        self.state = 'CLOSED'  # CLOSED, OPEN, HALF_OPEN
    
    def call(self, func, *args, **kwargs):
        if self.state == 'OPEN':
            if time.time() - self.last_failure_time > self.timeout:
                self.state = 'HALF_OPEN'
            else:
                raise CircuitBreakerOpenError()
        
        try:
            result = func(*args, **kwargs)
            self.on_success()
            return result
        except Exception as e:
            self.on_failure()
            raise
```

### 3. Enhanced Monitoring and Alerting

**Current Issues:**
- Limited metrics on error rates and patterns
- No proactive alerting for persistent issues
- Difficult to distinguish between temporary and serious problems

**Recommendations:**
```python
class ErrorMetrics:
    def __init__(self):
        self.error_counts = defaultdict(int)
        self.error_rates = defaultdict(list)
        
    def record_error(self, error_type, operation, timestamp=None):
        if timestamp is None:
            timestamp = time.time()
        
        self.error_counts[error_type] += 1
        self.error_rates[error_type].append(timestamp)
        
        # Alert if error rate exceeds threshold
        if self.get_error_rate(error_type, window=300) > 0.1:  # 10% error rate
            self.send_alert(error_type, operation)
```

### 4. Improved Configuration Management

**Current Issues:**
- Hardcoded timeout and retry values
- No runtime configuration updates
- Limited environment-specific tuning

**Recommendations:**
```python
@dataclass
class NetworkConfig:
    min_retry_wait: int = 5
    max_retry_wait: int = 86400  # 1 day
    max_retries: int = 20
    connection_timeout: int = 30
    read_timeout: int = 300
    circuit_breaker_threshold: int = 5
    circuit_breaker_timeout: int = 300
    
    @classmethod
    def from_environment(cls):
        return cls(
            min_retry_wait=int(os.getenv('MICALL_MIN_RETRY_WAIT', 5)),
            max_retry_wait=int(os.getenv('MICALL_MAX_RETRY_WAIT', 86400)),
            # ... other configurations
        )
```

### 5. Better Error Context and Recovery

**Current Issues:**
- Limited context in error messages
- No suggested recovery actions
- Difficult to correlate related failures

**Recommendations:**
```python
class OperationContext:
    def __init__(self, operation_type, sample_id=None, run_id=None):
        self.operation_type = operation_type
        self.sample_id = sample_id
        self.run_id = run_id
        self.start_time = time.time()
        self.attempt_count = 0
        
    def log_error(self, error, suggested_action=None):
        context = {
            'operation': self.operation_type,
            'sample_id': self.sample_id,
            'run_id': self.run_id,
            'duration': time.time() - self.start_time,
            'attempt': self.attempt_count,
            'error_type': type(error).__name__,
            'error_message': str(error),
            'suggested_action': suggested_action
        }
        logger.error("Operation failed", extra=context, exc_info=True)
```

### 6. Graceful Shutdown and Resource Cleanup

**Current Issues:**
- Abrupt shutdown may leave resources in inconsistent state
- No mechanism to finish in-progress operations
- Limited cleanup of temporary files and connections

**Recommendations:**
```python
class GracefulShutdownHandler:
    def __init__(self):
        self.shutdown_requested = False
        signal.signal(signal.SIGTERM, self.handle_shutdown)
        signal.signal(signal.SIGINT, self.handle_shutdown)
    
    def handle_shutdown(self, signum, frame):
        logger.info("Shutdown requested, finishing current operations...")
        self.shutdown_requested = True
        
    def should_continue(self):
        return not self.shutdown_requested
```

## Implementation Priority

### High Priority (Immediate)
1. **Network error classification** - Distinguish between different network failure types
2. **Circuit breaker pattern** - Prevent resource waste on persistent failures
3. **Enhanced logging context** - Better error correlation and debugging

### Medium Priority (Next Quarter)
1. **Configuration management** - Runtime-configurable timeouts and retry policies
2. **Monitoring and metrics** - Error rate tracking and alerting
3. **Graceful shutdown** - Clean resource cleanup

### Low Priority (Future)
1. **Predictive failure detection** - ML-based anomaly detection
2. **Auto-scaling retry policies** - Dynamic adjustment based on system load
3. **Distributed tracing** - End-to-end request tracking

## Testing Recommendations

### Unit Tests
- Mock network failures to test retry logic
- Verify circuit breaker state transitions
- Test configuration edge cases

### Integration Tests
- Simulate various network conditions (slow, intermittent, failed)
- Test end-to-end error propagation
- Verify data consistency after failures

### Chaos Engineering
- Randomly inject network failures during normal operations
- Test system behavior under sustained outages
- Verify recovery after infrastructure repairs

## Conclusion

The MiCall Watcher system demonstrates a solid foundation for error handling with comprehensive retry mechanisms and state management. However, there are significant opportunities for improvement in network error classification, circuit breaking, and observability. The recommendations provided would enhance system reliability, reduce resource waste, and improve operational visibility.

The proposed improvements follow industry best practices for distributed systems and would make the system more resilient to the inherent unreliability of network operations while providing better insights into system health and performance.
