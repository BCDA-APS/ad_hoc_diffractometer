# Diagnosing a Threading Object That Doesn't Terminate in IPython

## 1. List all live threads

The first step — see what threads are running and whether any are stuck:

```python
import threading
for t in threading.enumerate():
    print(t.name, t.daemon, t.is_alive())
```

A non-daemon thread that stays alive after its work is done will block interpreter
shutdown and can cause unstable IPython behavior (kernel hangs, cell re-execution
anomalies, etc.).

---

## 2. Get a stack trace of every live thread

```python
import sys
import threading
import traceback

for thread_id, frame in sys._current_frames().items():
    name = {t.ident: t.name for t in threading.enumerate()}.get(thread_id, "unknown")
    print(f"\n--- Thread: {name} (id={thread_id}) ---")
    traceback.print_stack(frame)
```

This shows exactly where each thread is blocked — a thread stuck in `time.sleep()`,
`queue.get()`, `socket.recv()`, or a `Lock.acquire()` is a strong signal.

---

## 3. Check for daemon vs non-daemon threads

Non-daemon threads prevent interpreter shutdown.  If you see an unexpected
non-daemon thread:

```python
for t in threading.enumerate():
    if not t.daemon and t is not threading.main_thread():
        print(f"Non-daemon thread still alive: {t.name}")
```

Setting `t.daemon = True` **before** `t.start()` is the fix at source.

---

## 4. Use `faulthandler` to catch a hang

If the kernel actually hangs (rather than just misbehaving), `faulthandler` dumps
all thread stacks on a signal:

```python
import faulthandler
import signal
faulthandler.register(signal.SIGUSR1)  # then from a shell: kill -USR1 <pid>
```

Or enable it unconditionally (dumps on SIGSEGV/crash too):

```python
faulthandler.enable()
```

---

## 5. Identify the source — what spawns the thread?

Common culprits in IPython sessions:

| Source | Symptom |
|---|---|
| `concurrent.futures.ThreadPoolExecutor` not shut down | Workers stay alive waiting for tasks |
| `threading.Timer` not cancelled | Timer thread runs after it should have stopped |
| Background `queue.Queue` consumer | Thread blocks on `queue.get()` indefinitely |
| `asyncio` event loop in a thread | Loop runs in a daemon thread but holds resources |
| Matplotlib interactive backend | GUI event loop thread |
| Third-party library (e.g. `watchdog`, `kafka`) | Background monitor threads |

---

## 6. Reproduce and isolate

In IPython, run progressively simpler code to find the minimum reproducer:

```python
# After the suspected operation:
import threading
print([t.name for t in threading.enumerate()])
```

Run this **before** and **after** the suspected call.  Any new thread that appears
and does not disappear is the candidate.

---

## 7. `%reset` and kernel restart as controls

- `%reset` clears the namespace but does **not** kill threads.  If instability
  persists after `%reset`, a thread is the likely cause.
- A full kernel restart (`Kernel > Restart`) kills all threads.  If instability
  clears after a restart but returns after a specific operation, that operation
  is spawning the thread.

---

## Summary of diagnostic order

1. `threading.enumerate()` — see what is alive
2. `sys._current_frames()` — see where each thread is blocked
3. `faulthandler` — catch hard hangs
4. Before/after comparison — identify which call spawns the thread
5. Check for non-daemon status and missing `shutdown()` / `cancel()` / `join()`
   calls at source

---

## Bluesky RunEngine — specific diagnostics

The RunEngine (RE) runs an `asyncio` event loop in a dedicated thread.  The
plan completing successfully but RE not returning control is a known symptom
with several distinct causes.

### Check RE state

```python
print(RE.state)
```

Expected after a completed plan: `'idle'`.  If it shows `'running'` or
`'paused'`, the RE has not finished its internal cleanup.

| State | Meaning |
|---|---|
| `'idle'` | Ready — no issue here |
| `'running'` | Still executing — a device callback or a `wait()` is blocking |
| `'paused'` | Waiting for `RE.resume()` or `RE.abort()` — most common cause of "hangs" |
| `'panicked'` | An unhandled exception in the event loop thread |

### Force exit from a paused state

```python
RE.abort()    # abandon the current plan cleanly
# or
RE.stop()     # request a stop (softer than abort)
# or
RE.halt()     # immediate halt, no cleanup callbacks
```

### Inspect the event loop thread

The RE's internal thread is named `'RunEngine'` or similar:

```python
import threading
re_threads = [t for t in threading.enumerate() if 'RunEngine' in t.name or 'run_engine' in t.name.lower()]
for t in re_threads:
    print(t.name, t.daemon, t.is_alive())
```

Then get its stack trace (see section 2 above) to see exactly where it is
blocked.

### Common bluesky-specific causes

| Cause | Diagnostic | Fix |
|---|---|---|
| Device `set()` or `trigger()` never completes its `Status` object | Stack shows RE waiting in `asyncio` event loop on a `Future` | Ensure the device marks its `Status` as done: `status.set_finished()` or `status.set_exception(...)` |
| `bps.wait()` waiting on a group that never completes | Stack shows `wait()` blocking | Check that all `bps.trigger()` / `bps.set()` calls in the group have matching completions |
| Callback (subscription) raises an exception silently | RE state goes to `'paused'` or hangs | Wrap callbacks in try/except; check `RE.exceptions` if available |
| `RE.subscribe()` callback blocks (e.g. slow file I/O) | Stack shows callback thread blocked | Make callbacks non-blocking; use a queue to offload work |
| `asyncio` coroutine in a device never yields | Event loop starved | Ensure device coroutines use `await asyncio.sleep(0)` to yield |
| Multiple RE instances sharing the same event loop | Unpredictable blocking | Use one RE per session; check `RE._loop` is not shared |

### Inspect pending asyncio tasks

```python
import asyncio

loop = RE._loop  # the RE's internal event loop
if loop.is_running():
    # From inside an async context or using run_coroutine_threadsafe:
    fut = asyncio.run_coroutine_threadsafe(
        _list_tasks(loop), loop
    )

async def _list_tasks(loop):
    for task in asyncio.all_tasks(loop):
        print(task.get_name(), task.done(), task.cancelled())
        task.print_stack()
```

Or more simply from IPython (which itself runs an event loop):

```python
for task in asyncio.all_tasks(RE._loop):
    print(task)
```

### Check device Status objects

A `Status` that never completes is the single most common cause.  Instrument
a device temporarily:

```python
from ophyd import Device
original_set = MyDevice.set

def patched_set(self, *args, **kwargs):
    status = original_set(self, *args, **kwargs)
    status.add_callback(lambda s: print(f"Status done: {s.success}"))
    return status

MyDevice.set = patched_set
```

### Recommended RE construction for IPython sessions

Use `call_returns_result=True` and a `context_managers` entry to ensure the
RE always cleans up:

```python
from bluesky import RunEngine
RE = RunEngine({})
RE.subscribe(best_effort_callback)
```

Always pair with:

```python
from bluesky.utils import install_kicker
install_kicker()   # prevents event loop starvation in IPython/Jupyter
```

`install_kicker()` (or its Jupyter equivalent `install_nb_kicker()`) runs
periodic event loop ticks so that the RE's `asyncio` loop does not stall
waiting for the IPython event loop to yield.

---

## Subscribers

A subscriber (added via `RE.subscribe()`) is called synchronously inside the
RE's event loop thread for every document the RE emits (`start`, `descriptor`,
`event`, `stop`).  If a subscriber blocks, raises silently, or accumulates
state incorrectly, the RE can appear to hang or produce wrong results.

### List all active subscribers

```python
# Each token maps to a (name, callable) pair
for token, (name, func) in RE._subscribes.items():
    print(f"token={token}  name={name!r}  func={func}")
```

The `token` is the integer returned by `RE.subscribe()` and can be passed to
`RE.unsubscribe(token)`.

### Check whether a subscriber is blocking

Add a timing wrapper around a suspected subscriber:

```python
import time
import functools

def timed(func):
    @functools.wraps(func)
    def wrapper(name, doc):
        t0 = time.perf_counter()
        result = func(name, doc)
        elapsed = time.perf_counter() - t0
        if elapsed > 0.1:   # flag anything slower than 100 ms
            print(f"[SLOW SUBSCRIBER] {func.__name__} took {elapsed:.3f}s on {name!r}")
        return result
    return wrapper

# Wrap and re-subscribe:
RE.unsubscribe(old_token)
new_token = RE.subscribe(timed(my_callback))
```

### Check whether a subscriber raises silently

By default the RE catches and swallows exceptions in subscribers.  Instrument
to surface them:

```python
def catching(func):
    @functools.wraps(func)
    def wrapper(name, doc):
        try:
            return func(name, doc)
        except Exception as exc:
            print(f"[SUBSCRIBER EXCEPTION] {func.__name__} on {name!r}: {exc!r}")
            raise
    return wrapper
```

### Isolate by unsubscribing one at a time

```python
# Save all tokens at a known-good point:
tokens = list(RE._subscribes.keys())

# Then unsubscribe each in turn and re-run the plan:
for token in tokens:
    RE.unsubscribe(token)
    RE(my_plan())   # does it return now?
    # re-subscribe before trying the next one
```

Bisecting this way identifies which subscriber is the culprit.

### Common subscriber problems

| Problem | Symptom | Fix |
|---|---|---|
| Slow file I/O (HDF5, TIFF write) | RE returns but very late | Move I/O to a background thread or `queue.Queue` |
| Unhandled exception swallowed by RE | Plan finishes but data is wrong or incomplete | Wrap in try/except and log; check databroker insert errors |
| Subscriber holds a lock also held by the main thread | Deadlock — RE never returns | Avoid shared locks between subscriber and main thread |
| Accumulating large in-memory structures | Memory growth; eventually OOM | Clear or flush subscriber state between runs |
| `stop` document never received by subscriber | Subscriber state never resets | Ensure plan always emits `stop`; check for plan that raises before stop |

---

## Suspenders

A suspender (added via `RE.install_suspender()`) monitors a signal and pauses
the RE when a condition is met, resuming when it clears.  The RE will print a
message when it suspends, so a suspender is usually visible — but it can be
subtle if the signal flickers or if the print is missed in a busy session.

### List all active suspenders

```python
for suspender in RE.suspenders:
    print(suspender)
    print(f"  tripped={suspender.tripped}  should_resume={suspender.should_resume}")
```

### Check whether any suspender is currently tripped

```python
tripped = [s for s in RE.suspenders if s.tripped]
if tripped:
    print("Tripped suspenders:", tripped)
else:
    print("No suspenders are currently tripped")
```

### Force-clear a tripped suspender (for diagnosis only)

```python
# Remove all suspenders temporarily:
RE.clear_suspenders()

# Re-run the plan — if it completes, a suspender was the cause.
RE(my_plan())
```

### Common suspender problems

| Problem | Symptom | Fix |
|---|---|---|
| Signal starts at a bad value before devices are ready | RE suspends immediately on first run | Add a pre-run settle time or check signal initialisation |
| `resume_cond` never becomes true | RE suspends and never resumes | Check that the signal that clears the suspender is actually changing |
| Suspender installed on a PV that is disconnected | `tripped` is True because the value is unknown | Check PV connection: `signal.connected` |
| Multiple suspenders with conflicting conditions | RE oscillates between suspended and running | Log each suspender's state separately |

---

## Other causes to check

If threading, subscribers, and suspenders are all ruled out, consider the
following.

### Plan generator never exhausts

A plan is a Python generator.  If it contains an infinite loop or a
`yield from` that never returns, the RE runs forever:

```python
# Instrument the plan to count yields:
def counting(plan):
    n = 0
    for msg in plan:
        n += 1
        if n % 100 == 0:
            print(f"Message {n}: {msg.command}")
        yield msg

RE(counting(my_plan()))
```

### A `Msg` with `command='wait'` targets a non-existent group

If a plan issues `bps.wait(group='mygroup')` but nothing was ever started
under that group name, the RE waits indefinitely:

```python
# Print every message the plan emits before running it:
for msg in my_plan():
    print(msg)
```

Look for `wait` messages and verify the `group` kwarg matches a prior
`trigger`, `set`, or `kickoff` message with the same group name.

### Callback chain depth (recursion in the event loop)

Deep recursion inside the `asyncio` event loop raises `RecursionError` which
the RE may catch internally.  Check:

```python
import sys
print(sys.getrecursionlimit())   # default 1000
# Temporarily raise if needed for diagnosis:
sys.setrecursionlimit(2000)
```

### Ophyd device not fully connected

A device whose PV is not connected will cause `set()` or `trigger()` to
return a `Status` that never completes.  Check before running:

```python
from ophyd.utils import make_all_devices
for name, dev in RE.devices.items() if hasattr(RE, 'devices') else []:
    if not dev.connected:
        print(f"NOT CONNECTED: {name}")

# Or check a specific device directly:
print(my_device.connected)
print(my_device.read_attrs)
my_device.wait_for_connection(timeout=5)
```

### Databroker / handler errors at `stop` document

If the databroker insert raises on the `stop` document, the RE may appear to
hang waiting for confirmation:

```python
# Subscribe a minimal stop-document logger before running:
def log_stop(name, doc):
    if name == 'stop':
        print(f"Stop document received: exit_status={doc.get('exit_status')}")

token = RE.subscribe(log_stop)
RE(my_plan())
RE.unsubscribe(token)
```

If "Stop document received" never prints, the plan is not completing at the
RE level — the issue is earlier in the execution, not in databroker.

### Summary checklist

Work through this list in order:

1. `RE.state` — is the RE actually stuck or just slow?
2. Unsubscribe all subscribers one at a time — does the plan return?
3. `RE.suspenders` — is any suspender tripped?
4. `RE.clear_suspenders()` then re-run — does it return?
5. Print all `Msg` objects from the plan — any `wait` with an unmatched group?
6. Check all device `.connected` — any PV disconnected?
7. Subscribe a `stop`-document logger — does stop ever arrive?
8. Thread stack traces — where is the RE event loop thread actually blocked?
