"""
Shared download machinery: retrying what is worth retrying, and writing through a `.part` file.

Two entry points, because the package fetches two different kinds of thing. `download_file` is for
bulk files: it is unpaced and holds no state, so many can run at once. `download_api` is for REST
endpoints: it paces its requests and remembers, per host, how long the pacing has had to grow.

Every response lands on `<out_fpath>.part` and is moved onto its final name only once it is
complete, so an interrupted download cannot leave a truncated file under a name a cache would
trust.

`urlcleanup()` is deliberately absent. Its only effects are unlinking the temp files `urlretrieve`
creates when no filename is given -- none, here, since every call passes one -- and discarding the
cached global opener, which costs a rebuilt handler chain per request and races between threads.
"""

import logging
import os
import threading
from http.client import HTTPException
from time import sleep
from urllib.error import HTTPError
from urllib.error import URLError
from urllib.parse import urlparse
from urllib.request import urlretrieve

# Timing statuses worth another attempt. Every other 4xx describes the request itself and will
# fail identically however many times it is repeated.
TRANSIENT_HTTP_STATUSES = frozenset({408, 425, 429})

# Seconds to wait before each request to a given host, keyed by host name. Raised by
# `download_api` whenever a request to that host has to be retried, and never lowered, so an API
# that has asked to be slowed down stays slowed for the life of the process.
_HOST_DELAYS = {}
_HOST_DELAYS_LOCK = threading.Lock()


def discard_partial(fpath):
    """
    Delete a partially written download, ignoring the case where it never got created.

    Args:
        fpath (str): Path to remove.

    Returns:
        None
    """
    if os.path.exists(fpath):
        os.remove(fpath)


def is_transient_error(error):
    """
    Decide whether a failed request is worth another attempt.

    An HTTP status counts as transient only when it describes the server or the moment: any 5xx,
    plus the 408, 425 and 429 timing statuses. Every other 4xx describes the request itself, so a
    missing entry costs one attempt rather than the whole retry budget. Failures below HTTP --
    name resolution, refused or reset connections, timeouts and short reads -- are all transient.
    Anything else, including a local filesystem failure while writing the file, is permanent.

    Args:
        error (BaseException): The exception raised by the failed attempt.

    Returns:
        bool: True if repeating the request could plausibly succeed.
    """
    # HTTPError first: it subclasses URLError, so the isinstance below would otherwise claim it.
    if isinstance(error, HTTPError):
        return error.code >= 500 or error.code in TRANSIENT_HTTP_STATUSES
    return isinstance(error, (URLError, HTTPException, TimeoutError, ConnectionError))


def reset_host_delays():
    """
    Forget the recorded per-host delays, returning every host to its caller's base delay.

    Returns:
        None
    """
    with _HOST_DELAYS_LOCK:
        _HOST_DELAYS.clear()


def get_host_delay(url, base_delay):
    """
    Report the delay currently being applied to a host, seeding it on first use.

    Args:
        url (str): Any URL on the host.
        base_delay (float): Delay to record for a host not seen before.

    Returns:
        float: Seconds to wait before the next request to that host.
    """
    host = urlparse(url).netloc
    with _HOST_DELAYS_LOCK:
        return _HOST_DELAYS.setdefault(host, base_delay)


def escalate_host_delay(url, max_delay):
    """
    Double the delay applied to a host, up to a ceiling.

    Args:
        url (str): Any URL on the host.
        max_delay (float): Ceiling on the returned delay.

    Returns:
        float: The host's new delay, in seconds.
    """
    host = urlparse(url).netloc
    with _HOST_DELAYS_LOCK:
        delay = min(_HOST_DELAYS.get(host, 0.0) * 2, max_delay)
        _HOST_DELAYS[host] = delay
        return delay


def write_through_part(url, out_fpath, transform, log_extra):
    """
    Fetch `url` and place it at `out_fpath`, writing through a `.part` file.

    The response is written beside the destination and moved onto it with `os.replace` only once
    it is complete, so a failed fetch leaves nothing under the final name. With a `transform`, the
    response goes to a second scratch file and only the transform's output is placed. Every
    scratch file is removed before returning, whether or not the fetch succeeded.

    Args:
        url (str): Address to fetch.
        out_fpath (str): Final path for the file.
        transform (callable): Called as `transform(src_fpath, dst_fpath)` to convert the response
            before it is placed, or None to place it unchanged.
        log_extra (dict): `extra` mapping for the log records.

    Returns:
        None

    Raises:
        Exception: Whatever the fetch or the transform raised.
    """
    part_fpath = f"{out_fpath}.part"
    raw_fpath = f"{out_fpath}.raw.part" if transform else part_fpath
    try:
        logging.debug(f"Fetching {url} -> {raw_fpath}", extra=log_extra)
        urlretrieve(url, raw_fpath)
        if transform:
            transform(raw_fpath, part_fpath)
        os.replace(part_fpath, out_fpath)
    except Exception:
        discard_partial(part_fpath)
        raise
    finally:
        if transform:
            discard_partial(raw_fpath)


def download_file(url, out_fpath, transform=None, max_retries=5, base_delay=0.25, max_delay=30.0, log_extra=None):
    """
    Download `url` to `out_fpath`, retrying transient failures with exponential backoff.

    Requests are not paced: a delay is waited only after a failed attempt, starting at `base_delay`
    and doubling up to `max_delay`. A permanent failure such as a 404 gives up after the first
    attempt. Nothing is shared between calls, so any number may run concurrently.

    Args:
        url (str): Address to fetch.
        out_fpath (str): Final path for the file.
        transform (callable, optional): Called as `transform(src_fpath, dst_fpath)` to convert the
            response before it is placed. Defaults to None, which places it unchanged.
        max_retries (int): Attempts allowed before giving up. Defaults to 5.
        base_delay (float): Seconds to wait after the first failed attempt. Defaults to 0.25.
        max_delay (float): Ceiling on the doubling delay. Defaults to 30.0.
        log_extra (dict, optional): `extra` mapping for the log records, which must carry a
            "stage" key. Defaults to None, which logs under a generic download stage.

    Returns:
        bool: True once the file is in place, False if every attempt failed.
    """
    extra = log_extra or {"stage": "Downloading"}
    delay = base_delay
    for attempt in range(1, max_retries + 1):
        try:
            write_through_part(url, out_fpath, transform, extra)
            return True
        except Exception as e:
            if attempt == max_retries or not is_transient_error(e):
                logging.warning(f"Giving up on {url} after {attempt} attempt(s): {e}", extra=extra)
                return False
            logging.debug(
                f"Attempt {attempt}/{max_retries} failed for {url} ({e}); retrying in {delay:.2f}s",
                extra=extra,
            )
            sleep(delay)
            delay = min(delay * 2, max_delay)
    return False


def download_api(url, out_fpath, max_retries=5, base_delay=0.25, max_delay=30.0, log_extra=None):
    """
    Download a REST response to `out_fpath`, pacing requests to the host it came from.

    The host's current delay is waited before every request, starting at `base_delay`. A transient
    failure doubles that delay, capped at `max_delay`, and it is never lowered again, so the
    backoff one request needed becomes the pace of every later request to the same host. The delay
    is shared by every caller in the process and outlives any one of them; `reset_host_delays`
    clears it.

    Args:
        url (str): Address to fetch.
        out_fpath (str): Final path for the response.
        max_retries (int): Attempts allowed before giving up. Defaults to 5.
        base_delay (float): Delay recorded for a host not seen before. Defaults to 0.25.
        max_delay (float): Ceiling on the host's delay. Defaults to 30.0.
        log_extra (dict, optional): `extra` mapping for the log records, which must carry a
            "stage" key. Defaults to None, which logs under a generic download stage.

    Returns:
        bool: True once the response is in place, False if every attempt failed.
    """
    extra = log_extra or {"stage": "Downloading"}
    for attempt in range(1, max_retries + 1):
        sleep(get_host_delay(url, base_delay))
        try:
            write_through_part(url, out_fpath, None, extra)
            return True
        except Exception as e:
            if attempt == max_retries or not is_transient_error(e):
                logging.warning(f"Giving up on {url} after {attempt} attempt(s): {e}", extra=extra)
                return False
            delay = escalate_host_delay(url, max_delay)
            logging.debug(
                f"Attempt {attempt}/{max_retries} failed for {url} ({e}); "
                f"pacing {urlparse(url).netloc} at {delay:.2f}s",
                extra=extra,
            )
    return False
