"""Small environment guards for third-party network clients."""

from contextlib import contextmanager
import os


@contextmanager
def without_invalid_generic_proxy():
    """Temporarily hide the nonstandard ``socks://`` generic proxy.

    HTTPX accepts ``socks5://`` when its optional SOCKS dependency is
    installed, but some HPC login environments publish ``socks://`` through
    ``ALL_PROXY`` alongside working HTTP and HTTPS proxies.  Removing only the
    invalid generic setting lets each client use those protocol-specific
    proxies.  The caller's environment is restored after the request.
    """
    removed = {}
    for name in ('ALL_PROXY', 'all_proxy'):
        value = os.environ.get(name, '')
        if value.strip().lower().startswith('socks://'):
            removed[name] = os.environ.pop(name)
    try:
        yield
    finally:
        os.environ.update(removed)
