"""Benchmarking utilities."""

import logging
from time import time

logging.basicConfig(filename='benchmark.log', level=logging.INFO)


def timeit(func):
    """Time function duration."""

    def function_timer(*args, **kwargs):
        t0 = time()
        value = func(*args, **kwargs)
        delta_t = time() - t0
        # logging.info(f"{func.__name__} -- {delta_t}s")
        return value

    return function_timer
