"""Benchmarking utilities."""

import functools
import logging
from statistics import mean
from time import time

DELIM_LENGTH = 15

logging.basicConfig(filename="benchmark.log", level=logging.INFO)


def timeit(_func=None, *, fname=None, n=1, delim=False):
    """Time function duration."""

    def decorator_timeit(func):
        @functools.wraps(func)
        def wrapper_timeit(*args, **kwargs):
            if delim:
                logging.info("-" * DELIM_LENGTH)

            times = []
            for _ in range(n):
                t0 = time()
                value = func(*args, **kwargs)
                times.append(time() - t0)

            logging.info(
                f"{func.__name__} -- Average over {n} repeats = {mean(times)}s"
            )

            if delim:
                logging.info("-" * DELIM_LENGTH)
            return value

        return wrapper_timeit

    fname = "benchmark.log" if not fname else fname
    logging.basicConfig(filename=fname)

    if _func is None:
        return decorator_timeit
    else:
        return decorator_timeit(_func)
