"""Benchmarking utilities."""

from __future__ import annotations

import functools
import logging
from statistics import mean
from time import time

DELIM_LENGTH = 15

logger = logging.getLogger(__name__)

logging.basicConfig(filename="benchmark.log", level=logging.INFO)


def timeit(_func=None, *, fname=None, n=1, delim=False):
    """Time function duration."""

    def decorator_timeit(func):
        @functools.wraps(func)
        def wrapper_timeit(*args, **kwargs):
            if delim:
                logger.info("-" * DELIM_LENGTH)

            times = []
            value = None
            for _ in range(n):
                t0 = time()
                value = func(*args, **kwargs)
                times.append(time() - t0)

            logger.info("%s -- Average over %s repeats = %ss", func.__name__, n, mean(times))

            if delim:
                logger.info("-" * DELIM_LENGTH)
            return value

        return wrapper_timeit

    fname = fname if fname else "benchmark.log"
    logging.basicConfig(filename=fname)

    if _func is None:
        return decorator_timeit
    return decorator_timeit(_func)
