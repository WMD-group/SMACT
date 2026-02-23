"""Compatibility utilities for supporting multiple Python versions."""

from __future__ import annotations

from enum import Enum

try:
    from enum import StrEnum
except ImportError:

    class StrEnum(str, Enum):
        """Backport of Python 3.11's StrEnum for Python 3.10."""

        def __str__(self):
            return str(self.value)
