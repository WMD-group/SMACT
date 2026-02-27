from __future__ import annotations

import os
import sys

import pytest

import smact

# openTSNE C extension crashes with SIGILL on Python 3.13 (ubuntu),
# which kills the entire pytest process during collection.
# Skip the elementembeddings test file to prevent this.
collect_ignore_glob: list[str] = []
if sys.version_info >= (3, 13):
    collect_ignore_glob.append("*test_elementembeddings*")


@pytest.fixture
def element_na():
    return smact.Element("Na")


@pytest.fixture
def element_cl():
    return smact.Element("Cl")


@pytest.fixture
def element_fe():
    return smact.Element("Fe")


@pytest.fixture
def element_o():
    return smact.Element("O")


@pytest.fixture
def files_dir():
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), "files")
