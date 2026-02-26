from __future__ import annotations

import os

import pytest

import smact


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
