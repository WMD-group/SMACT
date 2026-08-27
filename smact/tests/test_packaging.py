"""Checks that the data files SMACT needs at runtime are present in the installed package.

SMACT reads almost all of its reference data from files under ``smact/data``, which reach
users only if they are declared in ``[tool.setuptools.package-data]``. The rest of the test
suite cannot catch a missing declaration: CI installs the project with ``uv sync``, which is
an editable install, so ``smact.data_directory`` resolves to the checkout, where every
committed file is present whether or not it would ship. A dropped glob would therefore leave
the suite green while ``pip install smact`` produced a package that raises FileNotFoundError
on first use.

The ``packaging`` workflow runs this module against an installed wheel, in a directory that
does not contain the source tree, so that the checks below see the artefact users get. They
are written to pass in either context, so they also run as part of the normal suite.

See GH issue #643, where the default lambda table was read out of pymatgen's install
directory and vanished when pymatgen reorganised its own package data.
"""

from __future__ import annotations

import tomllib
from pathlib import Path

import pytest

import smact

# smact.data_directory is built from Path(__file__).resolve(), so resolve here too rather
# than comparing it against an unresolved smact.__file__: on macOS /tmp is a symlink to
# /private/tmp and the two spellings of the same directory would not compare equal.
PACKAGE_DIR = Path(smact.__file__).resolve().parent
DATA_DIR = Path(smact.data_directory)

# Mirrors the patterns in [tool.setuptools.package-data] in pyproject.toml. Each one must
# still match something once installed: a category losing all of its files is the failure
# mode a deleted glob produces. test_declared_globs_match_pyproject keeps the two in step.
DECLARED_DATA_GLOBS = (
    "data/*.txt",
    "data/*.csv",
    "data/*.json",
    "data/*.xlsx",
    "data/species_rep/*.json",
)

# Declared in pyproject.toml but matching nothing today, so it is exempt from the check that
# every pattern is populated.
UNUSED_DATA_GLOBS = ("data/*.data",)

# Data files that smact/ reads at runtime, so shipping without them breaks a documented
# code path rather than merely losing a convenience. Add to this when new data is read
# from disk; it is not meant to enumerate every file in smact/data.
REQUIRED_DATA_FILES = (
    "element_data.txt",
    "element_valence_modified.csv",
    "hhi.txt",
    "magpie.csv",
    "ordered_periodic.txt",
    "oxidation_state_probability_table.json",
    "oxidation_states.txt",
    "oxidation_states_SP.txt",
    "oxidation_states_icsd.txt",
    "oxidation_states_icsd24_counts.json",
    "oxidation_states_icsd24_filtered.txt",
    "oxidation_states_wiki.txt",
    "shannon_radii.csv",
    "shannon_radii_ML_extended.csv",
    "SSE.csv",
    "SSE_2015.csv",
    "SSE_Pauling.csv",
    "species_rep/ion_embedding_M3GNet-MP-2023.11.1-oxi-Eform_cosine_similarity.json",
    "species_rep/ion_embedding_M3GNet-MP-2023.11.1-oxi-band_gap_cosine_similarity.json",
    "species_rep/skipspecies_20221028_319ion_dim200_cosine_similarity.json",
)


def _repo_pyproject() -> Path | None:
    """Return the repo's pyproject.toml, or None when running against an installed package."""
    candidate = PACKAGE_DIR.parent / "pyproject.toml"
    return candidate if candidate.is_file() else None


def test_data_directory_is_inside_the_package():
    # smact.data_directory is derived from smact.__file__, so if it points anywhere other
    # than the installed package the checks below would be testing the wrong tree.
    assert DATA_DIR.is_dir(), f"{DATA_DIR} is missing from the installed package"
    assert DATA_DIR == PACKAGE_DIR / "data"


@pytest.mark.parametrize("pattern", DECLARED_DATA_GLOBS)
def test_declared_data_globs_are_populated(pattern):
    matches = sorted(PACKAGE_DIR.glob(pattern))
    assert matches, f"no installed file matches the declared package-data pattern {pattern!r}"


@pytest.mark.parametrize("relative_path", REQUIRED_DATA_FILES)
def test_required_data_file_is_installed(relative_path):
    path = DATA_DIR / relative_path
    assert path.is_file(), f"{relative_path} is missing from the installed package"
    # A truncated or placeholder file would still satisfy is_file().
    assert path.stat().st_size > 0, f"{relative_path} is installed but empty"


def test_py_typed_marker_is_installed():
    # Without the marker, type checkers ignore SMACT's annotations in downstream projects.
    assert (PACKAGE_DIR / "py.typed").is_file()


def test_declared_globs_match_pyproject():
    """Keep DECLARED_DATA_GLOBS in step with pyproject.toml.

    Skipped when running against an installed package, where pyproject.toml is absent.
    """
    pyproject = _repo_pyproject()
    if pyproject is None:
        pytest.skip("running against an installed package, so pyproject.toml is unavailable")

    with pyproject.open("rb") as f:
        config = tomllib.load(f)

    declared = config["tool"]["setuptools"]["package-data"]["smact"]
    data_globs = [pattern for pattern in declared if pattern.startswith("data/") and pattern not in UNUSED_DATA_GLOBS]

    assert sorted(data_globs) == sorted(DECLARED_DATA_GLOBS), (
        "pyproject.toml and DECLARED_DATA_GLOBS disagree about which data files ship. "
        "Update this test alongside [tool.setuptools.package-data]."
    )


def test_element_data_loads_from_the_installed_package():
    # Exercises element_data.txt, the oxidation state tables and shannon_radii.csv through
    # the public API, rather than only checking that the files are on disk.
    iron = smact.Element("Fe")
    assert iron.number == 26
    assert iron.oxidation_states

    fe2 = smact.Species("Fe", 2)
    assert fe2.oxidation == 2
    # None here would mean shannon_radii.csv shipped but held no rows for this species.
    assert fe2.average_ionic_radius is not None
    assert fe2.average_ionic_radius > 0


def test_species_embedding_tables_load():
    # The paths are module-level constants in doper, so this checks the files the shipped
    # code actually points at.
    from smact.dopant_prediction import doper

    for path in (
        doper.SKIPSPECIES_COSINE_SIM_PATH,
        doper.SPECIES_M3GNET_MP2023_EFORM_COSINE_PATH,
        doper.SPECIES_M3GNET_MP2023_GAP_COSINE_PATH,
    ):
        assert Path(path).is_file(), f"{path} is missing from the installed package"
