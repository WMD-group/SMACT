from __future__ import annotations

import os
import shutil
import sys
import tempfile
from importlib.util import find_spec
from typing import Any
from unittest.mock import patch

import pandas as pd
import pytest
import requests
from pymatgen.core import SETTINGS, Composition

from smact import Element
from smact.data_loader import (
    lookup_element_oxidation_states_custom,
    lookup_element_shannon_radius_data_extendedML,
    lookup_element_sse2015_data,
)
from smact.screening import SmactFilterOutputs, smact_filter
from smact.structure_prediction.database import _validate_table_name
from smact.utils.composition import comp_maker, composition_dict_maker, formula_maker, parse_formula
from smact.utils.crystal_space import download_compounds_with_mp_api, generate_composition_with_smact
from smact.utils.oxidation import ICSD24OxStatesFilter
from smact.utils.species import _parse_spec_old, parse_spec

MP_URL = "https://api.materialsproject.org"
MP_API_AVAILABLE = bool(find_spec("mp_api"))

# Lazily evaluated at first use to avoid network calls during test collection.
_skip_mprester_tests: bool | None = None


def _should_skip_mprester() -> bool:
    global _skip_mprester_tests  # noqa: PLW0603
    if _skip_mprester_tests is None:
        try:
            _skip_mprester_tests = requests.get(MP_URL, timeout=60).status_code != 200
        except (ModuleNotFoundError, ImportError, requests.exceptions.ConnectionError):
            _skip_mprester_tests = True
    return _skip_mprester_tests


files_dir = os.path.join(os.path.dirname(os.path.realpath(__file__)), "files")
TEST_ICSD_OX_STATES = os.path.join(files_dir, "oxidation_states_icsd24_consensus.txt")
TEST_ICSD_OX_STATES_W_ZERO = os.path.join(files_dir, "oxidation_states_icsd24_consensus_w_0.txt")


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture
def mock_filter_output():
    return [
        (("Fe", "O"), (2, -2), (1, 1)),
        (("Fe", "O"), (1, 1)),
        (("Fe", "Fe", "O"), (2, 3, -2), (1, 2, 4)),
    ]


@pytest.fixture
def smact_filter_output():
    return smact_filter(
        els=[Element("Li"), Element("Ge"), Element("P"), Element("S")],
        stoichs=[[10], [1], [2], [12]],
    )


@pytest.fixture
def ox_filter():
    return ICSD24OxStatesFilter()


@pytest.fixture
def test_ox_states_content():
    with open(TEST_ICSD_OX_STATES) as f:
        return f.read()


@pytest.fixture
def test_ox_states_w_zero_content():
    with open(TEST_ICSD_OX_STATES_W_ZERO) as f:
        return f.read()


# ---------------------------------------------------------------------------
# TestComposition
# ---------------------------------------------------------------------------


def test_parse_formula():
    formulas = ["Li10GeP2S12", "Mg0.5O0.5", "CaMg(CO3)2"]

    LGPS = parse_formula(formulas[0])
    assert isinstance(LGPS, dict)
    for el_sym, ammt in LGPS.items():
        assert isinstance(el_sym, str)
        assert isinstance(ammt, float)
    assert LGPS["Li"] == 10
    assert LGPS["Ge"] == 1
    assert LGPS["P"] == 2
    assert LGPS["S"] == 12

    MgO = parse_formula(formulas[1])
    assert isinstance(MgO, dict)
    assert MgO["Mg"] == 0.5
    assert MgO["O"] == 0.5

    dolomite = parse_formula(formulas[2])
    assert isinstance(dolomite, dict)
    assert dolomite["Ca"] == 1
    assert dolomite["Mg"] == 1
    assert dolomite["C"] == 2
    assert dolomite["O"] == 6


def test_comp_maker(mock_filter_output, smact_filter_output):
    comp1 = comp_maker(mock_filter_output[0])
    comp2 = comp_maker(mock_filter_output[1])
    comp3 = comp_maker(mock_filter_output[2])
    smact_output_item: Any = smact_filter_output[1]
    comp4 = comp_maker(smact_output_item)
    for comp in [comp1, comp2, comp3, comp4]:
        assert isinstance(comp, Composition)
    assert Composition("FeO") == comp2
    assert Composition({"Fe2+": 1, "O2-": 1}) == comp1
    assert Composition({"Fe2+": 1, "Fe3+": 2, "O2-": 4}) == comp3
    assert Composition({"Li+": 10, "Ge4+": 1, "P5+": 2, "S2-": 12}) == comp4


def test_formula_maker(mock_filter_output, smact_filter_output):
    form1 = formula_maker(mock_filter_output[0])
    form2 = formula_maker(mock_filter_output[1])
    form3 = formula_maker(mock_filter_output[2])
    smact_output_item: Any = smact_filter_output[1]
    form4 = formula_maker(smact_output_item)
    assert form1 == "FeO"
    assert form2 == "FeO"
    assert form1 == form2
    assert form3 == "Fe3O4"
    assert form4 == "Li10Ge(PS6)2"


def test_composition_dict_maker(mock_filter_output):
    # 3-element tuple (symbols, ox_states, stoichs) - returns species keys
    dict1 = composition_dict_maker(mock_filter_output[0])
    assert isinstance(dict1, dict)
    assert "Fe2+" in dict1
    assert "O2-" in dict1
    assert dict1["Fe2+"] == 1.0
    assert dict1["O2-"] == 1.0

    # 2-element tuple (symbols, stoichs) - returns element keys
    dict2 = composition_dict_maker(mock_filter_output[1])
    assert isinstance(dict2, dict)
    assert "Fe" in dict2
    assert "O" in dict2


def test_smact_filter_return_output_default():
    els = [Element("Li"), Element("O")]
    comps = smact_filter(els, threshold=2)
    assert isinstance(comps, list)
    assert len(comps) > 0
    # Default should return tuples
    assert isinstance(comps[0], tuple)


def test_smact_filter_return_output_formula():
    els = [Element("Li"), Element("O")]
    comps = smact_filter(els, threshold=2, return_output=SmactFilterOutputs.formula)
    assert isinstance(comps, list)
    assert len(comps) > 0
    for comp in comps:
        assert isinstance(comp, str)
    assert "Li2O" in comps


def test_smact_filter_return_output_dict():
    els = [Element("Li"), Element("O")]
    comps = smact_filter(els, threshold=2, return_output=SmactFilterOutputs.composition_dict)
    assert isinstance(comps, list)
    assert len(comps) > 0
    for comp in comps:
        assert isinstance(comp, dict)


def test_smact_filter_return_output_species_not_unique():
    els = [Element("Li"), Element("O")]
    # Formula output
    formulas = smact_filter(
        els,
        threshold=2,
        species_unique=False,
        return_output=SmactFilterOutputs.formula,
    )
    assert isinstance(formulas, list)
    assert len(formulas) > 0
    for f in formulas:
        assert isinstance(f, str)

    # Dict output
    dicts = smact_filter(
        els,
        threshold=2,
        species_unique=False,
        return_output=SmactFilterOutputs.composition_dict,
    )
    assert isinstance(dicts, list)
    assert len(dicts) > 0
    for d in dicts:
        assert isinstance(d, dict)


def test_smact_filter_invalid_return_output():
    els = [Element("Li"), Element("O")]
    bad_return_output: Any = "invalid"
    with pytest.raises(ValueError):
        smact_filter(els, threshold=2, return_output=bad_return_output)


# ---------------------------------------------------------------------------
# TestCrystalSpace
# ---------------------------------------------------------------------------


def test_convert_formula():
    combinations = [Element("Li"), Element("O")]
    expected_formulas = ["Li2O2", "LiO2", "Li2O", "Li2O2"]
    compounds = generate_composition_with_smact.convert_formula(combinations=combinations, num_elements=2, max_stoich=2)
    assert expected_formulas == compounds


@pytest.mark.parametrize("ox_states", ["smact14", "icsd24"])
def test_generate_composition_with_smact(ox_states, tmp_path):
    save_dir = str(tmp_path / "binary" / "df_binary_label.pkl")
    oxidation_states_sets_dict = {
        "smact14": {"smact_allowed": 388},
        "icsd24": {"smact_allowed": 342},
    }
    smact_df = generate_composition_with_smact.generate_composition_with_smact(
        num_elements=2,
        max_stoich=3,
        max_atomic_num=20,
        save_path=save_dir,
        oxidation_states_set=ox_states,
    )
    assert isinstance(smact_df, pd.DataFrame)
    assert len(smact_df) == 1330
    assert smact_df["smact_allowed"].sum() == oxidation_states_sets_dict[ox_states]["smact_allowed"]
    # Check if the data was saved to disk
    assert os.path.exists(save_dir)


@pytest.mark.parametrize("ox_states", ["smact14", "icsd24"])
def test_generate_composition_with_smact_custom(ox_states, tmp_path):
    save_dir = str(tmp_path / "binary" / "df_binary_label.pkl")
    oxidation_states_sets_dict = {
        "smact14": {"smact_allowed": 388},
        "icsd24": {"smact_allowed": 342},
    }
    smact_df = generate_composition_with_smact.generate_composition_with_smact_custom(
        num_elements=2,
        max_stoich=3,
        max_atomic_num=20,
        save_path=save_dir,
        oxidation_states_set=ox_states,
    )
    assert isinstance(smact_df, pd.DataFrame)
    assert len(smact_df) == 1330
    assert smact_df["smact_allowed"].sum() == oxidation_states_sets_dict[ox_states]["smact_allowed"]
    # Check if the data was saved to disk
    assert os.path.exists(save_dir)


@pytest.mark.integration
@pytest.mark.skipif(
    (
        sys.platform == "win32"
        or not (os.environ.get("MP_API_KEY") or SETTINGS.get("PMG_MAPI_KEY"))
        or not MP_API_AVAILABLE
    ),
    reason="Test requires MP_API_KEY and fails on Windows due to filepath issues.",
)
def test_download_compounds_with_mp_api():
    if _should_skip_mprester():
        pytest.skip("Materials Project API endpoint not reachable.")
    save_mp_dir = "data/binary/mp_data"
    if MP_API_AVAILABLE:
        download_compounds_with_mp_api.download_mp_data(
            mp_api_key=os.environ.get("MP_API_KEY"),
            num_elements=2,
            max_stoich=1,
            save_dir=save_mp_dir,
        )

    # Check if the data was downloaded
    assert os.path.exists(save_mp_dir)
    assert len(os.listdir(save_mp_dir)) > 0

    # Clean up
    shutil.rmtree(save_mp_dir)


# ---------------------------------------------------------------------------
# OxidationStatesTest
# ---------------------------------------------------------------------------


def test_oxidation_states_filter(ox_filter):
    assert isinstance(ox_filter.ox_states_df, pd.DataFrame)
    threshold = 10
    filtered_df = ox_filter.filter(consensus=threshold)

    assert isinstance(filtered_df, pd.DataFrame)
    assert filtered_df.columns.tolist() == ["element", "oxidation_state"]


def test_oxidation_states_write(ox_filter, test_ox_states_content, test_ox_states_w_zero_content):
    comment = "Testing writing of ICSD 24 oxidation states list."

    with tempfile.TemporaryDirectory() as tmpdir:
        filename = os.path.join(tmpdir, "test_ox_states")
        ox_filter.write(
            filename,
            comment=comment,
            consensus=3,
            include_zero=False,
            commonality="low",
        )

        assert os.path.exists(f"{filename}.txt")
        with open(f"{filename}.txt") as f:
            content = f.read()
        assert comment in content
        assert content == test_ox_states_content

        filename_w_zero = os.path.join(tmpdir, "test_ox_states_w_zero")
        ox_filter.write(
            filename_w_zero,
            comment=comment,
            consensus=3,
            include_zero=True,
            commonality="low",
        )

        assert os.path.exists(f"{filename_w_zero}.txt")
        with open(f"{filename_w_zero}.txt") as f:
            content = f.read()
        assert comment in content
        assert content == test_ox_states_w_zero_content


def test_get_species_list(ox_filter):
    # Test with default parameters
    species_list = ox_filter.get_species_list()
    assert isinstance(species_list, list)
    assert len(species_list) > 0  # Ensure the list is not empty

    # Test with include_zero=True
    species_list_with_zero = ox_filter.get_species_list(include_zero=True)
    assert isinstance(species_list_with_zero, list)
    assert len(species_list_with_zero) > 0

    # Test with include_one_oxidation_state=True
    species_list_with_one = ox_filter.get_species_list(include_one_oxidation_state=True)
    assert isinstance(species_list_with_one, list)
    assert len(species_list_with_one) > 0

    # Test with different commonality levels
    species_list_low = ox_filter.get_species_list(commonality="low")
    assert isinstance(species_list_low, list)
    assert len(species_list_low) > 0

    species_list_medium = ox_filter.get_species_list(commonality="medium")
    assert isinstance(species_list_medium, list)
    assert len(species_list_medium) > 0

    species_list_high = ox_filter.get_species_list(commonality="high")
    assert isinstance(species_list_high, list)
    assert len(species_list_high) > 0

    # Test with a specific consensus threshold
    species_list_threshold = ox_filter.get_species_list(consensus=5)
    assert isinstance(species_list_threshold, list)
    assert len(species_list_threshold) > 0

    # Test commonality="main" (should return species with max proportion for each element)
    species_list_main = ox_filter.get_species_list(commonality="main")
    assert isinstance(species_list_main, list)
    assert len(species_list_main) > 0

    # Get the original dataframe for comparison
    df_main = ox_filter.get_species_occurrences_df()
    assert isinstance(df_main, pd.DataFrame)
    grouped = df_main.groupby("element")["species_proportion (%)"]
    max_proportions = grouped.max()

    # Verify that only species with maximum proportion for each element are included
    for species in species_list_main:
        element = species.split("+")[0].split("-")[0].rstrip("0123456789")  # Extract element from species string
        filtered = df_main[df_main["species"] == species]
        assert isinstance(filtered, pd.DataFrame)
        species_proportion = filtered["species_proportion (%)"].iloc[0]
        assert species_proportion == max_proportions[element]

    # Test specific commonality threshold
    threshold = 50.0  # Testing with 50% threshold
    species_list_threshold = ox_filter.get_species_list(commonality=threshold)
    assert isinstance(species_list_threshold, list)

    # Verify that all species meet the threshold requirement
    df_threshold = ox_filter.get_species_occurrences_df()
    assert isinstance(df_threshold, pd.DataFrame)
    for species in species_list_threshold:
        filtered = df_threshold[df_threshold["species"] == species]
        assert isinstance(filtered, pd.DataFrame)
        proportion = filtered["species_proportion (%)"].iloc[0]
        assert proportion >= threshold

    # Test that "main" returns different results than threshold-based filtering
    species_list_main = ox_filter.get_species_list(commonality="main")
    species_list_high_threshold = ox_filter.get_species_list(commonality=90.0)
    assert set(species_list_main) != set(species_list_high_threshold)


def test_oxidation_states_filter_species_occurrences(ox_filter):
    species_occurrences_df = ox_filter.get_species_occurrences_df(consensus=1)
    assert isinstance(species_occurrences_df, pd.DataFrame)
    assert species_occurrences_df.columns.tolist() == [
        "element",
        "species",
        "results_count",
        "species_proportion (%)",
    ]
    assert species_occurrences_df.shape == (490, 4)
    first_row = species_occurrences_df.iloc[0]
    assert isinstance(first_row, pd.Series)
    assert first_row["species"] == "O2-"
    assert first_row["results_count"] == 116910


# ---------------------------------------------------------------------------
# TestCompositionEdgeCases
# ---------------------------------------------------------------------------


def test_parse_formula_invalid_chars():
    with pytest.raises(ValueError, match="invalid formula"):
        parse_formula("123invalid")


def test_parse_formula_empty_string():
    result = parse_formula("")
    assert result == {}


# ---------------------------------------------------------------------------
# TestOxidationEdgeCases
# ---------------------------------------------------------------------------


def test_commonality_type_error(ox_filter):
    bad_commonality: Any = [1, 2, 3]
    with pytest.raises(TypeError, match="commonality must be"):
        ox_filter.filter(commonality=bad_commonality)


def test_get_species_occurrences_unsorted(ox_filter):
    occurrences = ox_filter.get_species_occurrences_df(sort_by_occurrences=False)
    assert isinstance(occurrences, pd.DataFrame)
    assert len(occurrences) > 0


def test_filter_numeric_commonality(ox_filter):
    filtered = ox_filter.filter(commonality=25.0)
    assert isinstance(filtered, pd.DataFrame)
    assert len(filtered) > 0


def test_get_species_occurrences_include_zero(ox_filter):
    occurrences = ox_filter.get_species_occurrences_df(include_zero=True)
    assert isinstance(occurrences, pd.DataFrame)
    assert len(occurrences) > 0


def test_get_species_list_value_error_path(ox_filter):
    # Inject a corrupted row where oxidation_state is a non-integer string.
    # The code splits on spaces, so use a single token that fails int().
    corrupted_df = pd.DataFrame(
        {
            "element": ["Na", "Na"],
            "oxidation_state": ["bad", "1"],
        }
    )
    with patch.object(ox_filter, "filter", return_value=corrupted_df):
        species_list = ox_filter.get_species_list()

    # int("bad") -> ValueError -> skipped via continue
    # Only "Na1+" from the second row should appear
    assert len(species_list) == 1
    assert "Na+" in species_list


# ---------------------------------------------------------------------------
# TestDatabaseEdgeCases
# ---------------------------------------------------------------------------


def test_invalid_table_name():
    with pytest.raises(ValueError, match="Invalid table name"):
        _validate_table_name("DROP TABLE; --")


# ---------------------------------------------------------------------------
# TestSpeciesParsing
# ---------------------------------------------------------------------------


def test_parse_spec_old_no_element_symbol_raises():
    with pytest.raises(ValueError, match="Invalid species string"):
        _parse_spec_old("123")


def test_parse_spec_old_negative_ox_state():
    # "Fe(-2)": main regex fails (parens), _parse_spec_old finds digit '2' and '-'
    ele, charge = _parse_spec_old("Fe(-2)")
    assert ele == "Fe"
    assert charge == -2


def test_parse_spec_old_zero_with_digit_zero():
    # "Fe0": no sign -> main regex fails -> _parse_spec_old
    ele, charge = parse_spec("Fe0")
    assert ele == "Fe"
    assert charge == 0


def test_parse_spec_old_bare_plus():
    # "Fe(+)": parens block main regex match; no digit; '+' present
    ele, charge = _parse_spec_old("Fe(+)")
    assert ele == "Fe"
    assert charge == 1


def test_parse_spec_old_bare_minus():
    # "Fe(-)": parens block main regex; no digit; '-' present
    ele, charge = _parse_spec_old("Fe(-)")
    assert ele == "Fe"
    assert charge == -1


# ---------------------------------------------------------------------------
# TestDataLoaderWarnings
# ---------------------------------------------------------------------------


def test_warn_on_missing_logs_debug(test_ox_states_path, caplog):
    import logging

    with caplog.at_level(logging.DEBUG, logger="smact.data_loader"):
        result = lookup_element_oxidation_states_custom("Xx", test_ox_states_path)
    assert result is None
    assert any("not found" in msg for msg in caplog.messages)


def test_shannon_radii_extendedML_found_and_missing():
    # Valid element: exercises _load_shannon_radii_extendedML (lines 381-398)
    data = lookup_element_shannon_radius_data_extendedML("Fe")
    assert data is not None
    assert isinstance(data, list)

    # Unknown element: exercises missing-symbol warning path (lines 447-451)
    result = lookup_element_shannon_radius_data_extendedML("Xx")
    assert result is None


def test_sse2015_missing_symbol():
    result = lookup_element_sse2015_data("Xx")
    assert result is None
