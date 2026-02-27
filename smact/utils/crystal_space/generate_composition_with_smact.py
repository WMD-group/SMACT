"""Utility functions for SMACT generation of compositions."""

from __future__ import annotations

import itertools
import logging
import multiprocessing
import os
import warnings
from functools import partial
from pathlib import Path

import pandas as pd
from pymatgen.core import Composition
from tqdm import tqdm

from smact import Element, data_directory, ordered_elements
from smact.data_loader import lookup_element_oxidation_states_custom
from smact.screening import smact_filter

logger = logging.getLogger(__name__)

# Map short names used in smact_filter to their underlying data files so that
# generate_composition_with_smact_custom can accept the same named sets.
_NAMED_OX_SETS: dict[str, str] = {
    "smact14": os.path.join(data_directory, "oxidation_states.txt"),
    "icsd16": os.path.join(data_directory, "oxidation_states_icsd.txt"),
    "icsd24": os.path.join(data_directory, "oxidation_states_icsd24_filtered.txt"),
    "pymatgen_sp": os.path.join(data_directory, "oxidation_states_SP.txt"),
}

warnings.simplefilter(action="ignore", category=UserWarning)


def convert_formula(combinations: list, num_elements: int, max_stoich: int) -> list:
    """Convert combinations into chemical formula.

    Args:
        combinations (list): list of lists of smact.Element objects.
        num_elements (int): the number of elements in a compound.
        max_stoich (int): the maximum stoichiometric coefficient.

    Returns:
        local_compounds (list): A list of chemical formula.
    """
    symbols = [element.symbol for element in combinations]
    local_compounds = []
    for counts in itertools.product(range(1, max_stoich + 1), repeat=num_elements):
        formula_dict = dict(zip(symbols, counts, strict=True))
        formula = Composition(formula_dict).reduced_formula
        local_compounds.append(formula)
    return local_compounds


def _generate_unique_compounds(
    num_elements: int,
    max_stoich: int,
    max_atomic_num: int,
    num_processes: int | None,
) -> list[str]:
    """Steps 1 & 2: generate all unique reduced formulas.

    Args:
        num_elements: number of elements per compound.
        max_stoich: maximum stoichiometric coefficient.
        max_atomic_num: maximum atomic number to include.
        num_processes: number of worker processes (None = cpu_count).

    Returns:
        Deduplicated list of reduced chemical formulas.
    """
    logger.info("#1. Generating all possible combinations of elements...")
    elements = [Element(element) for element in ordered_elements(1, max_atomic_num)]
    combinations = list(itertools.combinations(elements, num_elements))
    logger.info("Number of generated combinations: %d", len(combinations))

    logger.info("#2. Generating all possible stoichiometric combinations...")
    pool = multiprocessing.Pool(processes=(multiprocessing.cpu_count() if num_processes is None else num_processes))
    compounds = list(
        tqdm(
            pool.imap_unordered(  # type: ignore[misc]
                partial(convert_formula, num_elements=num_elements, max_stoich=max_stoich),
                combinations,  # type: ignore[arg-type]
            ),
            total=len(combinations),
        )
    )
    pool.close()
    pool.join()

    compounds = [item for sublist in compounds for item in sublist]
    logger.info("Number of generated compounds: %d", len(compounds))
    compounds = sorted(set(compounds))
    logger.info("Number of generated compounds (unique): %d", len(compounds))
    return compounds


def _build_results_df(
    compounds: list[str],
    smact_results: list,
    save_path: str | None,
) -> pd.DataFrame:
    """Step 4: build and optionally persist the results DataFrame.

    Args:
        compounds: full list of candidate reduced formulas.
        smact_results: raw output from pool.imap_unordered over smact_filter.
        save_path: optional file path to pickle the DataFrame.

    Returns:
        DataFrame indexed by formula with a boolean ``smact_allowed`` column.
    """
    logger.info("#4. Making data frame of results...")
    smact_allowed = []
    for result in smact_results:
        for res in result:
            symbols_stoich = zip(res[0], res[2], strict=True)
            composition_dict = dict(symbols_stoich)
            smact_allowed.append(Composition(composition_dict).reduced_formula)
    smact_allowed = list(set(smact_allowed))
    logger.info("Number of compounds allowed by SMACT: %d", len(smact_allowed))

    results_df = pd.DataFrame(data=False, index=compounds, columns=["smact_allowed"])  # type: ignore[call-overload]
    results_df.loc[smact_allowed, "smact_allowed"] = True

    if save_path is not None:
        Path(save_path).parent.mkdir(parents=True, exist_ok=True)
        results_df.to_pickle(save_path)
        logger.info("Saved to %s", save_path)

    return results_df


def generate_composition_with_smact(
    num_elements: int = 2,
    max_stoich: int = 8,
    max_atomic_num: int = 103,
    num_processes: int | None = None,
    save_path: str | None = None,
    oxidation_states_set: str = "icsd24",
) -> pd.DataFrame:
    """
    Generate all possible compositions of a given number of elements and
    filter them with SMACT.

    Args:
        num_elements (int): the number of elements in a compound. Defaults to 2.
        max_stoich (int): the maximum stoichiometric coefficient. Defaults to 8.
        max_atomic_num (int): the maximum atomic number. Defaults to 103.
        num_processes (int): the number of processes to use. Defaults to None.
        save_path (str): the path to save the results. Defaults to None.
        oxidation_states_set (str): the oxidation states set to use. Options are "smact14", "icsd16", "icsd24", "pymatgen_sp". For reproducing the Faraday Discussions results, use "smact14". For custom oxidation states lists check generate_composition_with_smact_custom below.

    Returns:
        df (pd.DataFrame): A DataFrame of SMACT-generated compositions with boolean smact_allowed column.

    """
    compounds = _generate_unique_compounds(num_elements, max_stoich, max_atomic_num, num_processes)

    # 3. filter compounds with smact
    logger.info("#3. Filtering compounds with SMACT...")
    elements_pauling = [
        Element(element) for element in ordered_elements(1, max_atomic_num) if Element(element).pauling_eneg is not None
    ]  # omit elements without Pauling electronegativity (e.g., He, Ne, Ar, ...)
    compounds_pauling = list(itertools.combinations(elements_pauling, num_elements))

    pool = multiprocessing.Pool(processes=(multiprocessing.cpu_count() if num_processes is None else num_processes))
    results = list(
        tqdm(
            pool.imap_unordered(  # type: ignore[misc]
                partial(smact_filter, threshold=max_stoich, oxidation_states_set=oxidation_states_set),
                compounds_pauling,  # type: ignore[arg-type]
            ),
            total=len(compounds_pauling),
        )
    )
    pool.close()
    pool.join()

    return _build_results_df(compounds, results, save_path)


def generate_composition_with_smact_custom(
    num_elements: int = 2,
    max_stoich: int = 8,
    max_atomic_num: int = 103,
    num_processes: int | None = None,
    save_path: str | None = None,
    oxidation_states_set: str = "icsd24",
) -> pd.DataFrame:
    """
    Generate all possible compositions of a given number of elements and
    filter them with SMACT.

    Args:
        num_elements (int): the number of elements in a compound. Defaults to 2.
        max_stoich (int): the maximum stoichiometric coefficient. Defaults to 8.
        max_atomic_num (int): the maximum atomic number. Defaults to 103.
        num_processes (int): the number of processes to use. Defaults to None.
        save_path (str): the path to save the results. Defaults to None.
        oxidation_states_set (str): the path to the oxidation states file. Defaults to "icsd24".

    Returns:
        df (pd.DataFrame): A DataFrame of SMACT-generated compositions with boolean smact_allowed column.

    """
    compounds = _generate_unique_compounds(num_elements, max_stoich, max_atomic_num, num_processes)

    # 3. filter compounds with smact
    logger.info("#3. Filtering compounds with SMACT...")

    ox_filepath = _NAMED_OX_SETS.get(oxidation_states_set, oxidation_states_set)
    ox_states_custom = lookup_element_oxidation_states_custom("all", ox_filepath, copy=False)
    fr_eneg = Element("Fr").pauling_eneg
    elements_pauling = [
        Element(element)
        for element in ordered_elements(1, max_atomic_num)
        if ox_states_custom is not None
        and element in ox_states_custom  # type: ignore[operator]
        and Element(element).pauling_eneg is not None
        and (fr_eneg is None or Element(element).pauling_eneg >= fr_eneg)  # type: ignore[operator]
    ]
    compounds_pauling = list(itertools.combinations(elements_pauling, num_elements))

    pool = multiprocessing.Pool(processes=(multiprocessing.cpu_count() if num_processes is None else num_processes))
    results = list(
        tqdm(
            pool.imap_unordered(  # type: ignore[misc]
                partial(smact_filter, threshold=max_stoich, oxidation_states_set=ox_filepath),
                compounds_pauling,  # type: ignore[arg-type]
            ),
            total=len(compounds_pauling),
        )
    )
    pool.close()
    pool.join()

    return _build_results_df(compounds, results, save_path)
