"""Utility functions for SMACT generation of compositions."""

from __future__ import annotations

import itertools
import multiprocessing
import warnings
from functools import partial
from pathlib import Path

import pandas as pd
from pymatgen.core import Composition
from tqdm import tqdm

from smact import Element, ordered_elements
from smact.screening import smact_filter

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
        formula_dict = dict(zip(symbols, counts, strict=False))
        formula = Composition(formula_dict).reduced_formula
        local_compounds.append(formula)
    return local_compounds


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
        oxidation_states_set (str): the oxidation states set to use. Options are "smact14", "icsd16", "icsd24", "pymatgen_sp" or a filepath to a custom oxidation states list. For reproducing the Faraday Discussions results, use "smact14".

    Returns:
        df (pd.DataFrame): A DataFrame of SMACT-generated compositions with boolean smact_allowed column.

    """
    # 1. generate all possible combinations of elements
    print("#1. Generating all possible combinations of elements...")

    elements = [Element(element) for element in ordered_elements(1, max_atomic_num)]
    combinations = list(itertools.combinations(elements, num_elements))
    print(f"Number of generated combinations: {len(list(combinations))}")

    # 2. generate all possible stoichiometric combinations
    print("#2. Generating all possible stoichiometric combinations...")

    pool = multiprocessing.Pool(processes=multiprocessing.cpu_count() if num_processes is None else num_processes)
    compounds = list(
        tqdm(
            pool.imap_unordered(
                partial(
                    convert_formula,
                    num_elements=num_elements,
                    max_stoich=max_stoich,
                ),
                combinations,
            ),
            total=len(combinations),
        )
    )

    pool.close()
    pool.join()
    # Flatten the list of lists into a single list
    compounds = [item for sublist in compounds for item in sublist]

    print(f"Number of generated compounds: {len(compounds)}")
    compounds = list(set(compounds))
    print(f"Number of generated compounds (unique): {len(compounds)}")

    # 3. filter compounds with smact
    print("#3. Filtering compounds with SMACT...")
    elements_pauling = [
        Element(element) for element in ordered_elements(1, max_atomic_num) if Element(element).pauling_eneg is not None
    ]  # omit elements without Pauling electronegativity (e.g., He, Ne, Ar, ...)
    compounds_pauling = list(itertools.combinations(elements_pauling, num_elements))

    pool = multiprocessing.Pool(processes=multiprocessing.cpu_count() if num_processes is None else num_processes)
    results = list(
        tqdm(
            pool.imap_unordered(
                partial(smact_filter, threshold=max_stoich, oxidation_states_set=oxidation_states_set),
                compounds_pauling,
            ),
            total=len(compounds_pauling),
        )
    )
    pool.close()
    pool.join()

    # 4. make data frame of results
    print("#4. Making data frame of results...")
    # make dataframework with index is compound and columns are boolean smact results
    smact_allowed = []

    for result in results:
        for res in result:
            symbols_stoich = zip(res[0], res[2], strict=False)
            composition_dict = dict(symbols_stoich)
            smact_allowed.append(Composition(composition_dict).reduced_formula)
    smact_allowed = list(set(smact_allowed))
    print(f"Number of compounds allowed by SMACT: {len(smact_allowed)}")

    df = pd.DataFrame(data=False, index=compounds, columns=["smact_allowed"])
    df.loc[smact_allowed, "smact_allowed"] = True

    if save_path is not None:
        Path(save_path).parent.mkdir(parents=True, exist_ok=True)
        df.to_pickle(save_path)
        print(f"Saved to {save_path}")

    return df
