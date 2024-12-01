"""Utility functions for downloading compounds from the Materials Project."""

from __future__ import annotations

import itertools
import json
import string
import time
from collections import defaultdict
from pathlib import Path

from mp_api.client import MPRester
from pymatgen.core import Composition
from tqdm import tqdm


def download_mp_data(
    mp_api_key: str | None = None,
    num_elements: int = 2,
    max_stoich: int = 8,
    save_dir: str = "data/binary/mp_api",
    request_interval: float = 0.1,
):
    """
    Downloads Materials Project data all possible combinations of `num_elements` elements
    with atomic numbers.
    When chemical formula is same, the one with lowest energy above hull is saved.
    The data is saved to a specified directory.

    Args:
    ----
        mp_api_key (str, optional): the API key for Materials Project.
        num_elements (int, optional): the number of elements in each compound to consider.
            Defaults to 2.
        max_stoich (int, optional): the maximum integer of stoichiometric coefficient
            in chemical formula. Defaults to 8.
        save_dir (str, optional): the directory to save the downloaded data to.
            Defaults to "data/mp_api".
        request_interval (float, optional): the time interval between API requests, in seconds.
            Defaults to 1.

    Returns:
    -------
        None

    """
    # check if MP_API_KEY is set
    if mp_api_key is None:
        raise ValueError("Please set your MP_API_KEY in the environment variable.")
    # set save directory
    save_dir = Path(save_dir)
    save_dir.mkdir(parents=True, exist_ok=True)

    # make a list for all possible combinartions of formula anonymous
    symbols = string.ascii_uppercase
    formula_anonymous_list = []
    for stoichs in itertools.combinations_with_replacement(range(1, max_stoich + 1), num_elements):
        formula_dict = {symbols[i]: stoich for i, stoich in enumerate(stoichs)}
        formula_anonymous_list.append(Composition(formula_dict).reduced_formula)
    formula_anonymous_list = sorted(set(formula_anonymous_list))

    e_hull_dict = defaultdict(lambda: float("inf"))

    for formula_anonymous in tqdm(formula_anonymous_list):
        print(f"Downloading data for {formula_anonymous}...")
        # download data from MP
        with MPRester(mp_api_key) as mpr:
            docs = mpr.materials.summary.search(
                formula=formula_anonymous,
                fields=[
                    "formula_pretty",
                    "material_id",
                    "formula_anonymous",
                    "volume",
                    "density",
                    "density_atomic",
                    "energy_per_atom",
                    "formation_energy_per_atom",
                    "energy_above_hull",
                    "is_stable",
                    "band_gap",
                    "efermi",
                    "total_magnetization",
                    "structure",
                ],
            )
        # save data with lowest energy above hull
        for doc in docs:
            formula_pretty = doc.formula_pretty
            energy_above_hull = doc.energy_above_hull

            if (energy_above_hull) < e_hull_dict[formula_pretty]:
                e_hull_dict[formula_pretty] = energy_above_hull
                with open(save_dir / f"{formula_pretty}.json", "w") as f:
                    json.dump(doc.dict(), f)
        time.sleep(request_interval)
