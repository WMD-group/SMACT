"""
Structure prediction implementation.

Todo:
----
    * Test with a fully populated database.
    * Implement n-ary substitution probabilities;
      at the moment, only zero- and single-species
      substitutions are considered.

"""

from __future__ import annotations

import itertools
from typing import TYPE_CHECKING

import numpy as np

from .utilities import parse_spec, unparse_spec

if TYPE_CHECKING:
    from collections.abc import Generator

    import pandas as pd

    from .database import StructureDB
    from .mutation import CationMutator
    from .structure import SmactStructure


class StructurePredictor:
    """
    Provides structure prediction functionality.

    Implements a statistically-based model for determining
    likely structures of a given composition, based on a
    database of known compositions and a lambda table of
    weights.

    Based on the algorithm presented in:
        Hautier, G., Fischer, C., Ehrlacher, V., Jain, A., and Ceder, G. (2011)
        Data Mined Ionic Substitutions for the Discovery of New Compounds.
        Inorganic Chemistry, 50(2), 656-663.
        `doi:10.1021/ic102031h <https://pubs.acs.org/doi/10.1021/ic102031h>`_

    """

    def __init__(self, mutator: CationMutator, struct_db: StructureDB, table: str) -> None:
        """
        Initialize class.

        Args:
        ----
            mutator: A :class:`CationMutator` for probability calculations.
            struct_db: A :class:`StructureDB` from which to read structures
                to attempt to mutate.
            table: The table to reference within the database

        """
        self.cm = mutator
        self.db = struct_db
        self.table = table

    def _try_substitution(
        self,
        parent: SmactStructure,
        diff_spec: tuple[str, int],
        diff_spec_str: str,
        diff_sub_probs: pd.Series,
        species: list[tuple[str, int]],
        thresh: float | None,
    ) -> tuple[SmactStructure, float, SmactStructure] | None:
        """Attempt a single unary substitution on a parent structure.

        Returns a (mutated_structure, probability, parent) tuple on success,
        or None if the substitution should be skipped.
        """
        # Filter out any structures with identical species
        if parent.has_species(diff_spec):
            return None

        # Ensure parent has as many species as target
        if len(parent.species) != len(species):
            return None

        # Get species to be substituted; ensure only 1 species is obtained
        extra = set(parent.get_spec_strs()) - set(map(unparse_spec, species)) - {diff_spec_str}
        if len(extra) > 1:
            return None
        (alt_spec,) = extra

        if parse_spec(alt_spec)[1] != diff_spec[1]:
            # Different charge
            return None

        try:
            p = diff_sub_probs.loc[alt_spec]
        except KeyError:
            # Not in the Series
            return None

        if thresh is not None and p <= thresh:
            return None

        try:
            self.cm._mutate_structure(parent, alt_spec, diff_spec_str)
        except ValueError:
            # Poorly decorated
            return None

        return (
            self.cm._mutate_structure(parent, alt_spec, diff_spec_str),
            p,
            parent,
        )

    def predict_structs(
        self,
        species: list[tuple[str, int]],
        thresh: float | None = 1e-3,
        include_same: bool | None = True,
    ) -> Generator[tuple[SmactStructure, float, SmactStructure], None, None]:
        """
        Predict structures for a combination of species.

        Args:
        ----
            species: A list of (element, charge). The constituent species
                of the target compound.
            thresh: The probability threshold, below which to discard
                predictions.
            include_same: Whether to include unmodified structures
                from the database, i.e. structures containing all the
                same species. Defaults to True.

        Yields:
        ------
            Potential structures, as tuples of (structure, probability, parent).

        """
        # For now, consider just structures with the same species, and unary substitutions.
        # This means we need only consider structures with a difference of 0 or 1 species.

        if include_same:
            for identical in self.db.get_with_species(species, self.table):
                yield (identical, 1.0, identical)

        sub_spec = itertools.combinations(species, len(species) - 1)
        sub_spec = list(map(list, sub_spec))

        potential_unary_parents: list[list[SmactStructure]] = [
            self.db.get_with_species(specs, self.table) for specs in sub_spec
        ]

        for spec_idx, parents in enumerate(potential_unary_parents):
            # Get missing ion
            # Ensure a different ion is obtained
            if len(set(species) - set(sub_spec[spec_idx])) < 1:
                continue
            (diff_spec,) = set(species) - set(sub_spec[spec_idx])
            diff_spec_str = unparse_spec(diff_spec)

            # Determine conditional substitution likelihoods
            diff_sub_probs = self.cm.cond_sub_probs(diff_spec_str)

            for parent in parents:
                result = self._try_substitution(parent, diff_spec, diff_spec_str, diff_sub_probs, species, thresh)
                if result is not None:
                    yield result

    def _try_nary_substitution(
        self,
        parent: SmactStructure,
        diff_species: list[tuple[str, int]],
        diff_spec_str: list[str],
        diff_sub_probs: list[pd.Series],
        species: list[tuple[str, int]],
        n_ary: int,
        thresh: float | None,
    ) -> tuple[SmactStructure, float, SmactStructure] | None:
        """Attempt an n-ary substitution on a parent structure.

        Returns a (mutated_structure, probability, parent) tuple on success,
        or None if the substitution should be skipped.
        """
        # Filter out structures where the parent already has all the diff species
        if all(parent.has_species(ds) for ds in diff_species):
            return None

        # Ensure parent has as many species as target
        if len(parent.species) != len(species):
            return None

        # Get species to be substituted; ensure n species are obtained
        alt_spec_set = set(parent.get_spec_strs()) - set(map(unparse_spec, species)) - set(diff_spec_str)
        if len(alt_spec_set) != n_ary:
            return None
        alt_spec = list(alt_spec_set)

        try:
            p = [diff_sub_probs[i].loc[alt_spec[i]] for i in range(n_ary)]
        except KeyError:
            # Not in the Series
            return None

        p_prod = float(np.prod(p))

        if thresh is not None and p_prod <= thresh:
            return None

        try:
            self.cm._nary_mutate_structure(parent, alt_spec, diff_spec_str)
        except ValueError:
            # Poorly decorated
            return None

        return (
            self.cm._nary_mutate_structure(parent, alt_spec, diff_spec_str),
            p_prod,
            parent,
        )

    def nary_predict_structs(
        self,
        species: list[tuple[str, int]],
        n_ary: int | None = 2,
        thresh: float | None = 1e-3,
        include_same: bool | None = True,
    ) -> Generator[tuple[SmactStructure, float, SmactStructure], None, None]:
        """
        Predicts structures for a combination of species.

        Args:
        ----
            species: A list of (element, charge). The constituent species
             of the target compound.
            thresh: The probability threshold, below which to discard predictions.
            n_ary: The number of species in a parent compound to replace.
            include_same: Whether to include unmodified structures from the database,
             i.e. structures containing all the same species.

        Yields:
        ------
            Potential structures, as tuples of (structure, probability, parent).

        """
        if include_same:
            for identical in self.db.get_with_species(species, self.table):
                yield (identical, 1.0, identical)

        # Ensure that we can obtain a subset of species of the target compound
        if n_ary is None:
            return
        if n_ary < 1 or n_ary >= len(species):
            return
        sub_species = itertools.combinations(species, len(species) - n_ary)
        sub_species = list(map(list, sub_species))

        potential_nary_parents: list[list[SmactStructure]] = [
            self.db.get_with_species(specs, self.table) for specs in sub_species
        ]

        for spec_idx, parents in enumerate(potential_nary_parents):
            # Get missing ions
            diff_species = list(set(species) - set(sub_species[spec_idx]))
            diff_spec_str = [unparse_spec(i) for i in diff_species]

            diff_sub_probs = [self.cm.cond_sub_probs(i) for i in diff_spec_str]

            for parent in parents:
                result = self._try_nary_substitution(
                    parent, diff_species, diff_spec_str, diff_sub_probs, species, n_ary, thresh
                )
                if result is not None:
                    yield result
