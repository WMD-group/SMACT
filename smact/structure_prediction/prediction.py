"""Structure prediction implementation.

Todo:
    * Test with a fully populated database.
    * Implement n-ary substitution probabilities;
      at the moment, only zero- and single-species
      substitutions are considered.

"""

import itertools
from typing import Generator, List, Optional, Tuple

import numpy as np

from .database import StructureDB
from .mutation import CationMutator
from .structure import SmactStructure
from .utilities import parse_spec, unparse_spec


class StructurePredictor:
    """Provides structure prediction functionality.

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

    def __init__(
        self, mutator: CationMutator, struct_db: StructureDB, table: str
    ):
        """Initialize class.

        Args:
            mutator: A :class:`CationMutator` for probability calculations.
            struct_db: A :class:`StructureDB` from which to read strucutures
                to attempt to mutate.
            table: The table to reference within the database

        """
        self.cm = mutator
        self.db = struct_db
        self.table = table

    def predict_structs(
        self,
        species: List[Tuple[str, int]],
        thresh: Optional[float] = 1e-3,
        include_same: Optional[bool] = True,
    ) -> Generator[Tuple[SmactStructure, float, SmactStructure], None, None]:
        """Predict structures for a combination of species.

        Args:
            species: A list of (element, charge). The constituent species
                of the target compound.
            thresh: The probability threshold, below which to discard
                predictions.
            include_same: Whether to include unmodified structures
                from the database, i.e. structures containing all the
                same species. Defaults to True.

        Yields:
            Potential structures, as tuples of (structure, probability, parent).

        """
        # For now, consider just structures with the same species, and unary substitutions.
        # This means we need only consider structures with a difference of 0 or 1 species.

        if include_same:
            for identical in self.db.get_with_species(species, self.table):
                yield (identical, 1.0, identical)

        sub_spec = itertools.combinations(species, len(species) - 1)
        sub_spec = list(map(list, sub_spec))

        potential_unary_parents: List[List[SmactStructure]] = list(
            self.db.get_with_species(specs, self.table) for specs in sub_spec
        )

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
                # print("Testing parent")
                # Filter out any structures with identical species
                if parent.has_species(diff_spec):
                    continue

                # Ensure parent has as many species as target
                if len(parent.species) != len(species):
                    continue

                # Determine probability
                # Get species to be substituted
                # Ensure only 1 species is obtained
                if (
                    len(
                        set(parent.get_spec_strs())
                        - set(map(unparse_spec, species))
                        - {diff_spec_str}
                    )
                    > 1
                ):
                    continue
                (alt_spec,) = (
                    set(parent.get_spec_strs())
                    - set(map(unparse_spec, species))
                    - {diff_spec_str}
                )

                if parse_spec(alt_spec)[1] != diff_spec[1]:
                    # Different charge
                    continue

                try:
                    p = diff_sub_probs.loc[alt_spec]
                except:
                    # Not in the Series
                    continue

                if p > thresh:
                    try:
                        mutated = self.cm._mutate_structure(
                            parent, alt_spec, diff_spec_str
                        )
                    except ValueError:
                        # Poorly decorated
                        continue
                    yield (
                        self.cm._mutate_structure(
                            parent, alt_spec, diff_spec_str
                        ),
                        p,
                        parent,
                    )

    def nary_predict_structs(
        self,
        species: List[Tuple[str, int]],
        n_ary: Optional[int] = 2,
        thresh: Optional[float] = 1e-3,
        include_same: Optional[bool] = True,
    ) -> Generator[Tuple[SmactStructure, float, SmactStructure], None, None]:
        """Predicts structures for a combination of species.

        Args:
            species: A list of (element, charge). The constituent species
             of the target compound.
            thresh: The probability threshold, below which to discard predictions.
            n_ary: The number of species in a parent compound to replace.
            include_same: Whether to include unmodified structures from the database,
             i.e. structures containing all the same species.

        Yields:
            Potential structures, as tuples of (structure, probability, parent).
        """

        if include_same:
            for identical in self.db.get_with_species(species, self.table):
                yield (identical, 1.0, identical)

        # Ensure that we can obtain a subset of species of the target compound
        if len(species) - n_ary == 0:
            return None
        sub_species = itertools.combinations(species, len(species) - n_ary)
        sub_species = list(map(list, sub_species))

        potential_nary_parents: List[List[SmactStructure]] = list(
            self.db.get_with_species(specs, self.table)
            for specs in sub_species
        )

        for spec_idx, parents in enumerate(potential_nary_parents):
            # Get missing ions
            # Ensure we get the correct number of species
            # Ensure the ions obtained are different
            # if len(set(species) - set(sub_species[spec_idx])) !=2:
            # continue
            diff_species = list(set(species) - set(sub_species[spec_idx]))
            diff_spec_str = [unparse_spec(i) for i in diff_species]

            diff_sub_probs = [self.cm.cond_sub_probs(i) for i in diff_spec_str]

        for parent in parents:
            # print("testing parent")
            # Filter out any structures with identical species
            if n_ary == 1:
                if parent.has_species(diff_species[0]):
                    continue
            elif n_ary == 2:
                if parent.has_species(diff_species[0]) and parent.has_species(
                    diff_species[1]
                ):
                    continue
            elif n_ary == 3:
                if (
                    parent.has_species(diff_species[0])
                    and parent.has_species(diff_species[1])
                    and parent.has_species(diff_species[2])
                ):
                    continue

            # Ensure parent has as many species as target
            if len(parent.species) != len(species):
                continue

            # Determine probability
            # Get species to be substituted
            # Ensure n species are obtained

            if (
                len(
                    set(parent.get_spec_strs())
                    - set(map(unparse_spec, species))
                    - set(diff_species)
                )
                != n_ary
            ):
                continue
            alt_spec = list(
                set(parent.get_spec_strs())
                - set(map(unparse_spec, species))
                - set(diff_species)
            )

            # Need to consider p(A,X)p(B,Y) and p(A,Y)p(B,X)
            # if utilities.parse_spec(alt_spec_1)[1] != diff_species_1[1] and utilities.parse_spec(alt_spec_2)[1] != diff_species_2[1] :
            # Different charge
            # continue

            try:
                p = []
                for i in range(n_ary):
                    p.append(diff_sub_probs[i].loc[alt_spec[i]])
            except:
                # Not in the Series
                continue

            p = np.prod(p)

            if p > thresh:
                try:
                    mutated = self.cm._nary_mutate_structure(
                        parent, alt_spec, diff_spec_str
                    )

                except ValueError:
                    # Poorly decorated
                    continue
                yield (
                    self.cm._nary_mutate_structure(
                        parent, alt_spec, diff_spec_str
                    ),
                    p,
                    parent,
                )
