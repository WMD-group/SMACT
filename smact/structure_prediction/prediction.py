"""Structure prediction implementation."""

from typing import Generator, List, Tuple, Optional

from .database import StructureDB
from .mutation import CationMutator
from .structure import SmactStructure


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

    def __init__(self, mutator: CationMutator, struct_db: StructureDB):
        """Initialize class."""
        self.cm = mutator
        self.db = struct_db

    def predict_structs(
      self,
      species: List[Tuple[str, int]],
      thresh: Optional[float] = 1e-3, ) -> Generator[Tuple[SmactStructure, float], None, None]:
        """Predict structures for a combination of species.

        Args:
            species: A list of (element, charge). The constituent species
                of the target compound.
            thresh: The probability threshold, below which to discard
                predictions.

        Yields:
            Potential structures, as tuples of (structure, probability).

        """
        raise NotImplementedError
